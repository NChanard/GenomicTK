#' Bin a GRanges.
#'
#' BinGRanges
#' Bin a GRanges and allow to apply a summary method (e.g: 'mean', 'median', 'sum', 'max, 'min' ...) to a chossen numericals variables of ranges in a same bin.
#' @param gRange.gnr <GRanges>: a GRanges to bin.
#' @param chromSize.dtf <data.frame>: a data.frame where first colum correspond to the chromosomes names, and the second column correspond to the chromosomes lengths in base pairs.
#' @param binSize.num <numerical>: width of the bins.
#' @param method.chr <character>: name of a summary method as 'mean', 'median', 'sum', 'max, 'min'. (Default 'mean')
#' @param variablesName.chr_vec <character> : a character vector that specify the metadata columns of GRanges on which apply the summary method.
#' @param na.rm <logical> : a logical value indicating whether 'NA' values should be stripped before the computation proceeds. (Default TRUE)
#' @param cores.num <numerical> : the number of cores. (Default 1)
#' @param reduce.bln <logical> : whether duplicated Bin must been reduced with de summary method. (Default TRUE)
#' @param verbose.bln <logical>: if TRUE show the progression in console. (Default TRUE)
#' @return A binned GRanges.
#' @examples
#' GRange.gnr <- GenomicRanges::GRanges(
#'     seqnames = S4Vectors::Rle(c("chr1", "chr2"), c(3, 1)),
#'     ranges = IRanges::IRanges(c(1,201,251,1), end = c(200,250,330,100), names = letters[1:4]),
#'     strand = S4Vectors::Rle(BiocGenerics::strand(c("*")), 4),
#'     score = c(50,NA,100,30)
#'     )
#' GRange.gnr
#' chromSize.dtf = data.frame(c("chr1","chr2"),c(350,100))
#' binSize.num <- 100
#' binnedGRanges.gnr <- BinGRanges(
#'     gRange.gnr = GRange.gnr,
#'     chromSize.dtf=chromSize.dtf,
#'     binSize.num=binSize.num,
#'     method.chr ="mean",
#'     variablesName.chr_vec="score",
#'     na.rm=TRUE
#' )
#' binnedGRanges.gnr
BinGRanges = function (gRange.gnr=NULL, chromSize.dtf=NULL, binSize.num=NULL, method.chr="mean", variablesName.chr_vec=NULL, na.rm=TRUE, cores.num=1, reduce.bln=TRUE, verbose.bln=TRUE){
    if(is.null(chromSize.dtf)){
        seqlengths.lst <- GenomeInfoDb::seqlengths(gRange.gnr)
    }else{
        seqlengths.lst <- dplyr::pull(chromSize.dtf,2) |>
            stats::setNames(dplyr::pull(chromSize.dtf,1))
        seqlengths.lst <- seqlengths.lst[intersect(names(seqlengths.lst), levels(GenomeInfoDb::seqnames(gRange.gnr)@values))]
        gRange.gnr <- GenomeInfoDb::keepSeqlevels(gRange.gnr,value=names(seqlengths.lst), "coarse")
        GenomeInfoDb::seqlengths(gRange.gnr) <- seqlengths.lst
    }
    binnedGenome.gnr <- GenomicRanges::tileGenome(seqlengths.lst, tilewidth=binSize.num, cut.last.tile.in.chrom=TRUE)
    ovlp.dtf <- GenomicRanges::findOverlaps(binnedGenome.gnr, gRange.gnr)
    binnedGRanges.gnr <- binnedGenome.gnr[ovlp.dtf@from]
    S4Vectors::mcols(binnedGRanges.gnr) <- S4Vectors::mcols(gRange.gnr[ovlp.dtf@to])
    binnedGRanges.gnr$bin <-  paste0(GenomeInfoDb::seqnames(binnedGRanges.gnr),":", ceiling(BiocGenerics::start(binnedGRanges.gnr)/binSize.num))
    dupplicated.lgk <- duplicated(binnedGRanges.gnr$bin)
    dupplicated.id  <- binnedGRanges.gnr$bin[dupplicated.lgk]
    if(reduce.bln && length(dupplicated.id)){
        binnedGRange.tbl <- tibble::tibble(data.frame(binnedGRanges.gnr))
        nodup_binnedGRange.tbl <- dplyr::slice(binnedGRange.tbl, which(DevTK::NotIn(binnedGRange.tbl$bin, dupplicated.id)) )
        dup_binnedGRange.tbl   <- dplyr::slice(binnedGRange.tbl, which(binnedGRange.tbl$bin %in% dupplicated.id) ) 
        dup_binnedGRange.tbl   <- dplyr::group_by(dup_binnedGRange.tbl, bin=dup_binnedGRange.tbl$bin) |>
            tidyr::nest()
        jobLenght.num <- nrow(dup_binnedGRange.tbl)
        if(cores.num==1){
            start.tim <- Sys.time()
            if(verbose.bln){cat("\n")}
            dup_binnedGRange.lst <- lapply(seq_len(jobLenght.num),function(row.ndx){
                if(verbose.bln){DevTK::ShowLoading(start.tim,row.ndx, jobLenght.num)}
                rowName.chr <- dup_binnedGRange.tbl$bin[[row.ndx]]
                row <- dup_binnedGRange.tbl$data[[row.ndx]]
                row.lst <- lapply(seq_along(row),function(col.ndx){
                    col <- dplyr::pull(row,col.ndx)
                    colName.chr <- names(row)[col.ndx]
                    if(is.numeric(col) & colName.chr %in% variablesName.chr_vec){
                        return(as.numeric(eval(parse(text=method.chr))(as.numeric(col),na.rm=na.rm)))
                    }else if(length(unique(col))==1){
                        return(unique(col))
                    }else{
                        return(list(col))
                    }
                }) |>
                    stats::setNames(names(row))
                row.tbl <- tibble::as_tibble(row.lst) |>
                    tibble::add_column(bin = rowName.chr)
                return(row.tbl)
            })
            dup_binnedGRange.dtf <- plyr::rbind.fill(dup_binnedGRange.lst)
            dup_binnedGRange.tbl <- tibble::tibble(dup_binnedGRange.dtf)
            dup_binnedGRange.tbl <- dplyr::mutate(dup_binnedGRange.tbl, strand=forcats::as_factor(dup_binnedGRange.tbl$strand))
            if(verbose.bln){cat("\n")}
        }else if(cores.num>=2){
            parCl <- parallel::makeCluster(cores.num, type ="PSOCK")
            dup_binnedGRange.lst <- parallel::parLapply(parCl,seq_len(jobLenght.num),function(row.ndx){
                rowName.chr <- dup_binnedGRange.tbl$bin[[row.ndx]]
                row <- dup_binnedGRange.tbl$data[[row.ndx]]
                row.lst <- lapply(seq_along(row),function(col.ndx){
                    col <- dplyr::pull(row,col.ndx)
                    colName.chr <- names(row)[col.ndx]
                    if(is.numeric(col) & colName.chr %in% variablesName.chr_vec){
                        return(as.numeric(eval(parse(text=method.chr))(as.numeric(col),na.rm=na.rm)))
                    }else if(length(unique(col))==1){
                        return(unique(col))
                    }else{
                        return(list(col))
                    }
                }) |>
                    stats::setNames(names(row))
                row.tbl <- tibble::as_tibble(row.lst) |>
                    tibble::add_column(bin =  rowName.chr)
                return(row.tbl)
            }) 
            dup_binnedGRange.dtf <- plyr::rbind.fill(dup_binnedGRange.lst)
            dup_binnedGRange.tbl <- tibble::tibble(dup_binnedGRange.dtf)
            dup_binnedGRange.tbl <- dplyr::mutate(dup_binnedGRange.tbl, strand=forcats::as_factor(dup_binnedGRange.tbl$strand))
            parallel::stopCluster(parCl)
        }
        for(colName.chr in names(dup_binnedGRange.tbl)){
            method.chr <- dplyr::pull(dup_binnedGRange.tbl,dplyr::all_of(colName.chr)) |> class()
            method.fun <- eval(parse(text=paste0("as.",method.chr)))
            nodup_binnedGRange.tbl <- nodup_binnedGRange.tbl |> dplyr::mutate(dplyr::across(dplyr::all_of(colName.chr),method.fun))
        }
        binnedGRange.tbl <- dplyr::bind_rows(dup_binnedGRange.tbl,nodup_binnedGRange.tbl)
        binnedGRanges.gnr <- methods::as(binnedGRange.tbl,'GRanges')
    }
    binnedGRanges.gnr <- sort(binnedGRanges.gnr)
    GenomeInfoDb::seqinfo(binnedGRanges.gnr) <- GenomeInfoDb::seqinfo(gRange.gnr)
    return(binnedGRanges.gnr)
}