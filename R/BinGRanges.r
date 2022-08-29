#' BinGRanges
#' 
#' Bin a GRanges and allow to apply a summary method ('mean', 'median', 'sum', 'max, 'min' ...) to chossen numerical variables of ranges in a same bin.
#' @param gRange.gnr <GRanges>: A GRange to bin
#' @param chromSize.dtf <data.frame>: A data.frame of genome where first colum correspond to the chromosomes names, and the second column correspond to the chromosomes lengths in base pairs.
#' @param binSize.int <integer>: A number that specify the width bins.
#' @param method.chr <character>: A string of a summary method name as 'mean', 'median', 'sum', 'max, 'min'. (Default 'mean'')
#' @param variablesName.chr_vec <character> : A character vector that specify the metadata columns of GRanges on which apply the summary method.
#' @param na.rm <logical> : A logical value indicating whether 'NA' values should be stripped before the computation proceeds. (Default TRUE)
#' @param cores.num <integer> : An integer to specify the number of cores. (Default 1)
#' @param reduce.bln <logical> : A logical value indicating whether duplicated Bin must been reduced with de summary method. (Default TRUE)
#' @param verbose.bln <logical>: A logical value. If TRUE show the progression in console. (Default TRUE)
#' @return A GRange object
#' @examples
#' GRange.gnr <- GenomicRanges::GRanges(
#'     seqnames = S4Vectors::Rle(c("chr1", "chr2"), c(3, 1)),
#'     ranges = IRanges::IRanges(c(1,201,251,1), end = c(200,250,330,100), names = letters[1:4]),
#'     strand = S4Vectors::Rle(BiocGenerics::strand(c("*")), 4),
#'     score = c(50,NA,100,30)
#'     )
#' GRange.gnr
#' chromSize.dtf = data.frame(c("chr1","chr2"),c(350,100))
#' binSize.int <- 100
#' binnedGRanges.gnr <- BinGRanges(
#'     gRange.gnr = GRange.gnr,
#'     chromSize.dtf=chromSize.dtf,
#'     binSize.int=binSize.int,
#'     method.chr ="mean",
#'     variablesName.chr_vec="score",
#'     na.rm=TRUE
#' )
#' binnedGRanges.gnr
BinGRanges = function (gRange.gnr=NULL, chromSize.dtf=NULL, binSize.int=NULL, method.chr="mean", variablesName.chr_vec=NULL, na.rm=TRUE, cores.num=1, reduce.bln=TRUE, verbose.bln=TRUE){
        if(is.null(chromSize.dtf)){
            seqlengths.lst <- GenomeInfoDb::seqlengths(gRange.gnr)
        }else{
            seqlengths.lst <- chromSize.dtf %>% dplyr::pull(2) %>% magrittr::set_names({chromSize.dtf%>% dplyr::pull(1)})
            GenomeInfoDb::seqlengths(gRange.gnr) <- seqlengths.lst
        }
        binnedGenome.gnr <- GenomicRanges::tileGenome(seqlengths.lst, tilewidth=binSize.int, cut.last.tile.in.chrom=TRUE)
        ovlp.dtf <- GenomicRanges::findOverlaps(binnedGenome.gnr , gRange.gnr)
        binnedGRanges.gnr <- binnedGenome.gnr[ovlp.dtf@from]
        S4Vectors::mcols(binnedGRanges.gnr) <- S4Vectors::mcols(gRange.gnr[ovlp.dtf@to])
        binnedGRanges.gnr$bin <-  paste0(GenomeInfoDb::seqnames(binnedGRanges.gnr),":", ceiling(BiocGenerics::start(binnedGRanges.gnr)/binSize.int))
        dupplicated.id <-  binnedGRanges.gnr$bin %>% duplicated() %>% which(binnedGRanges.gnr$bin %in% .) %>% binnedGRanges.gnr$bin[.]
        if(reduce.bln && length(dupplicated.id)){
            binnedGRange.tbl <- tibble::tibble(data.frame(binnedGRanges.gnr))
            nodup_binnedGRange.tbl <- binnedGRange.tbl %>% dplyr::slice(., which(DevTK::NotIn(binnedGRange.tbl$bin, dupplicated.id)) )
            dup_binnedGRange.tbl <- binnedGRange.tbl %>% dplyr::slice(., which(binnedGRange.tbl$bin %in% dupplicated.id) ) %>% dplyr::group_by(bin) %>% tidyr::nest()
            jobLenght.num <- nrow(dup_binnedGRange.tbl)
            if(cores.num==1){
                start.tim <- Sys.time()
                if(verbose.bln){cat("\n")}
                dup_binnedGRange.tbl <- lapply(seq_len(jobLenght.num),function(row.ndx){
                    if(verbose.bln){DevTK::ShowLoading(start.tim,row.ndx, jobLenght.num)}
                    rowName.chr <- dup_binnedGRange.tbl$bin[[row.ndx]]
                    row <- dup_binnedGRange.tbl$data[[row.ndx]]
                    lapply(seq_along(row),function(col.ndx){
                        col <- dplyr::pull(row,col.ndx)
                        colName.chr <- names(row)[col.ndx]
                        if(is.numeric(col) & colName.chr %in% variablesName.chr_vec){
                            return(as.numeric(eval(parse(text=method.chr))(as.numeric(col),na.rm=na.rm)))
                        }else if(length(unique(col))==1){
                            return(unique(col))
                        }else{
                            return(list(col))
                        }
                    }) %>% magrittr::set_names(.,names(row)) %>% tibble::as_tibble(.) %>% tibble::add_column(bin =  rowName.chr) %>% return(.)
                }) %>% plyr::rbind.fill(.) %>% tibble::tibble(.) %>% dplyr::mutate(strand=forcats::as_factor(strand))
                if(verbose.bln){cat("\n")}
            }else if(cores.num>=2){
                parCl <- parallel::makeCluster(cores.num, type ="FORK")
                doParallel::registerDoParallel(parCl)
                dup_binnedGRange.tbl <- parallel::parLapply(parCl,seq_len(jobLenght.num),function(row.ndx){
                    rowName.chr <- dup_binnedGRange.tbl$bin[[row.ndx]]
                    row <- dup_binnedGRange.tbl$data[[row.ndx]]
                    lapply(seq_along(row),function(col.ndx){
                        col <- dplyr::pull(row,col.ndx)
                        colName.chr <- names(row)[col.ndx]
                        if(is.numeric(col) & colName.chr %in% variablesName.chr_vec){
                            return(as.numeric(eval(parse(text=method.chr))(as.numeric(col),na.rm=na.rm)))
                        }else if(length(unique(col))==1){
                            return(unique(col))
                        }else{
                            return(list(col))
                        }
                    }) %>% magrittr::set_names(.,names(row)) %>% tibble::as_tibble(.) %>% tibble::add_column(bin =  rowName.chr) %>% return(.)
                }) %>% plyr::rbind.fill(.) %>% tibble::tibble(.) %>% dplyr::mutate(strand=forcats::as_factor(strand))
                parallel::stopCluster(parCl)
                DevTK::KillZombies()
            }
            for(colName.chr in names(dup_binnedGRange.tbl)){
                method.chr <- dup_binnedGRange.tbl %>% dplyr::pull(.,dplyr::all_of(colName.chr)) %>% class %>% paste0("as.",.) %>% parse(text=.) %>% eval(.)
                nodup_binnedGRange.tbl %<>% dplyr::mutate(dplyr::across(dplyr::all_of(colName.chr),method.chr))
            }
            binnedGRange.tbl <- dplyr::bind_rows(dup_binnedGRange.tbl,nodup_binnedGRange.tbl)
            binnedGRanges.gnr <- methods::as(binnedGRange.tbl,'GRanges')
        }
        binnedGRanges.gnr %<>% sort
        GenomeInfoDb::seqinfo(binnedGRanges.gnr) <- GenomeInfoDb::seqinfo(gRange.gnr)
        return(binnedGRanges.gnr)
}