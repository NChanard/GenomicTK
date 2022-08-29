#' IndexFeatures
#'
#' Bin multiple GRanges and summary all in one GRange. Could overlap ranges with constraints regions
#' @param gRange.gnr_lst <GRanges or GRangesList or list[GRanges]>: some GRanges or a list of GRanges or a GRangesList
#' @param constraint.gnr <GRanges>: A GRange of constraint regions. If NULL chromosomes in chromSize.dtf are used (Default NULL)
#' @param chromSize.dtf <data.frame>: A data.frame of genome where first colum correspond to the chromosomes names, and the second column correspond to the chromosomes lengths in base pairs.
#' @param binSize.int <integer>: A number that specify the width bins.
#' @param method.chr <character>: A string of a summary method name as 'mean', 'median', 'sum', 'max, 'min'. (Default 'mean'')
#' @param variablesName.chr_vec <character> : A character vector that specify the metadata columns of GRanges on which apply the summary method.
#' @param cores.num <integer> : An integer to specify the number of cores. (Default 1)
#' @param verbose.bln <logical>: A logical value. If TRUE show the progression in console. (Default TRUE)
#' @return A GRange object.
#' @examples
#' GRanges_1.gr <- GenomicRanges::GRanges(
#'     seqnames = S4Vectors::Rle( c("chr1", "chr2"), c(4, 2) ),
#'     ranges = IRanges::IRanges(
#'         start = c(5, 220, 260, 540, 50, 540),
#'         end =  c(80, 240, 280, 560, 150, 580)),
#'     strand = S4Vectors::Rle(BiocGenerics::strand(c("-","+","*","-","+","+")), c(1,1,1,1,1,1)),
#'     score = 1:6,
#'     name = letters[1:6]
#' )
#' GRanges_1.gr
#' GRanges_2.gr <- GenomicRanges::GRanges(
#'     seqnames = S4Vectors::Rle( c("chr1", "chr2"), c(3, 2) ),
#'     ranges = IRanges::IRanges(
#'         start = c(140, 220, 450, 320, 360),
#'         end =  c(160, 280, 550, 350, 380)),
#'     strand = S4Vectors::Rle(BiocGenerics::strand(c("-","+","*")), c(1,2,2)),
#'     score = seq(10,50,10),
#'     name = letters[7:11]
#' )
#' GRanges_2.gr
#' GRanges_Constraint <- GenomicRanges::GRanges(
#'     seqnames = S4Vectors::Rle( c("chr1", "chr2"), c(2, 2) ),
#'     ranges = IRanges::IRanges(
#'         start = c(1, 251, 1, 50),
#'         end =  c(250, 490, 400, 190),
#'         names = paste0("C_",1:4)),
#'     strand = S4Vectors::Rle(BiocGenerics::strand(c("*")), c(4))
#' )
#' GRanges_Constraint
#' index.tbl = IndexFeatures(
#'     gRange.gnr_lst=list(GR_1 = GRanges_1.gr, GR_2 = GRanges_2.gr),
#'     constraint.gnr=GRanges_Constraint,
#'     chromSize.dtf=data.frame(
#'         seqnames=c("chr1", "chr2"),
#'         seqlengths=c(3500,2000)),
#'     binSize.int=100,
#'     method.chr ="mean",
#'     variablesName.chr_vec = c("score")
#'     )
#' index.tbl

IndexFeatures <- function(gRange.gnr_lst=NULL, constraint.gnr=NULL, chromSize.dtf=NULL, binSize.int=NULL, method.chr="mean", variablesName.chr_vec=NULL,cores.num=1, verbose.bln=TRUE){
    # Constraint Informations
        if (is.null(constraint.gnr)){
            constraint.gnr <- GenomicRanges::GRanges(
                seqnames = S4Vectors::Rle(chromSize.dtf[, 1]),
                ranges = IRanges::IRanges(start = rep(1, length(chromSize.dtf[, 2])), end = chromSize.dtf[, 2]),
                strand = S4Vectors::Rle(BiocGenerics::strand('*')),
                name =  chromSize.dtf[, 1]
            )
        }else{
            if(is.null(constraint.gnr$name) | length(which(!is.na(constraint.gnr$name)))==0 ){
                constraint.gnr$name = paste0("Constraint_", seq_along(constraint.gnr))
            }
        }
        seqLevelsStyle.chr <- GenomeInfoDb::seqlevelsStyle(constraint.gnr)
        if(length(seqLevelsStyle.chr)>1){
            seqLevelsStyle.chr <- seqLevelsStyle.chr[[1]]
            GenomeInfoDb::seqlevelsStyle(constraint.gnr) <- seqLevelsStyle.chr
        }
        binnedConstraint.gnr <- GenomicTK::BinGRanges(gRange.gnr=constraint.gnr, chromSize.dtf=chromSize.dtf, binSize.int=binSize.int, verbose.bln=verbose.bln, reduce.bln=FALSE, cores.num=cores.num)
    # Feature Names
        if (inherits(gRange.gnr_lst,"GRanges")){
            gRange.gnr_lst %<>% list(.) %>% magrittr::set_names(., "Feature")
        }else if (inherits(gRange.gnr_lst,"GRangesList")){
            gRange.gnr_lst  %<>% as.list(.)
        }
        if (gRange.gnr_lst %>% names %>% is.null){
            gRange.gnr_lst %<>% length %>% seq(1, .) %>% lapply(., function(i){paste0("Feature_", i)}) %>% unlist %>% magrittr::set_names(gRange.gnr_lst, .)
        }
        gRange.gnr_lst %<>% lapply(., length) %>% unlist %>% order(., decreasing=TRUE) %>% magrittr::extract(gRange.gnr_lst, .)
        feature.chr_vec <- names(gRange.gnr_lst)
    # GRanges Binning
        jobLenght.num <- length(gRange.gnr_lst)
        start.tim <- Sys.time()
        binnedFeature.lst <- lapply(seq_len(jobLenght.num), function(feature.ndx){
            feature.chr <- feature.chr_vec[[feature.ndx]]
            feature.gnr <- gRange.gnr_lst[[feature.chr ]] %>% IRanges::subsetByOverlaps(.,constraint.gnr)
            GenomeInfoDb::seqlevelsStyle(feature.gnr) <- seqLevelsStyle.chr
            binnedFeature.gnr <- GenomicTK::BinGRanges(gRange.gnr=feature.gnr, chromSize.dtf=chromSize.dtf, binSize.int=binSize.int,  method.chr=method.chr, variablesName.chr_vec=variablesName.chr_vec, verbose.bln=verbose.bln, reduce.bln=TRUE, cores.num=cores.num)
            binnedFeat.tbl = tibble::tibble(BinnedFeature.ndx = seq_along(binnedFeature.gnr),Feature.name = binnedFeature.gnr$name) %>%
                tidyr::unnest(.,cols = c(Feature.name)) %>%
                dplyr::group_by(Feature.name) %>%
                tidyr::nest(.) %>%
                magrittr::set_names(c("Feature.name","BinnedFeature.ndx"))
            binnedConstraint.tbl = tibble::tibble(BinnedConstraint.ndx = seq_along(binnedConstraint.gnr),Constraint.name = binnedConstraint.gnr$name) %>%
                dplyr::group_by(Constraint.name) %>%
                tidyr::nest(.) %>%
                magrittr::set_names(c("Constraint.name","BinnedConstraint.ndx"))
            featConstOvlp.ovlp <- GenomicRanges::findOverlaps(feature.gnr, constraint.gnr)
            featConstOvlp.tbl <- tibble::tibble(Feature.name = feature.gnr$name[featConstOvlp.ovlp@from],Constraint.name = constraint.gnr$name[featConstOvlp.ovlp@to]) %>%
                dplyr::left_join(., binnedFeat.tbl, by="Feature.name") %>%
                dplyr::select(-Feature.name) %>%
                tidyr::unnest(.,cols = c(BinnedFeature.ndx)) %>%
                unique %>%
                dplyr::group_by(Constraint.name) %>%
                tidyr::nest(.) %>%
                magrittr::set_names(c("Constraint.name","BinnedFeature.ndx")) %>%
                dplyr::left_join(.,binnedConstraint.tbl, by="Constraint.name")
            subJobLenght.num <- featConstOvlp.tbl %>% nrow
            start.tim <- Sys.time()
            if(cores.num==1){
                binnedFeature.gnr_lst <- lapply(seq_len(subJobLenght.num),function(row.ndx){
                    if(verbose.bln){DevTK::ShowLoading(start.tim, row.ndx+(feature.ndx-1)*subJobLenght.num,(subJobLenght.num*jobLenght.num))}
                    ranges.ndx <- featConstOvlp.tbl$BinnedFeature.ndx[row.ndx] %>% unlist(.,use.names=FALSE)
                    constraint.ndx <- featConstOvlp.tbl$BinnedConstraint.ndx[row.ndx] %>% unlist(.,use.names=FALSE)
                    subBinnedFeature.gnr <- IRanges::subsetByOverlaps(binnedFeature.gnr[ranges.ndx], binnedConstraint.gnr[constraint.ndx])
                    subBinnedFeature.gnr$constraint <- featConstOvlp.tbl$Constraint.name[row.ndx]
                    subBinnedFeature.gnr %>% return(.)
                })
            }else if(cores.num>=2){
                parCl <- parallel::makeCluster(cores.num, type ="FORK")
                doParallel::registerDoParallel(parCl)
                binnedFeature.gnr_lst <- parallel::parLapply(parCl,seq_len(subJobLenght.num),function(row.ndx){
                    ranges.ndx <- featConstOvlp.tbl$BinnedFeature.ndx[row.ndx] %>% unlist(.,use.names=FALSE)
                    constraint.ndx <- featConstOvlp.tbl$BinnedConstraint.ndx[row.ndx] %>% unlist(.,use.names=FALSE)
                    subBinnedFeature.gnr <- IRanges::subsetByOverlaps(binnedFeature.gnr[ranges.ndx], binnedConstraint.gnr[constraint.ndx])
                    subBinnedFeature.gnr$constraint <- featConstOvlp.tbl$Constraint.name[row.ndx]
                    subBinnedFeature.gnr %>% return(.)
                })
                parallel::stopCluster(parCl)
                DevTK::KillZombies()
            }
            binnedFeature.gnr <- GenomicTK::MergeGRanges(binnedFeature.gnr_lst, sort.bln=FALSE, reduce.bln=FALSE)
            binnedFeature.gnr$bln <- 1
            names(S4Vectors::mcols(binnedFeature.gnr)) %<>% paste0(feature.chr, ".", .)
            names(S4Vectors::mcols(binnedFeature.gnr))[which(names(S4Vectors::mcols(binnedFeature.gnr)) == paste0(feature.chr, ".bin"))] <- "bin"
            names(S4Vectors::mcols(binnedFeature.gnr))[which(names(S4Vectors::mcols(binnedFeature.gnr)) == paste0(feature.chr, ".constraint"))] <- "constraint"
            binnedFeature.gnr$name <- paste0(binnedFeature.gnr$bin, ":", binnedFeature.gnr$constraint)
            metadataBinnedFeature.dtf <- binnedFeature.gnr %>% S4Vectors::mcols(.) %>% data.frame
            S4Vectors::mcols(binnedFeature.gnr) <- NULL
            return(list(binnedFeature.gnr=binnedFeature.gnr, featureMetadata.dtf=metadataBinnedFeature.dtf))
        })

        binnedIndex.gnr <- binnedFeature.lst %>% lapply(., "[[", "binnedFeature.gnr") %>% GenomicTK::MergeGRanges(., sort.bln=FALSE, reduce.bln=FALSE)
        S4Vectors::mcols(binnedIndex.gnr) <- binnedFeature.lst %>% lapply(., "[[", "featureMetadata.dtf") %>% plyr::rbind.fill(.)
        ids.lst <- binnedIndex.gnr$name
        dupplicatedIds.lst <- ids.lst %>% duplicated %>% which %>% magrittr::extract(ids.lst, .) %>% unique
        idDuplicated.ndx <- ids.lst %>% magrittr::is_in(., dupplicatedIds.lst) %>% which
    # Merge GRanges into one index and duplicated Bin handle
        if(length(idDuplicated.ndx)){
            binnedIndexDuplicated.tbl <- binnedIndex.gnr[idDuplicated.ndx] %>% data.frame(.) %>% tibble::tibble(.) %>% dplyr::group_by(name) %>% tidyr::nest(.)
            binnedIndexNoDuplicated.tbl <- binnedIndex.gnr[-idDuplicated.ndx] %>% data.frame(.) %>% tibble::tibble(.)
            start.tim <- Sys.time()
            jobLenght.num <- nrow(binnedIndexDuplicated.tbl)
            if(cores.num==1){
                binnedIndexDuplicated.lst <- lapply(seq_len(jobLenght.num), function(row.ndx){
                    if(verbose.bln){DevTK::ShowLoading(start.tim, row.ndx, jobLenght.num)}
                    rowName.chr <- binnedIndexDuplicated.tbl$name[[row.ndx]]
                    row =  binnedIndexDuplicated.tbl$data[[row.ndx]]
                    lapply(seq_along(row),function(col.ndx){
                        col=dplyr::pull(row,col.ndx)
                        # col.name = names(row)[col.ndx]
                            if(length(unique(stats::na.omit(col)))==0){
                                return(NA)
                            }else if(length(unique(stats::na.omit(col)))==1) {
                                return(unique(stats::na.omit(col)))
                            }else {
                                return(list(unlist(col)))
                            }
                        })  %>% magrittr::set_names(.,names(row)) %>% tibble::as_tibble(.) %>% tibble::add_column(name = rowName.chr) %>% return(.)
                })
            }else if(cores.num>=2){
                parCl <- parallel::makeCluster(cores.num, type ="FORK")
                doParallel::registerDoParallel(parCl)
                binnedIndexDuplicated.lst <- parallel::parLapply(parCl,seq_len(jobLenght.num), function(row.ndx){
                    rowName.chr <- binnedIndexDuplicated.tbl$name[[row.ndx]]
                    row =  binnedIndexDuplicated.tbl$data[[row.ndx]]
                    lapply(seq_along(row),function(col.ndx){
                        col=dplyr::pull(row,col.ndx)
                        # col.name = names(row)[col.ndx]
                            if(length(unique(stats::na.omit(col)))==0){
                                return(NA)
                            }else if(length(unique(stats::na.omit(col)))==1) {
                                return(unique(stats::na.omit(col)))
                            }else {
                                return(list(unlist(col)))
                            }
                        })  %>% magrittr::set_names(.,names(row)) %>% tibble::as_tibble(.) %>% tibble::add_column(name = rowName.chr) %>% return(.)
                })
                parallel::stopCluster(parCl)
                DevTK::KillZombies()
            }
            binnedIndexDuplicated.tbl <- dplyr::bind_rows(binnedIndexDuplicated.lst)
            binnedIndex.gnr <- rbind(binnedIndexDuplicated.tbl, binnedIndexNoDuplicated.tbl) %>% data.frame(.) %>% methods::as(., "GRanges")
        }
        for(featureName.chr in feature.chr_vec){
            colname.chr =  paste0(featureName.chr, ".bln")
            S4Vectors::mcols(binnedIndex.gnr)[which(is.na(S4Vectors::mcols(binnedIndex.gnr)[, colname.chr])),colname.chr] <- 0
            S4Vectors::mcols(binnedIndex.gnr)[,colname.chr] <- methods::as(S4Vectors::mcols(binnedIndex.gnr)[, colname.chr], "Rle")
        }
        S4Vectors::mcols(binnedIndex.gnr) %<>% as.data.frame %>% dplyr::select(name,bin,constraint,tidyselect::everything()) 
    return(sort(binnedIndex.gnr))
}