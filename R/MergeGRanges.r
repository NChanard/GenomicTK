#' Merge GRanges or a list of GRanges
#' 
#' MergeGRanges
#' Merge GRanges or a list of GRanges
#' @param ... <GRanges or GRangesList or list[GRanges]>: some GRanges or a list of GRanges or a GRangesList.
#' @param sort.bln <logical>: whether the result should be sorted. (Default FALSE)
#' @param reduce.bln <logical>: whether the result should be reduce. See GenomicRanges::reduce for more details. (Default FALSE)
#' @return a GRange object.
#' @examples
#' GRange_1.grn <- GenomicRanges::GRanges(
#'     seqnames = S4Vectors::Rle(c("chr1", "chr2", "chr1"), c(1, 3, 1)),
#'     ranges = IRanges::IRanges(101:105, end = 111:115, names = letters[1:5]),
#'     strand = S4Vectors::Rle(BiocGenerics::strand(c("-", "+", "*", "+")), c(1, 1, 2, 1)),
#'     score = 1:5
#' )
#' GRange_2.grn <- GenomicRanges::GRanges(
#'     seqnames = S4Vectors::Rle(c("chr1", "chr3"), c(1, 4)),
#'     ranges = IRanges::IRanges(106:110, end = 116:120, names = letters[6:10]),
#'     strand = S4Vectors::Rle(BiocGenerics::strand(c("*", "+", "-")), c(2, 1, 2)),
#'     score = 6:10
#' )
#' GRange_1.grn
#' GRange_2.grn
#' MergeGRanges(GRange_1.grn,GRange_2.grn)
#' GRange.lst = list(GRange_1.grn,GRange_2.grn)
#' MergeGRanges(GRange.lst)
#' MergeGRanges(GRange.lst, reduce.bln=TRUE)
MergeGRanges = function(...,sort.bln=FALSE, reduce.bln=FALSE){
        mergedGrange.grn <- unlist(GenomicRanges::GRangesList(...))
        if(sort.bln){mergedGrange.grn <- sort(mergedGrange.grn )}
        if(reduce.bln){mergedGrange.grn <- GenomicRanges::reduce(mergedGrange.grn)}
        return(mergedGrange.grn)
}