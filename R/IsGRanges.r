#' IsGRange
#'
#' Reports whether x is a GRanges object
#' @param x.any <any>: An object to test
#' @return a boolean that say if x is a GRange object
#' @examples
#' GRange.grn <- GenomicRanges::GRanges(
#'     seqnames = S4Vectors::Rle(c("chr1", "chr2", "chr1"), c(1, 3, 1)),
#'     ranges = IRanges::IRanges(101:105, end = 111:115, names = letters[1:5]),
#'     strand = S4Vectors::Rle(BiocGenerics::strand(c("-", "+", "*", "+")), c(1, 1, 2, 1)),
#'     seqinfo = c(chr1=200, chr2=300),
#'     score = 1:5
#' )
#' IsGRange(GRange.grn)
IsGRange <- function(GRange.grn){inherits(GRange.grn, "GRanges")}