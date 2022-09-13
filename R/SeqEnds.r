#' Get all sequences lengths.
#'
#' SeqEnds
#' Get all sequences lengths for each ranges of a GRanges object.
#' @param x.gnr <GRanges>: a GRanges object.
#' @return an integer vector.
#' @examples
#' GRange.grn <- GenomicRanges::GRanges(
#'     seqnames = S4Vectors::Rle(c("chr1", "chr2", "chr1"), c(1, 3, 1)),
#'     ranges = IRanges::IRanges(101:105, end = 111:115, names = letters[1:5]),
#'     strand = S4Vectors::Rle(BiocGenerics::strand(c("-", "+", "*", "+")), c(1, 1, 2, 1)),
#'     seqinfo = c(chr1=200, chr2=300),
#'     score = 1:5
#' )
#' SeqEnds(GRange.grn)
SeqEnds <- function(x.gnr){GenomeInfoDb::seqlengths(x.gnr)[as.character(GenomeInfoDb::seqnames(x.gnr))]}