#' Convert String in GRanges object
#'
#' StrToGRanges
#' Convert ranges describe with string (i.e seqname:start-end:strand) in GRanges object.
#' @param x.chr_vec <character>: strings to convert on GRanges.
#' @return a GRanges object.
#' @examples
#' StrToGRanges("chr1:1-100:+")
#' StrToGRanges(c("chr1:1-100:+","chr2:400-500:-","chr1:10-50:*"))
StrToGRanges <- function(x.chr_vec){
    x.gnr <- lapply(x.chr_vec, function(x.chr){
        x.chr <- unlist(strsplit(x.chr ,":"))
        seqnames.chr <- x.chr[1]
        ranges.num <- strsplit(x.chr[2],"-") |> unlist() |> as.numeric()
        start.num <- ranges.num[1]
        end.num <- ifelse(is.na(ranges.num[2]),yes=start.num, no=ranges.num[2]) 
        strand.chr <- ifelse(is.na(x.chr[3]),yes="*", no=x.chr[3])
        x.gnr <- GenomicRanges::GRanges(
            seqnames=seqnames.chr,
            ranges=IRanges::IRanges(start=start.num, end=end.num),
            strand=strand.chr)
        return(x.gnr)
    }) |> GenomicTK::MergeGRanges() 
    S4Vectors::mcols(x.gnr)$names <- names(x.chr_vec)
    return(x.gnr)
}