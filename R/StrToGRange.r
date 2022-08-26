#' StrToGRange
#'
#' Convert ranges written in Genomic system (i.e seqname:start-end:strand) in GRanges object.
#' @param x.chr_vec <character>: strings to convert on GRanges
#' @return A GRanges object
#' @examples
#' StrToGRange("chr1:1-100:+")
#' StrToGRange(c("chr1:1-100:+","chr2:400-500:-","chr1:10-50:*"))
StrToGRange <- function(x.chr_vec){
    x.grn <- lapply(x.chr_vec, function(x.chr){
        x.chr %<>% stringr::str_split(":") %>% unlist
        seqnames.chr <- x.chr[1]
        ranges.num <- x.chr[2] %>% stringr::str_split("-") %>% unlist %>% as.numeric
        start.num <- ranges.num[1]
        end.num <- ifelse(is.na(ranges.num[2]),yes=start.num, no=ranges.num[2]) 
        strand.chr <- ifelse(is.na(x.chr[3]),yes="*", no=x.chr[3])
        GenomicRanges::GRanges(
            seqnames=seqnames.chr,
            ranges=IRanges::IRanges(start=start.num, end=end.num),
            strand=strand.chr) %>%
        return(.)
    }) %>% GenomicTK::MergeGRanges(.) 
    S4Vectors::mcols(x.grn)$names <- names(x.chr_vec)
    return(x.grn)
}