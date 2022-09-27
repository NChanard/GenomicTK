#' Convert numbers of base pairs into string and vice versa.
#'
#' GenomicSystem
#' Convert numbers of base into string with order of magnitude (Kbp, Mbp, Gbp) and vice versa.
#' @param x <character or numeric>: the number to convert or string to convert.
#' @param digits.num <integer>: The number of significant digits to be used. See signif() for more informations. (Default 3)
#' @return the converted number or string.
#' @examples
#' GenomicSystem(1540,3)
#' GenomicSystem(1540,2)
#' GenomicSystem("1Mbp")
#' GenomicSystem("1Kbp")
#' GenomicSystem("1k")
GenomicSystem <- function(x, digits.num=3) {
    if(is.numeric(x)){
        dplyr::case_when(
            x >= 1e9   ~ paste0(signif(x * 10**(-9), digits.num), "Gbp"),
            x >= 1e6   ~ paste0(signif(x * 10**(-6), digits.num), "Mbp"),
            x >= 1e3   ~ paste0(signif(x * 10**(-3), digits.num), "Kbp"),
            x >= 0     ~ paste0(signif(x * 10**( 0), digits.num), "Bp")
        )
    }else if(is.character(x)){
        x <- toupper(x)
        dplyr::case_when(
            grepl(x=x, pattern="G") ~ as.numeric(gsub(x=x, pattern="[a-z]", ignore.case=TRUE, replacement=""))*(10**9),
            grepl(x=x, pattern="M") ~ as.numeric(gsub(x=x, pattern="[a-z]", ignore.case=TRUE, replacement=""))*(10**6),
            grepl(x=x, pattern="K") ~ as.numeric(gsub(x=x, pattern="[a-z]", ignore.case=TRUE, replacement=""))*(10**3),
            suppressWarnings(!is.na(as.numeric(x))) ~ suppressWarnings(as.numeric(x))
        )
    }
}