#' GenomicSystem
#'
#' Convert numbers into writing with order of magnitude (Kilo, Mega, Giga) and vice versa.
#' @param x <character or numeric>: the number to convert
#' @param digits.num <numeric>: integer indicating the number of significant digits to be used. See ?signif() for more informations. (Default 3)
#' @return the converted number
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
            stringr::str_detect(x,"G") ~ stringr::str_remove_all(x,"[:alpha:]*") %>% as.numeric %>% magrittr::multiply_by(10**9),
            stringr::str_detect(x,"M") ~ stringr::str_remove_all(x,"[:alpha:]*") %>% as.numeric %>% magrittr::multiply_by(10**6),
            stringr::str_detect(x,"K") ~ stringr::str_remove_all(x,"[:alpha:]*") %>% as.numeric %>% magrittr::multiply_by(10**3),
            suppressWarnings(!is.na(as.numeric(x))) ~ suppressWarnings(as.numeric(x))
        )
    }
}