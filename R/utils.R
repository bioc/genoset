##' Take vector or matrix of log2 ratios, convert to copynumber

##' Utility function for converting log2ratio units (zero is normal) to copynumber units (two is normal)
##' @param x numeric data in log2ratio values
##' @return data of same type as "x" transformed into copynumber units
##' @export
##' @seealso cn2lr
##' @author Peter M. Haverty \email{phaverty@@gene.com}
lr2cn <- function(x) {
  return( 2 ^ (x + 1) )
}

##' Take vector or matrix of copynumber values, convert to log2ratios
##' Utility function for converting copynumber units (2 is normal) to log2ratio units (two is normal)
##' @param x numeric data in copynumber units
##' @return data of same type as "x" transformed into log2ratio units
##' @export
##' @seealso lr2cn
##' @author Peter M. Haverty \email{phaverty@@gene.com}
cn2lr <- function(x) {
  return ( log2(x) - 1 )
}
