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

##' Correct copy number for GC content
##'
##' Copy number estimates from various platforms show "Genomic Waves" (Diskin et al.,
##' Nucleic Acids Research, 2008) where copy number trends with local GC content.
##' This function regresses copy number on GC percentage and removes the effect
##' (returns residuals). GC content should be smoothed along the genome in wide
##' windows >= 100kb.
##'
##' @param ds numeric matrix of copynumber or log2ratio values, samples in columns
##' @param gc numeric vector, GC percentage for each row of ds, must not have NAs
##' @param retain.mean logical, center on zero or keep same mean?
##' @return numeric matrix, residuals of ds regressed on gc
##' @export gcCorrect
##' @family "gc content"
##' @examples
##'   gc = runif(n=100, min=1, max=100)
##'   ds = rnorm(100) + (0.1 * gc)
##'   gcCorrect(ds, gc)
##' @author Peter M. Haverty
gcCorrect <- function(ds, gc, retain.mean=TRUE) {
  if (!requireNamespace("stats",quietly=TRUE)) {
    stop("Failed to require stats package.\n")
  }
  ds = na.exclude(ds)
  if ("na.action" %in% names(attributes(ds))) {
    gc = gc[ - attr(ds, "na.action") ]
  }
  mm = cbind(rep.int(1, length(gc)), gc)
  fit = stats::lm.fit(mm, ds)
  fit$na.action = attr(ds, "na.action")
  ds.fixed = stats::residuals(fit)
  if (retain.mean == TRUE) {
    if (is.null(dim(ds))) {
      ds.fixed = ds.fixed + (sum(ds,na.rm=TRUE)/length(ds))
    } else {
      ds.fixed = sweep(ds.fixed,2,colMeans(ds,na.rm=TRUE),FUN="+")
    }
  }
  return(ds.fixed)
}
