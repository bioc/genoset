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

##' Calculate mBAF from BAF
##'
##' Calculate Mirrored B-Allele Frequence (mBAF) from B-Allele Frequency (BAF) as in
##' Staaf et al., Genome Biology, 2008.  BAF is converted to mBAF by folding around 0.5 so
##' that is then between 0.5 and 1. HOM value are then made NA to leave only HET values that
##' can be easily segmented. Values > hom.cutoff are made NA. Then, if genotypes (usually from
##' a matched normal) are provided as the matrix 'calls' additional HOMs can be set to NA. The
##' argument 'call.pairs' is used to match columns in 'calls' to columns in 'baf'.
##'
##' @param baf numeric matrix of BAF values
##' @param hom.cutoff numeric, values above this cutoff to be made NA (considered HOM)
##' @param calls matrix of NA, CT, AG, etc. genotypes to select HETs (in normals). Dimnames must match baf matrix.
##' @param call.pairs list, names represent target samples for HOMs to set to NA. Values represent columns in "calls" matrix.
##' @return numeric matix of mBAF values
##' @examples
##'    data(genoset)
##'    mbaf = baf2mbaf( baf(baf.ds), hom.cutoff=0.9 )
##'    calls = matrix(sample(c("AT","AA","CG","GC","AT","GG"),(nrow(baf.ds) * 2),replace=TRUE),ncol=2,dimnames=list(featureNames(baf.ds),c("K","L")))
##'    mbaf = baf2mbaf( baf(baf.ds), hom.cutoff=0.9, calls = calls, call.pairs = list(K="L",L="L") ) # Sample L is matched normal for tumor sample K, M only uses hom.cutoff
##'    assayDataElement(baf.ds,"mbaf") = baf2mbaf( baf(baf.ds), hom.cutoff=0.9 ) # Put mbaf back into the BAFSet object as a new element
##' @export
##' @author Peter M. Haverty
baf2mbaf <- function(baf, hom.cutoff=0.95, calls=NULL, call.pairs=NULL) {
  mbaf = abs(baf[,] - 0.5) + 0.5
  is.na(mbaf) <- mbaf > hom.cutoff
  
  if (!is.null(calls) && !is.null(call.pairs)) {

    # Use genotypes for/from samples specified by call.pairs to NA some HOMs
    if (! all(names(call.pairs) %in% colnames(baf)) ) {
      stop("call.pairs names and baf colnames mismatch\n")
    }
    if (! all(call.pairs %in% colnames(calls)) ) {
      stop("call.pairs values and calls colnames mismatch\n")
    }
    if ( ! identical( rownames(calls), rownames(baf) ) ) {
      stop("featureNames mismatch between calls and baf.")
    }
    
    # Check row matching between baf and calls.
    # Some calls rows will be missing from mbaf because PennCNV threw those features out.
    # Some rows of mbaf will not be in calls because some arrays have copy-only probes without calls
    # Can't subset both to intersection. Can only subset calls to those in mbaf
    # All baf values HET in calls will be NA, so no false positives
    # False negatives very possible, percent row overlap output will warn at least

    # NA all mbaf data for which there is a genotype and it is not HET

    hom.genotypes = c("AA","CC","GG","TT","AA","BB")
    if (nlevels(calls) > 0) { hom.genotypes = which( levels(calls) %in% hom.genotypes ) }
    is.na(mbaf[rownames(calls),names(call.pairs)]) <- calls[,unlist(call.pairs)] %in% hom.genotypes
  }
  return(mbaf)
}
