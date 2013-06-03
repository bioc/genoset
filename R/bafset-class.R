##### Class definition for BAFSet, which will extend GenoSet and have matrices for baf and lrr ###

################
# Class BAFSet
################
##' @include genoset-class.R
##' @include cnset-class.R

##' @exportClass BAFSet
setClass("BAFSet", contains=c("GenoSet"))

setValidity("BAFSet", function(object) {
  .Defunct("", msg="The BAFSet class is defunct. Please use GenoSet. BAFSet only added the baf/lrr getter/setter functions, which are redundant with x[, , 'baf'] and x[, , 'lrr'] now.")
})

setMethod("show","BAFSet",
          function(object) {
            .Defunct("", msg="The BAFSet class is defunct. Please use GenoSet. BAFSet only added the baf/lrr getter/setter functions, which are redundant with x[, , 'baf'] and x[, , 'lrr'] now.")
          })

##' Create a BAFSet object
##'
##' This function is the preferred method for creating a new BAFSet object. Users are
##' generally discouraged from calling "new" directly. This BAFSet function enforces
##' the requirement for "lrr" and "baf" matrices.  These and any other "..." arguments will
##' become part of the assayData slot of the resulting object. "..." can be matrices or
##' DataFrame objects (from the IRanges package). This function passes
##' control to the "initGenoSet" method which performs argument checking including
##' dimname matching among relevant slots and sets everything to genome order. Genome
##' order can be disrupted by "[" or "[[" calls and will be checked by methods that
##' require it.
##'
##' The BAFSet class is defunct. Please use GenoSet. BAFSet only added the baf/lrr getter/setter functions, which are redundant with x[, , 'baf'] and x[, , 'lrr'] now.
##' 
##' @param locData A GRanges or RangedData object specifying feature chromosome
##' locations. featureNames (names or rownames) are required to match featureNames of assayData.
##' @param lrr numeric matrix of copy number data with rownames
##' matching featureNames and colnames matching sampleNames
##' @param baf numeric matrix of B-Allele Frequency data with rownames
##' matching featureNames and colnames matching sampleNames
##' @param pData A data frame with rownames matching all data matrices
##' @param annotation character, string to specify chip/platform type
##' @param universe character, a string to specify the genome universe for locData. Overrides any universe/genome data in locData.
##' @param assayData assayData, usually an environment
##' @param ... More matrix or DataFrame objects to include in assayData slot
##' @return A BAFSet object
##' @export 
##' @author Peter M. Haverty
##' @aliases BAFSet-defunct
##' @examples
##'   test.sample.names = LETTERS[11:13]
##'   probe.names = letters[1:10]
##'   locData.rd = RangedData(ranges=IRanges(start=c(1,4,3,2,5:10),width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18")
##'   bs = BAFSet(
##'     locData=locData.rd,
##'     lrr=matrix(1:30,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
##'     baf=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
##'     pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
##'     annotation="SNP6"
##' )
##' @seealso bafset-class, genoset-class
BAFSet <- function(locData, lrr=NULL, baf=NULL, pData=NULL, annotation="", universe, assayData=NULL, ...) {
  .Defunct("GenoSet", msg="The BAFSet class is defunct. Please use GenoSet. BAFSet only added the baf/lrr getter/setter functions, which are redundant with x[, , 'baf'] and x[, , 'lrr'] now.")
  if (!is.null(assayData)) {
    if (! all(c("lrr","baf") %in% assayDataElementNames(assayData))) {
      stop("If assayData is specified, it must contain elements called 'lrr' and 'baf'.")
    }
    object = initGenoSet(type="BAFSet", locData=locData, pData=pData, annotation=annotation, universe=universe, assayData=assayData, ...)
  } else {
    object = initGenoSet(type="BAFSet", locData=locData, pData=pData, annotation=annotation, universe=universe, lrr=lrr, baf=baf, ...)
  }
  return(object)
}

##########
## Methods
##########

#####################
# Getters and Setters
#####################

##' Get or Set the baf assayData slot
##'
##' Get or Set the baf assayData slot
##'
##' @title Get baf data
##' @param object A BAFset object
##' @return matrix
##' @author Peter M. Haverty
##' @export baf
##' @examples
##'   data(genoset)
##'   baf(baf.ds)  # Returns assayDataElement called "baf"
##'   baf(baf.ds) <- baf2mbaf( baf(baf.ds) )
##' @rdname baf
setGeneric("baf", function(object) standardGeneric("baf"))
##' @rdname baf
##' @aliases baf,BAFSet-method
setMethod("baf", "BAFSet", function(object) { return(object@assayData$baf) } )

##' Get or Set the lrr assayData slot
##'
##' Get or Set the lrr assayData slot
##'
##' @title Get lrr data
##' @param object A BAFset object
##' @return matrix
##' @author Peter M. Haverty
##' @export lrr
##' @rdname lrr
##' @examples
##'   data(genoset)
##'   lrr(baf.ds)  # Returns assayDataElement called "lrr"
##'   lrr(baf.ds) <- lrr(baf.ds) + 0.1
setGeneric("lrr", function(object) standardGeneric("lrr"))
##' @rdname lrr
##' @aliases lrr,BAFSet-method
setMethod("lrr", "BAFSet", function(object) { return(object@assayData$lrr) } )

##' @export "baf<-"
##' @rdname baf
setGeneric("baf<-", function(object,value) standardGeneric("baf<-") )
##' @rdname baf
##' @aliases baf<-,BAFSet,matrix-method
setMethod("baf<-", signature(object="BAFSet", value="matrix"),
                 function(object,value) assayDataElementReplace(object, "baf", value))

##' @export "lrr<-"
##' @rdname lrr
setGeneric("lrr<-", function(object,value) standardGeneric("lrr<-") )
##' @rdname lrr
##' @aliases lrr<-,BAFSet,matrix-method
setMethod("lrr<-", signature(object="BAFSet", value="matrix"),
                 function(object,value) assayDataElementReplace(object, "lrr", value))

###########
# Coercions
###########
setAs("BAFSet","CNSet", def=
      function(from, to) {
        new(to, cn=lrr(from), phenoData=phenoData(from), locData=locData(from),
            experimentData=experimentData(from), annotation=annotation(from))
      })
setAs("BAFSet","ExpressionSet", def=
      function(from, to) {
        features = data.frame(chr=chr(from), start=start(from), end=end(from), row.names=featureNames(from))
        features = new("AnnotatedDataFrame",data=features)
        new("ExpressionSet", exprs=lrr(from), baf=baf(from), phenoData=phenoData(from),
            experimentData=experimentData(from), annotation=annotation(from),
            featureData=features)
      })

##' Make a pair of ExpressionSets from a BAFSet
##'
##' Often it is convenient to have a more standard "ExpressionSet" rather than
##' a BAFSet.  For example, when using infrastructure dependent on the
##' ExpressionSet slots, like limma or ExpressionSetOnDisk. This will create a
##' list of two ExpressionSets, one each for the baf and lrr data. To make a single
##' ExpressionSet, with the lrr data in the exprs slot and the baf data as an additional
##' member of assayData, use the standard coercion eset = as(bafset,"ExpressionSet").
##'
##' BAFSEt.toExpressionSets has been defunct. Please use as(x, 'ExpressionSet').
##' 
##' @param bs A BAFset object
##' @return A list with one ExpressionSet each for the baf and lrr data in the BAFSet object
##' @export
##' @examples
##'   data(genoset)
##'   eset.list = BAFSet.to.ExpressionSets(baf.ds)
##' @author Peter M. Haverty
##' @aliases BAFSet.to.ExpressionSets-defunct
BAFSet.to.ExpressionSets <- function(bs) {
  .Defunct("as", msg="BAFSEt.toExpressionSets has been defunct. Please use as(x, 'ExpressionSet').")
}

