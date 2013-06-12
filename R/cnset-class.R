##### Class definition for CNSet, which will extend GenoSet and have matrices for baf and cn ###

#############
# Class CNSet
#############
##' @include genoset-class.R

##' @exportClass CNSet
setClass("CNSet", contains=c("GenoSet"))

setValidity("CNSet", function(object) {
  .Defunct("GenoSet", msg="The CNSet class is defunct. Please use GenoSet. CNSet only added the cn getter/setter functions, which are redundant with x[, , 'cn'] now.")
})

##' @exportMethod show
##' @rdname show
##' @aliases show,CNSet-method
setMethod("show","CNSet",
          function(object) {
            .Defunct("GenoSet", msg="The CNSet class is defunct. Please use GenoSet. CNSet only added the cn getter/setter functions, which are redundant with x[, , 'cn'] now.")
          })

##' Create a CNSet object
##'
##' This function is the preferred method for creating a new CNSet object. Users are
##' generally discouraged from calling "new" directly. This CNSet function enforces
##' the requirement for a "cn" matrix.  This and any other "..." arguments will
##' become part of the assayData slot of the resulting object. "..." can be matrices
##' or DataFrame objects (from the IRanges package). This function passes
##' control to the "initGenoSet" method which performs argument checking including
##' dimname matching among relevant slots and sets everything to genome order. Genome
##' order can be disrupted by "[" or "[[" calls and will be checked by methods that
##' require it.
##'
##' The CNSet class is defunct. Please use GenoSet. CNSet only added the cn getter/setter functions, which are redundant with x[, , 'cn'] now.
##' 
##' @param locData A GRanges or RangedData object specifying feature chromosome
##' locations. featureNames (names or rownames) are required to match featureNames of matrices.
##' @param cn numeric matrix of copy number data with rownames
##' matching featureNames and colnames matching sampleNames
##' @param pData A data frame with rownames matching all data matrices
##' @param annotation character, string to specify chip/platform type
##' @param universe character, string to specify genome universe for locData. Overrides any universe/genome data in locData.
##' @param assayData assayData, usually an environment
##' @param ... More matrix or DataFrame objects to include in assayData
##' @return A CNSet object
##' @export
##' @author Peter M. Haverty
##' @aliases CNSet-defunct
CNSet <- function(locData, cn=NULL, pData=NULL, annotation="", universe, assayData=NULL, ...) {
  .Defunct("GenoSet", msg="The CNSet class is defunct. Please use GenoSet. CNSet only added the cn getter/setter functions, which are redundant with x[, , 'cn'] now.")
}

#########
# Methods
#########

#####################
# Getters and Setters
#####################


##' Get or Set the cn assayData slot
##'
##' Get or Set the cn assayData slot
##' 
##' @param object A BAFset object
##' @return matrix
##' @author Peter M. Haverty
##' @export cn
##' @rdname cn
setGeneric("cn", function(object) standardGeneric("cn"))
##' @rdname cn
##' @aliases cn,CNSet-method
setMethod("cn", "CNSet", function(object) { return(object@assayData$cn) } )

##' @export "cn<-"
##' @rdname cn
setGeneric("cn<-", function(object,value) standardGeneric("cn<-") )
##' @rdname cn
##' @aliases cn<-,CNSet,matrix-method
setMethod("cn<-", signature(object="CNSet", value="matrix"),
                 function(object,value) assayDataElementReplace(object, "cn", value))

###########
# Coersions
###########
setAs("ExpressionSet","CNSet", def=
      function(from, to) {
        if (! all( c("space","start","end") %in% colnames(fData(from)))) {
          stop("Feature Data of supplied ExpressionSet must have space (chromosome), start, and end columns")
        }
        locData = RangedData(
          ranges=IRanges(start=fData(from)$start, end=fData(from)$end, names=rownames(fData(from))),
          space=fData(from)$space)
        
        new(to, cn=exprs(from), phenoData=phenoData(from), locData=locData,
            experimentData=experimentData(from), annotation=annotation(from))
      })

