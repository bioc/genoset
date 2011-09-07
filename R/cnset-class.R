##### Class definition for CNSet, which will extend GenoSet and have matrices for baf and cn ###

#############
# Class CNSet
#############
##' @include genoset-class.R

##' @exportClass CNSet
setClass("CNSet", contains=c("GenoSet"))

setValidity("CNSet", function(object) {
  return(is.element("cn", assayDataElementNames(object)))
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
##' @param locData A RangedData object specifying feature chromosome
##' locations. Rownames are required to match featureNames.
##' @param cn numeric matrix of copy number data with rownames
##' matching sampleNames and colnames matching sampleNames
##' @param pData A data frame with rownames matching all data matrices
##' @param annotation character, string to specify chip/platform type
##' @param universe character, string to specify genome universe for locData
##' @param ... More matrix or DataFrame objects to include in assayData
##' @return A CNSet object
##' @export
##' @examples
##' test.sample.names = LETTERS[11:13]
##' probe.names = letters[1:10]
##' joe = CNSet(
##'    locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18"),
##'    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
##'    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
##'    annotation="SNP6"
##'    )
##' @author Peter M. Haverty
CNSet <- function(locData, cn, pData=NULL, annotation="", universe=NULL, ...) {
  object = initGenoSet(type="CNSet", locData=locData, pData=pData, annotation=annotation, universe=universe, cn=cn, ...)
  return(object)
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
##' @examples
##'   data(genoset)
##'   cn(cn.ds)  # Returns assayDataElement called "cn"
##'   cn(cn.ds) <- cn(cn.ds) + 5
##' @rdname cn
setGeneric("cn", function(object) standardGeneric("cn"))
##' @rdname cn
setMethod("cn", "CNSet", function(object) { return(object@assayData$cn) } )

##' @export "cn<-"
##' @rdname cn
setGeneric("cn<-", function(object,value) standardGeneric("cn<-") )
##' @rdname cn
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


#######
# Plots
#######

##' @author Peter M. Haverty
##' @rdname genoPlot
setMethod("genoPlot", signature(x="CNSet", y="ANY"),
          function(x, y, element="cn", chr=NULL, add=FALSE, ...) {
            # Plot copynumber
            if (element == "cn") {
              callNextMethod(x, y, element, chr, add, ylab="Log2 Ratio", ...)
              abline(h=0,lty=2,lwd=2)
            } else {
              callNextMethod(x, y, element, chr, add, ...)
            }
          } )

############
# Processing
############
