####  Class definition for GenoSet, which will extend eSet
######   GenoSet will provide a locData slot containing a RangedData object from the IRanges
######   package to hold genome locations of the features and allow for easy subsetting
######   by location.
######   Intended to be subset by other classes to add one or more data matrices to
######   the assayData slot.

##' GenoSet: An eSet for data with genome locations
##' 
##' Load, manipulate, and plot copynumber and BAF data. GenoSet class
##' extends eSet by adding a "locData" slot for a GenomicRanges object.
##' This object contains feature genome location data and
##' provides for efficient subsetting on genome location.
##' Genoset also implements an number of  convenience functions for processing of copy number and B-Allele Frequency data and for working with segmented data.
##'
##' @docType package
##' @name genoset-package
##' @aliases genoset genoset-package
##' @seealso genoset-datasets GenoSet
##' 
##' @importClassesFrom Biobase AnnotatedDataFrame AssayData eSet ExpressionSet MIAME Versioned VersionedBiobase
##' @importClassesFrom IRanges DataFrame RangedData Rle
##' @importClassesFrom GenomicRanges GRanges
##'
##' @importMethodsFrom GenomicRanges seqnames seqlevels names "names<-" length width genome "genome<-"
##' @importMethodsFrom Biobase annotation experimentData exprs fData featureNames "featureNames<-" phenoData sampleNames "sampleNames<-"
##' @importMethodsFrom IRanges as.data.frame as.list as.matrix cbind colnames "colnames<-" elementLengths end findOverlaps gsub
##' @importMethodsFrom IRanges intersect is.unsorted lapply levels mean na.exclude nrow order paste ranges Rle rownames
##' @importMethodsFrom IRanges "rownames<-" runLength runValue sapply space start unlist universe "universe<-"
##'
##' @importFrom Biobase assayDataElement assayDataElementNames assayDataElementReplace assayDataNew annotatedDataFrameFrom
##' @importFrom graphics abline axis axTicks box mtext plot.new plot.window points segments
##' @importFrom IRanges DataFrame IRanges RangedData "%over%"
##' @importFrom GenomicRanges seqlengths GRanges
##'
##' @import methods
##' @import BiocGenerics
##' 
##' @useDynLib genoset
NULL

##' Deprecated genoset features
##'
##' GenoSet is moving towards using GenomicRanges instead of RangedData. We are also getting rid of dependencies on eSet for a potential switch to an underlying SummarizedExperiment.
##' @name genoset-deprecated
##' @aliases genoset-deprecated
NULL

##' Defunct genoset features
##'
##' The CNSet and BAFSet classes are defunct.  They only really added getter/setter methods for specific assayDataElements,
##' so they are now redundant with the preferred method of using the assayDataElement name as the third argument to bracket, e.g.
##' \code{x[i, j, "lrr"]}. Accordingly \code{BAFSet.to.ExpressionSets} is also defunct.
##'
##' Additionally, names, ranges, and space on a GenoSet are also defunct. In an effort to make a consistent API for either RangedData or
##' GRanges in the locData slot, we recommend using \code{chrNames} for \code{names} and \code{chr} for \code{space}.
##' @name genoset-defunct
##' @aliases genoset-defunct
NULL

###############
# Class GenoSet
###############

##' @exportClass RangedDataOrGenomicRanges
##' @exportClass GenoSet
##' @exportClass RangedDataOrGenoSet
##' @exportClass RangedDataOrGenoSetOrGenomicRanges
setClassUnion("RangedDataOrGenomicRanges",c("RangedData","GenomicRanges"))
setClass("GenoSet", contains=c("eSet"), representation=representation(locData="RangedDataOrGenomicRanges"))
setClassUnion("RangedDataOrGenoSet",c("RangedData","GenoSet"))
setClassUnion("RangedDataOrGenoSetOrGenomicRanges",c("RangedData","GenoSet","GenomicRanges"))

setValidity("GenoSet", function(object) {
  return( all( rownames(locData(object)) == rownames(assayData(object)) ) )
})

##' Create a GenoSet or derivative object
##'
##' This function is the preferred method for creating a new GenoSet object. Users are
##' generally discouraged from calling "new" directly. The "..." argument is for any number of matrices of matching size that will
##' become part of the assayData slot of the resulting object. This function passes
##' control to the "genoSet" object which performs argument checking including
##' dimname matching among relevant slots and sets everything to genome order. Genome
##' order can be disrupted by "[" calls and will be checked by methods that
##' require it.
##' 
##' @param type character, the type of object (e.g. GenoSet, BAFSet, CNSet) to be created
##' @param locData A GRanges or RangedData object specifying feature chromosome
##' locations. rownames are required to match assayData.
##' @param pData A data frame with rownames matching colnames of all assayDataElements
##' @param annotation character, string to specify chip/platform type
##' @param universe character, a string to specify the genome universe for locData, overrides universe/genome data in locData
##' @param assayData assayData, usually an environment
##' @param ... More matrix or DataFrame objects to include in assayData
##' @return A GenoSet object or derivative as specified by "type" arg
##' @examples 
##'   test.sample.names = LETTERS[11:13]
##'   probe.names = letters[1:10]
##'   gs = GenoSet(
##'      locData=GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4))),
##'      cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
##'      pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
##'      annotation="SNP6"
##'   )
##' @author Peter M. Haverty
initGenoSet <- function(type, locData, pData=NULL, annotation="", universe, assayData=NULL, ...) {
  # Function to clean up items for slots and call new for GenoSet and its children
  # ... will be the matrices that end up in assayData

  if (! missing(universe)) {
    warning("Use of the universe argument for GenoSet creation is deprecated. Please note the genome in the location data object.")
    universe(locData) = universe
  }

  # Check/set genome order of locData
  if ( ! isGenomeOrder(locData, strict=TRUE) ) {
    locData = toGenomeOrder(locData, strict=TRUE )
  }
  clean.loc.rownames = rownames(locData)

  # Create assayData
  if (is.null(assayData)) {
    # Crib most of assayDataNew, skipping unnaming of dimnames to keep BigMatrix happy
    ad = new.env(parent=emptyenv())
    arglist <- list(...)
    if ((length(arglist) > 0L) && ((is.null(names(arglist))) || any(names(arglist) == ""))) { stop("all arguments must be named") }
    for (nm in names(arglist)) {
      ad[[nm]] <- arglist[[nm]]
    }
    msg <- assayDataValidMembers(ad)
    if (!is.logical(msg)) { stop(msg) }
  } else {
    ad = assayData
  }
  clean.rownames = featureNames(ad)

  if (length(clean.rownames) != length(clean.loc.rownames)) {
    stop("Row number mismatch for assayData and locData")
  }

  # Set row order to match locData, already know all ad elements have same row names
  if ( ! all(clean.loc.rownames == clean.rownames) ) {
    if (! setequal(clean.loc.rownames, clean.rownames)) {
      stop("Row name set mismatch for locData and assayData")
    } else {
      for (  ad.name in assayDataElementNames(ad) ) {
        ad[[ad.name]] = ad[[ad.name]][clean.loc.rownames,]        
      }
    }
  }

  # Check colnames of all data matrices identical and set to same order if necessary
  # AssayDataValidMembers does not do this for some reason
  first.name = assayDataElementNames(ad)[1]
  for (mat.name in assayDataElementNames(ad)[-1]) {
    if ( ! all( colnames(ad[[mat.name]]) == colnames(ad[[first.name]])) ) {
      if (! setequal(colnames(ad[[mat.name]]), colnames(ad[[first.name]]) ) ) {
        stop(paste("Mismatch between rownames of first data matrix and", mat.name))
      } else {
        ad[[mat.name]] == ad[[mat.name]][,colnames(ad[[first.name]])]
      }
    }
  }

  # Done editing assayData members, lock
  lockEnvironment(ad, bindings=TRUE)
  
  # Create or check phenoData
  if (is.null(pData)) {
    pData = data.frame(row.names=sampleNames(ad),check.names=FALSE)
  } else {
    if ( ! setequal( rownames(pData), sampleNames(ad) ) ) {
      stop( "Mismatch between sampleNames and rownames of pData" )
    }
    if ( any( sampleNames(ad) != rownames(pData) ) ) {
      pData = pData[ sampleNames(ad), ]
    }
  }
  pd = new("AnnotatedDataFrame",data=pData)

  # Create object
  rownames(locData) = NULL
  object = new(type, locData=locData, annotation=annotation, phenoData=pd, assayData=ad)
  return(object)
}

##' Create a GenoSet object
##'
##' This function is the preferred method for creating a new GenoSet object. Users are
##' generally discouraged from calling "new" directly. Any "..." arguments will
##' become part of the assayData slot of the resulting object. "..." can be matrices
##' or DataFrame objects (from IRanges). This function passes
##' control to the "initGenoSet" method which performs argument checking including
##' dimname matching among relevant slots and sets everything to genome order. Genome
##' order can be disrupted by "[" calls and will be checked by methods that
##' require it.
##' 
##' @param locData A RangedData object specifying feature chromosome
##' locations. Rownames are required to match featureNames.
##' @param pData A data frame with rownames matching all data matrices
##' @param annotation character, string to specify chip/platform type
##' @param universe character, a string to specify the genome universe for locData
##' @param assayData assayData, usually an environment
##' @param ... More matrix or DataFrame objects to include in assayData
##' @return A GenoSet object
##' @examples
##' test.sample.names = LETTERS[11:13]
##' probe.names = letters[1:10]
##' gs = GenoSet(
##'    locData=GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4))),
##'    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
##'    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
##'    annotation="SNP6"
##' )
##' @export GenoSet
##' @author Peter M. Haverty
GenoSet <- function(locData, pData=NULL, annotation="", universe, assayData=NULL, ...) {
  object = initGenoSet(type="GenoSet", locData=locData, pData=pData, annotation=annotation, universe=universe, assayData=assayData,...)
  return(object)
}

#########
# Methods
#########

#####################
# Getters and Setters
#####################

##' Genome version
##'
##' The genome positions of the features in locData. The UCSC notation (e.g. hg18, hg19, etc.) should be used.
##' @title Get and set the genome universe annotation.
##' @param x GenoSet
##' @return character, e.g. hg19
##' @author Peter M. Haverty
##' @exportMethod universe
##' @exportMethod genome
##' @examples
##'   data(genoset)
##'   genome(genoset.ds)
##'   genome(genoset.ds) = "hg19"
setMethod("genome", "GenoSet", function(x) {
  return(genome(x@locData))
})
setMethod("genome<-", "GenoSet", function(x, value) {
  genome(x@locData) = value
  return(x)
})
setMethod("universe", "GenoSet", function(x) {
  .Deprecated(new="genome", msg="RangedData is being replaced with GenomicRanges. Please use the genome method.")
  return(universe(x@locData))
} )
setMethod("universe", "GRanges", function(x) {
  .Deprecated(new="genome", msg="RangedData is being replaced with GenomicRanges. Please use the genome method.")
  if (length(unique(genome(x))) != 1) {
    warning("Taking first element of GRanges genome as universe.")
  }
  return(unname(genome(x)[1]))
} )
setMethod("universe<-", signature(x="GenoSet"),
                 function(x,value) {
                   .Deprecated(new="genome", msg="RangedData is being replaced with GenomicRanges. Please use the genome method.")
                   universe(x@locData) = value
                   return(x)
                   })
setMethod("universe<-", signature(x="GRanges"),
          function(x,value) {
            .Deprecated(new="genome", msg="RangedData is being replaced with GenomicRanges. Please use the genome method.")
            genome(x) = value
            return(x)
          })

##' Get colnames from a GenoSet
##'
##' Get colnames from a GenoSet
##' 
##' @param object GenoSet
##' @return character vector with names of samples
##' @examples
##'   data(genoset)
##'   head(colnames(genoset.ds))
##' @exportMethod sampleNames
##' @exportMethod colnames
setMethod("colnames", signature(x="GenoSet"),
          function(x) {
            rownames(pData(x))
          })
#setMethod("colnames<-", signature(x="GenoSet"),
#          function(x, value) {
#            rownames(pData(x)) = value
#            return(object)
#          })
setMethod("sampleNames", signature(object="GenoSet"),
          function(object) {
            .Deprecated(new="colnames", msg="Please use colnames. We are switching away from eSet-specific methods.")
            rownames(pData(object))
          })
setMethod("sampleNames<-", signature(object="GenoSet"),
          function(object, value) {
            .Deprecated(new="colnames<-", msg="Please use colnames. We are switching away from eSet-specific methods.")
            object = callNextMethod(object, value)
            return(object)
          })

##' Get rownames from RangedData, GRanges, or GenoSet
##'
##' Get rownames from RangedData, GRanges, or GenoSet.
##' 
##' @param object GRanges, RangedData, or GenoSet
##' @return character vector with names rows/features
##' @author Peter M. Haverty
##' @examples
##'   data(genoset)
##'   head(rownames(locData.gr))
##'   head(rownames(genoset.ds))
##' @exportMethod featureNames
##' @exportMethod rownames
##' @rdname rownames-methods
setMethod("featureNames", signature(object="GenoSet"),
          function(object) {
            .Deprecated(new="rownames", msg="Please use rownames. We are switching away from eSet-specific methods.")
            return(unname(featureNames(featureData(object))))
          })
##' @rdname rownames-methods
setMethod("featureNames", signature(object="GRanges"),
          function(object) {
            .Deprecated(new="rownames", msg="Please use rownames. We are switching away from eSet-specific methods.")
            rownames(object)
          })
##' @rdname rownames-methods
setMethod("featureNames", signature(object="RangedData"),
          function(object) {
            .Deprecated(new="rownames", msg="Please use rownames. We are switching away from eSet-specific methods.")
            rownames(object)
          })
##' @rdname rownames-methods
setMethod("rownames", signature(x="GRanges"),
          function(x) {
            names(x)
          })
##' @rdname rownames-methods
setMethod("rownames", signature(x="GenoSet"),
          function(x) {
            return(unname(featureNames(featureData(x))))
          })
##' @rdname rownames-methods
setMethod("featureNames<-",
          signature=signature(object="GenoSet", value="ANY"),
          function(object, value) {
            .Deprecated(new="rownames<-", msg="Please use rownames. We are switching away from eSet-specific methods.")
            object = callNextMethod(object, value)
            return(object)
          })
##' @rdname rownames-methods
setMethod("featureNames<-",
          signature=signature(object="GRanges", value="ANY"),
          function(object, value) {
            .Deprecated(new="rownames<-", msg="Please use rownames. We are switching away from eSet-specific methods.")
            rownames(object) = value
            return(object)
          })
##' @rdname rownames-methods
setMethod("featureNames<-",
          signature=signature(object="RangedData", value="ANY"),
          function(object, value) {
            .Deprecated(new="rownames<-", msg="Please use rownames. We are switching away from eSet-specific methods.")
            rownames(object) = value
            return(object)
          })
##' @rdname rownames-methods
setMethod("rownames<-",
                 signature=signature(x="GRanges", value="ANY"),
                 function(x, value) {
                   names(x) = value
                   return(x)
                 })

##' Access the feature genome position info
##'
##' The position information for each probe/feature is stored as an IRanges RangedData object.
##' The locData functions allow this data to be accessed or re-set.
##'
##' @param object GenoSet
##' @param value RangedData describing features
##' @section Methods:
##' \describe{
##' \item{\code{signature(object = "GenoSet")}}{
##' Get location data.
##' }
##' \item{\code{signature(object = "GenoSet", value = "RangedData")}}{
##' Set location data.
##' }}
##' @export locData
##' @export "locData<-"
##' @author Peter M. Haverty
##' @docType methods
##' @examples
##' data(genoset)
##' rd = locData(genoset.ds)
##' locData(genoset.ds) = rd
##' @aliases locData-methods
##' @aliases locData<--methods
##' @aliases locData,GenoSet-method
##' @aliases locData<-,GenoSet,RangedDataOrGenomicRanges-method
##' @aliases locData
##' @aliases locData<-
##' @return A GenoSet object
setGeneric("locData", function(object) standardGeneric("locData"))
setMethod("locData", "GenoSet", function(object) {
  locs = slot(object,"locData")
  rownames(locs) = rownames(object)
  return(locs)
} )
setGeneric("locData<-", function(object,value) standardGeneric("locData<-") )
setMethod("locData<-", signature(object="GenoSet", value="RangedDataOrGenomicRanges"),
                 function(object,value) {
                   if (! all( rownames(value) %in% rownames(object))) {
                       stop("Can not replace locData using rownames not in this GenoSet")
                     }
                   if (! all(rownames(value) == rownames(object))) {
                     object = object[rownames(value), ]
                   }
                   rownames(value) = NULL
                   slot(object,"locData") = value
                   return(object)
                   })

###########################################
# Shared API between GenoSet and RangedData
###########################################

##' Get start of location for each feature
##'
##' Get start of location for each feature
##' @param x GenoSet
##' @return integer
##' @author Peter M. Haverty
##' @aliases start,GenoSet-method
setMethod("start", "GenoSet", function(x) { return(start(locData(x))) } )

##' Get end of location for each feature
##'
##' Get end of location for each feature
##' @param x GenoSet
##' @return integer
##' @author Peter M. Haverty
##' @aliases end,GenoSet-method
setMethod("end", "GenoSet", function(x) { return(end(locData(x))) } )

##' Get width of location for each feature
##'
##' Get width of location for each feature
##' @param x GenoSet
##' @return integer
##' @author Peter M. Haverty
setMethod("width", "GenoSet", function(x) { return(width(locData(x))) } )

##' Get data matrix names
##'
##' Get names of data matrices. For the time being, this is \code{assayDataElementNames}. This function used to do \code{chrNames}.
##' @param x GenoSet
##' @return character
##' @author Peter Haverty
##' @exportMethod names
##' @aliases names,GenoSet-method
setMethod("names", "GenoSet", function(x) {
  return( assayDataElementNames(x) )
} )

setMethod("ranges", "GenoSet", function(x) {
  .Defunct(new="ranges",package="genoset",msg="The ranges method on a GenoSet is defunct. Please use ranges(locData(x)).")
})

setMethod("space", "GenoSet", function(x) {
  .Defunct(new="space",package="genoset",msg="The ranges method on a GenoSet is defunct. Please use space(locData(x)) or seqnames(locData(x)) as appropriate for RangedData or GRanges.")
} )

##' Get elementLengths from locData slot
##'
##' Get elementLengths from locData slot
##' @title ElementLengths for chromosome
##' @param x GenoSet
##' @return character
##' @author Peter Haverty
##' @exportMethod elementLengths
##' @aliases elementLengths,GenoSet-method
setMethod("elementLengths", "GenoSet", function(x) { return( elementLengths(locData(x)) ) } )
##' @aliases elementLengths,GRanges-method
setMethod("elementLengths", "GRanges", function(x) {
  if ( any(duplicated(runValue(seqnames(x)))) ) {  stop("GRanges not ordered by chromosome.") }
  return( structure(runLength(seqnames(x)),names=as.character(runValue(seqnames(x)))) )
})

##' @exportMethod nrow
##' @aliases nrow,GRanges-method
setMethod("nrow", "GRanges", function(x) { length(x) })

##' @exportMethod dim
##' @aliases dim,GenoSet-method
setMethod("dim", "GenoSet", function(x) { c(nrow(unname(featureData(x))),nrow(unname(phenoData(x))))})

#############
# Sub-setters
#############

##' Subset a GenoSet
##'
##' Subset a GenoSet
##' @exportMethod "["
##' @param x GenoSet
##' @param i character, RangedData, logical, integer
##' @param j character, RangedData, logical, integer
##' @param k character or integer
##' @param drop logical drop levels of space factor?
##' @param ... additional subsetting args
##' @examples
##'   data(genoset)
##'   genoset.ds[1:5,2:3]  # first five probes and samples 2 and 3
##'   genoset.ds[ , "K"]  # Sample called K
##'   rd = RangedData(ranges=IRanges(start=seq(from=15e6,by=1e6,length=7),width=1),names=letters[8:14],space=rep("chr17",7))
##'   genoset.ds[ rd, "K" ]  # sample K and probes overlapping those in rd, which overlap specifed ranges on chr17
##' @aliases [,GenoSet,ANY,ANY,ANY-method
##' @aliases [,GenoSet,ANY-method
##' @aliases [,GenoSet,RangedDataOrGenomicRanges-method
##' @aliases [,GenoSet,character-method
setMethod("[", signature=signature(x="GenoSet",i="ANY",j="ANY"),
          function(x,i,j,k,...,drop=FALSE) {
            if (! missing(k)) {
              if (is.numeric(k)) {
                if (k > length(assayDataElementNames(x))) {
                  stop("Numeric index k exceeds the number of assayDataElements.\n")
                }
                k = assayDataElementNames(x)[k]
              }
              if (!k %in% assayDataElementNames(x)) {
                stop("Index k is not a member of assayDataElementNames.\n")
              }
              if (missing(i) && missing(j)) {
                return(assayDataElement(x,k)) # Necessary to get whole big.matrix object
              } else if (missing(i)) {
                return(assayDataElement(x,k)[,j])
              } else if (missing(j)) {
                return(assayDataElement(x,k)[i,])
              } else {
                return(assayDataElement(x,k)[i,j])
              }
            }
            if ( ! missing(i) ) {
              # Re-ordering of RangedData can silently disobey in order to keep its desired order of chromosomes
              locs = locData(x)[i,,drop=TRUE]
              x@locData = locs
              i = match(rownames(locs),rownames(x))
            }
            callNextMethod(x,i,j,...,drop=drop)
          })

##' @aliases [,GenoSet,character,ANY,ANY-method
setMethod("[", signature=signature(x="GenoSet",i="character",j="ANY"),
          function(x,i,j,...,drop=FALSE) {
            if ( ! missing(i) ) {
              indices = match(i,rownames(x))
            }
            callNextMethod(x,indices,j,...,drop=drop)
          })

##' @aliases [,GenoSet,RangedDataOrGenomicRanges,ANY,ANY-method
setMethod("[", signature=signature(x="GenoSet", i="RangedDataOrGenomicRanges", j="ANY"),
          function(x,i,j,...,drop=FALSE) {
            indices = unlist(x@locData %over% i)
            callNextMethod(x,indices,j,...,drop=drop)
          })

##' @aliases [<-,GenoSet,ANY,ANY,ANY-method
setMethod("[<-", signature=signature(x="GenoSet", i="ANY", j="ANY"),          
          function(x,i,j,k,value) {
            if ( missing(k)) {
              stop("Must specify k to replace data in the GenoSet")
            }
            if (is.numeric(k)) {
                if (k > length(assayDataElementNames(x))) {
                  stop("Numeric index k exceeds the number of assayDataElements.\n")
                }
                k = assayDataElementNames(x)[k]
              }
            if (missing(i) && missing(j)) {
              if (! all( colnames(x) == colnames(value)) || ! all( rownames(x) == rownames(value))) {
                stop("Dimnames for incoming assayDataElement must match this genoset.\n")
              }
              return(assayDataElementReplace(x,k,value))
            }
            if (!k %in% assayDataElementNames(x)) {
              stop("Index k is not a member of assayDataElementNames.\n")
            }
            if (missing(i)) {
              assayDataElement(x,k)[,j] = value
              return(x)
            }
            if (is(i,"RangedData") || is(i,"GRanges")) {
              i = unlist(locData(x) %over% i)
            }
            if (missing(j)) {
              assayDataElement(x,k)[i,] = value
            } else {
              assayDataElement(x,k)[i,j] = value
            }
            return(x)
          })

#######
# Other
#######

##' Print a GenoSet
##'
##' Prints out a description of a GenoSet object
##' @exportMethod show
##' @aliases show,GenoSet-method
setMethod("show","GenoSet",
          function(object) {
            callNextMethod(object)
            cat("Feature Locations:\n")
            show(locData(object))
          })

########################
# Get genome information
########################

##' Chromosome name for each feature
##'
##' Get chromosome name for each feature.  Returns character, not the factor 'space'.
##' @title Look up chromosome for each feature 
##' @param object GRanges, RangedData or GenoSet
##' @return character vector of chromosome positions for each feature
##' @examples
##'   data(genoset)
##'   chr(genoset.ds)  # c("chr1","chr1","chr1","chr1","chr3","chr3","chrX","chrX","chrX","chrX")
##'   chr(locData(genoset.ds))  # The same
##' @author Peter Haverty
##' @export chr
setGeneric("chr", function(object) standardGeneric("chr"))
##' @aliases chr,RangedData-method
setMethod("chr", "RangedData", function(object) { return(as.character(space(object))) } )
##' @aliases chr,GenoSet-method
setMethod("chr", "GenoSet", function(object) { return(chr(slot(object,"locData"))) } )
##' @aliases chr,GRanges-method
setMethod("chr", "GRanges", function(object) { return(as.character(seqnames(object))) })

##' Chromosome position of features
##'
##' Get chromosome position of features/ranges. Defined as floor of mean of start and end.
##' @title Positions for features
##' @param object GRanges, RangedData or GenoSet
##' @return numeric vector of feature positions within a chromosome
##' @author Peter Haverty
##' @export pos
##' @examples
##'   data(genoset)
##'   pos(genoset.ds)  # 1:10
##'   pos(locData(genoset.ds))  # The same
##' @aliases pos,RangedDataOrGenoSetOrGenomicRanges-method
setGeneric("pos", function(object) standardGeneric("pos"))
setMethod("pos", "RangedDataOrGenoSetOrGenomicRanges",
          function(object) { return( start(object) + (width(object) - 1L) %/% 2L) } )

##' Get list of unique chromosome names
##'
##' Get list of unique chromosome names
##' 
##' @param object RangedData or GenoSet
##' @return character vector with names of chromosomes
##' @author Peter M. Haverty
##' @export chrNames
##' @examples
##'   data(genoset)
##'   chrNames(genoset.ds) # c("chr1","chr3","chrX")
##'   chrNames(locData(genoset.ds))  # The same
##'   chrNames(genoset.ds) = sub("^chr","",chrNames(genoset.ds))
setGeneric("chrNames", function(object) standardGeneric("chrNames") )
##' @aliases chrNames,GenoSet-method
setMethod("chrNames", signature(object="GenoSet"),
          function(object) {
            chrNames(locData(object))
          })
##' @aliases chrNames,RangedData-method
setMethod("chrNames", signature(object="RangedData"),
          function(object) {
            names(object)
          })
##' @aliases chrNames,GRanges-method
setMethod("chrNames", signature(object="GRanges"),
          function(object) {
            as.character(unique(seqnames(object)))
          })

##' @export "chrNames<-"
setGeneric("chrNames<-", function(object,value) standardGeneric("chrNames<-") )
##' @aliases chrNames<-,GenoSet-method
setMethod("chrNames<-", signature(object="GenoSet"),
          function(object,value) {
            chrNames(locData(object)) = value
            return(object)
          })
##' @aliases chrNames<-,RangedData-method
setMethod("chrNames<-", signature(object="RangedData"),
          function(object,value) {
            names(object) = value
            return(object)
          })
##' @aliases chrNames<-,GRanges-method
setMethod("chrNames<-", signature(object="GRanges"),
          function(object,value) {
            seqlevels(object) = value
            return(object)
          })

##' Get chromosome start and stop positions
##'
##' Provides a matrix of start, stop and offset, in base numbers for each chromosome.
##' 
##' @title Chromosome Information
##' @param object A GenoSet object or similar
##' @return list with start and stop position, by ordered chr
##' @author Peter Haverty
##' @export chrInfo
##' @examples
##'   data(genoset)
##'   chrInfo(genoset.ds)
##'   chrInfo(locData(genoset.ds))  # The same
setGeneric("chrInfo", function(object) standardGeneric("chrInfo") )
##' @aliases chrInfo,RangedDataOrGenoSetOrGenomicRanges-method
setMethod("chrInfo", signature(object="RangedDataOrGenoSetOrGenomicRanges"),
          function(object) {
            # Get max end value for each chr
            if (class(object) == "GRanges" && !any(is.na(seqlengths(object)))) {
              max.val = seqlengths(object)
            } else {
              chr.ind=chrIndices(object)
              max.val = aggregate(end(object), start=chr.ind[,1], end=chr.ind[,2], FUN=max)
            }
            if (length(max.val) == 1) {
              names(max.val) = chrNames(object)
            } else {
              max.val = max.val[ chrOrder(chrNames(object)) ]
            }

            chr.info = matrix(ncol=3,nrow=length(max.val), dimnames=list(names(max.val),c("start","stop","offset")))
            chr.info[,"stop"]    = cumsum(as.numeric(max.val))
            chr.info[,"offset"]  = c(0, chr.info[- nrow(chr.info),"stop"])
            chr.info[,"start"]   = chr.info[,"offset"] + 1
            
            return(chr.info)
          })


##' Get a matrix of first and last index of features in each chromosome
##'
##' Sometimes it is handy to know the first and last index for each chr.
##' This is like chrInfo but for feature indices rather than chromosome
##' locations. If chr is specified, the function will return a sequence
##' of integers representing the row indices of features on that chromosome.
##' 
##' @param object GenoSet, RangedData, or GRanges
##' @param chr character, specific chromosome name
##' @return data.frame with "first" and "last" columns
##' @author Peter M. Haverty
##' @export chrIndices
##' @examples
##'   data(genoset)
##'   chrIndices(genoset.ds)
##'   chrIndices(locData(genoset.ds))  # The same
setGeneric("chrIndices", function(object,chr=NULL) standardGeneric("chrIndices") )
##' @aliases chrIndices,RangedDataOrGenoSetOrGenomicRanges-method
setMethod("chrIndices", signature(object="RangedDataOrGenoSetOrGenomicRanges"),
          function(object,chr=NULL) {
            object.lengths = elementLengths(object)
            object.lengths = object.lengths[ object.lengths > 0 ]
            chr.last = cumsum(object.lengths)
            chr.last = chr.last[ chr.last > 0 ]
            chr.names = names(chr.last)
            chr.first = c(1,chr.last[- length(chr.last)] + 1)
            chr.info = matrix(c(chr.first,chr.last, chr.first-1), ncol=3,nrow=length(chr.names), dimnames=list(chr.names,c("first","last","offset")))
            if (!is.null(chr)) {
              if (! chr %in% rownames(chr.info)) { stop("Must specify a valid chromosome name in chrIndices.\n") }
              return( seq.int( chr.info[chr,"first"], chr.info[chr,"last"]) )
            } else {
              return(chr.info)
            }
        })

##' Get base positions of features in genome-scale units
##'
##' Get base positions of array features in bases counting from the start of the
##' genome. Chromosomes are ordered numerically, when possible, then lexically.
##'
##' @title Convert chromosome positions to positions from start of genome
##' @param object A GenoSet object or a RangedData object
##' @return numeric position of each feature in whole genome units, in original order
##' @author Peter M. Haverty
##' @examples
##'   data(genoset)
##'   head(genoPos(genoset.ds))
##'   head(genoPos(locData(genoset.ds)))  # The same
##' @export genoPos
setGeneric("genoPos", function(object) standardGeneric("genoPos") )
##' @aliases genoPos,RangedDataOrGenoSetOrGenomicRanges-method
setMethod("genoPos", signature(object="RangedDataOrGenoSetOrGenomicRanges"),
          function(object) {

            # For single chr objects, just return pos
            if ( length(chrNames(object)) == 1 ) {
              return(pos(object))
            }

            ### Add offset to pos by chr
            offset = chrInfo(object)[,"offset"]
            genopos = pos(object) + unlist(offset[chr(object)])
            
            return(genopos)
          })

###########
# Functions
###########

##' Subset or re-order assayData
##'
##' Subset or re-order assayData locked environment, environment, or list. Shamelessly stolen
##' from "[" method in Biobase version 2.8 along with guts of assayDataStorageMode()
##' @title Subset assayData
##' @param orig assayData environment
##' @param i row indices
##' @param j col indices
##' @param ... Additional args to give to subset operator
##' @param drop logical, drop dimensions when subsetting with single value?
##' @return assayData data structure
##' @export subsetAssayData
##' @examples
##'   data(genoset)
##'   ad = assayData(genoset.ds)
##'   small.ad = subsetAssayData(ad,1:5,2:3)
##' @author Peter M. Haverty
subsetAssayData <- function(orig, i, j, ..., drop=FALSE) {
  if (is(orig, "list")) {
    if (missing(i))                     # j must be present
      return(lapply(orig, function(obj) obj[, j, ..., drop = drop]))
    else {                              # j may or may not be present
      if (missing(j))
        return(lapply(orig, function(obj) obj[i,, ..., drop = drop]))
      else
        return(lapply(orig, function(obj) obj[i, j, ..., drop = drop]))
    }
  } else {
    aData <- new.env(parent=emptyenv())
    if (missing(i))                     # j must be present
      for(nm in ls(orig)) aData[[nm]] <- orig[[nm]][, j, ..., drop = drop]
    else {                              # j may or may not be present
      if (missing(j))
        for(nm in ls(orig)) aData[[nm]] <- orig[[nm]][i,, ..., drop = drop]
      else
        for(nm in ls(orig)) aData[[nm]] <- orig[[nm]][i, j, ..., drop = drop]
    }
    if (environmentIsLocked(orig)) {
      lockEnvironment(aData, bindings=TRUE)
    }
    return(aData)
  }
}
