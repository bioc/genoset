####  Class definition for GenoSet, which will extend eSet
######   GenoSet will provide a locData slot containing a RangedData object from the IRanges
######   package to hold genome locations of the features and allow for easy subsetting
######   by location.
######   Intended to be subset by other classes to add one or more data matrices to
######   the assayData slot.

##' GenoSet: An eSet for data with genome locations
##' 
##' Load, manipulate, and plot copynumber and BAF data. GenoSet class
##' extends eSet by adding a "locData" slot for a GRanges or RangedData object.
##' This object contains feature genome location data and
##' provides for efficient subsetting on genome location. CNSet and BAFSet extend
##' GenoSet and require assayData matrices for Copy Number (cn) or Log-R Ratio
##' (lrr) and B-Allele Frequency (baf) data. Implements and provides
##' convenience functions for processing of copy number and B-Allele Frequency
##' data.
##'
##' @docType package
##' @name genoset-package
##' @aliases genoset genoset-package
##' @seealso genoset-datasets GenoSet CNSet BAFSet
##' 
##' @importClassesFrom Biobase AnnotatedDataFrame AssayData eSet ExpressionSet MIAME Versioned VersionedBiobase
##' @importClassesFrom IRanges DataFrame RangedData Rle
##' @importClassesFrom GenomicRanges GRanges
##'
##' @importMethodsFrom GenomicRanges seqnames seqlevels names "names<-" length width
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
##' The CNSet and BAFSet classes have been deprecated.  They only really added getter/setter methods for specific assayDataElements,
##' so they are now redundant with the preferred method of using the assayDataElement name as the third argument to bracket, e.g.
##' \code{x[i, j, "lrr"]}. Accordingly \code{BAFSet.to.ExpressionSets} is also deprecated.
##'
##' Additionally, names, ranges, and space on a GenoSet are also deprecated. In an effort to make a consistent API for either RangedData or
##' GRanges in the locData slot, we recommend using \code{chrNames} for \code{names} and \code{chr} for \code{space}.
##' @name genoset-deprecated
##' @aliases genoset-deprecated
NULL
###############
# Class GenoSet
###############

##' @exportClass GenoSet
setClassUnion("RangedDataOrGRanges",c("RangedData","GRanges"))
setClass("GenoSet", contains=c("eSet"), representation=representation(locData="RangedDataOrGRanges"))
setClassUnion("RangedDataOrGenoSet",c("RangedData","GenoSet"))
setClassUnion("RangedDataOrGenoSetOrGRanges",c("RangedData","GenoSet","GRanges"))

setValidity("GenoSet", function(object) {
  return( all( featureNames(locData(object)) == featureNames(assayData(object)) ) )
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
##' locations. featureNames (names or rownames) are required to match featureNames.
##' @param pData A data frame with rownames matching sampleNames (colnames of all assayDataElements)
##' @param annotation character, string to specify chip/platform type
##' @param universe character, a string to specify the genome universe for locData, overrides universe/genome data in locData
##' @param assayData assayData, usually an environment
##' @param ... More matrix or DataFrame objects to include in assayData
##' @return A GenoSet object or derivative as specified by "type" arg
##' @examples 
##'   test.sample.names = LETTERS[11:13]
##'   probe.names = letters[1:10]
##'   gs = GenoSet(
##'      locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18"),
##'      cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
##'      pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
##'      annotation="SNP6"
##'   )
##' @author Peter M. Haverty
initGenoSet <- function(type, locData, pData=NULL, annotation="", universe, assayData=NULL, ...) {
  # Function to clean up items for slots and call new for GenoSet and its children
  # ... will be the matrices that end up in assayData

  if (! missing(universe)) {
    universe(locData) = universe
  }

  # Check/set genome order of locData
  if ( ! isGenomeOrder(locData, strict=TRUE) ) {
    locData = toGenomeOrder(locData, strict=TRUE )
  }
  clean.loc.rownames = featureNames(locData)

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
  clean.featureNames = featureNames(ad)

  if (length(clean.featureNames) != length(clean.loc.rownames)) {
    stop("Row number mismatch for assayData and locData")
  }

  # Set row order to match locData, already know all ad elements have same row names
  if ( ! all(clean.loc.rownames == clean.featureNames) ) {
    if (! setequal(clean.loc.rownames, clean.featureNames)) {
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
  featureNames(locData) = NULL
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
##'    locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18"),
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

##' Genome universe for locData
##'
##' The genome positions of the features in locData. The UCSC notation (e.g. hg18, hg19, etc.) should be used. For a
##' GRanges, the first value is returned if there are multiple.
##'
##' @title Get and set the genome universe annotation.
##' @param x GenoSet or GRanges
##' @return character, e.g. hg19
##' @author Peter M. Haverty
##' @exportMethod universe
##' @rdname genoset-methods
##' @examples
##'   data(genoset)
##'   universe(locData.rd)
##'   universe(locData.rd) = "hg19"
##' @aliases universe,GenoSet-method
##' @aliases universe,GRanges-method
setMethod("universe", "GenoSet", function(x) { return(universe(x@locData)) } )
setMethod("universe", "GRanges", function(x) {
  if (length(unique(genome(x))) != 1) {
    warning("Taking first element of GRanges genome as universe.")
  }
  return(unname(genome(x)[1]))
} )

##' Set genome universe
##'
##' Set genome universe
##' 
##' @param x GenoSet or GRanges
##' @param value character, new universe string, e.g. hg19
##' @return updated copy of x
##' @author Peter Haverty
##' @exportMethod "universe<-"
##' @rdname genoset-methods
##' @aliases universe<-,GenoSet-method
##' @aliases universe<-,GRanges-method
setMethod("universe<-", signature(x="GenoSet"),
                 function(x,value) {
                   universe(x@locData) = value
                   return(x)
                   })
setMethod("universe<-", signature(x="GRanges"),
          function(x,value) {
            genome(x) = value
            return(x)
          })

##' Get rownames from RangedData, GRanges, or GenoSet
##'
##' Get rownames from RangedData, GRanges, or GenoSet
##' 
##' @param object GRanges, RangedData, or GenoSet
##' @return character vector with names rows/features
##' @author Peter M. Haverty
##' @exportMethod featureNames
##' @examples
##'   data(genoset)
##'   head(featureNames(locData.rd))
##'   head(featureNames(as(locData.rd,"GRanges")))
##'   head(featureNames(cn.ds))
##' @exportMethod featureNames
##' @rdname featureNames
##' @aliases featureNames,GRanges-method
setMethod("featureNames", signature(object="GRanges"),
          function(object) {
            names(object)
          })

##' @rdname featureNames
##' @aliases featureNames,RangedData-method
setMethod("featureNames", signature(object="RangedData"),
          function(object) {
            rownames(object)
          })

##' @rdname featureNames
##' @aliases featureNames,GenoSet-method
setMethod("featureNames", signature(object="GenoSet"),
          function(object) {
            return(unname(featureNames(featureData(object))))
          })

##' Get sampleNames from a GenoSet
##'
##' Get sampleNames from a GenoSet
##' 
##' @param object GenoSet
##' @return character vector with names of samples
##' @exportMethod sampleNames
##' @examples
##'   data(genoset)
##'   head(sampleNames(cn.ds))
##' @exportMethod sampleNames
##' @rdname sampleNames
##' @aliases sampleNames,GenoSet-method
setMethod("sampleNames", signature(object="GenoSet"),
          function(object) {
            rownames(pData(object))
          })

##' Set featureNames
##'
##' Set featureNames of a GenoSet, GRanges, or RangedData (rownames, names, or rownames respectively).
##' @title Set featureNames
##' @param object GenoSet, RangedData, or GRanges
##' @param value ANY
##' @return A new object of the class of supplied object
##' @exportMethod "featureNames<-"
##' @author Peter M. Haverty
##' @rdname featureNames-set
##' @aliases featureNames<-,GenoSet-method
setMethod("featureNames<-",
                 signature=signature(object="GenoSet", value="ANY"),
                 function(object, value) {
                   object = callNextMethod(object,value)
                   featureNames(slot(object,"featureData")) = value
                   return(object)
                 })
##' @rdname featureNames-set
##' @aliases featureNames<-,GRanges-method
setMethod("featureNames<-",
                 signature=signature(object="GRanges", value="ANY"),
                 function(object, value) {
                   names(object) = value
                   return(object)
                 })
##' @rdname featureNames-set
##' @aliases featureNames<-,RangedData-method
setMethod("featureNames<-",
                 signature=signature(object="RangedData", value="ANY"),
                 function(object, value) {
                   rownames(object) = value
                   return(object)
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
##' @aliases locData<-,GenoSet,RangedDataOrGRanges-method
##' @aliases locData
##' @aliases locData<-
##' @return A GenoSet object
##' @rdname locData-methods
setGeneric("locData", function(object) standardGeneric("locData"))
setMethod("locData", "GenoSet", function(object) {
  locs = slot(object,"locData")
  featureNames(locs) = featureNames(object)
  return(locs)
} )
setGeneric("locData<-", function(object,value) standardGeneric("locData<-") )
setMethod("locData<-", signature(object="GenoSet", value="RangedDataOrGRanges"),
                 function(object,value) {
                   if (! all( featureNames(value) %in% featureNames(object))) {
                       stop("Can not replace locData using rownames not in this GenoSet")
                     }
                   if (! all(featureNames(value) == featureNames(object))) {
                     object = object[featureNames(value), ]
                   }
                   featureNames(value) = NULL
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
##' @rdname genoset-methods
##' @aliases start,GenoSet-method
setMethod("start", "GenoSet", function(x) { return(start(locData(x))) } )

##' Get end of location for each feature
##'
##' Get end of location for each feature
##' @param x GenoSet
##' @return integer
##' @author Peter M. Haverty
##' @rdname genoset-methods
##' @aliases end,GenoSet-method
setMethod("end", "GenoSet", function(x) { return(end(locData(x))) } )

##' Get width of location for each feature
##'
##' Get width of location for each feature
##' @param x GenoSet
##' @return integer
##' @author Peter M. Haverty
##' @rdname genoset-methods
setMethod("width", "GenoSet", function(x) { return(width(locData(x))) } )

##' Get data matrix names
##'
##' Get names of data matrices. For the time being, this is \code{assayDataElementNames}. This function used to do \code{chrNames}.
##' @param x GenoSet
##' @return character
##' @author Peter Haverty
##' @exportMethod names
##' @rdname genoset-methods
##' @aliases names,GenoSet-method
setMethod("names", "GenoSet", function(x) {
  return( assayDataElementNames(x) )
} )

##' Get ranges from locData slot
##'
##' Get ranges from locData slot. The ranges method on a GenoSet is deprecated. Please use ranges(locData(x)).
##' @title Ranges for chromosome
##' @param x GenoSet
##' @return character
##' @author Peter Haverty
##' @exportMethod ranges
##' @rdname genoset-methods
##' @aliases ranges,GenoSet-method
setMethod("ranges", "GenoSet", function(x) {
  .Deprecated(old="ranges",package="genoset",msg="The ranges method on a GenoSet is deprecated. Please use ranges(locData(x)).")
  return( ranges(locData(x)) )
})


##' Get space factor for GenoSet
##'
##' locData slot holds a RangedData, which keeps the chromosome of each
##' feature in a factor names 'space'. The ranges method on a GenoSet is deprecated. Please use space(locData(x)) or seqnames(locData(x)) as appropriate for RangedData or GRanges.
##' @param x GenoSet
##' @return factor
##' @author Peter M. Haverty
##' @rdname genoset-methods
##' @examples
##' data(genoset)
##' chr(genoset.ds)
##' start(genoset.ds)
##' end(genoset.ds)
##' chrNames(genoset.ds)
##' elementLengths(genoset.ds) # Returns the number of probes per chromosome
##' @aliases space,GenoSet-method
setMethod("space", "GenoSet", function(x) {
  .Deprecated(old="space",package="genoset",msg="The ranges method on a GenoSet is deprecated. Please use space(locData(x)) or seqnames(locData(x)) as appropriate for RangedData or GRanges.")
  return(space(locData(x)))
} )

##' Get elementLengths from locData slot
##'
##' Get elementLengths from locData slot
##' @title ElementLengths for chromosome
##' @param x GenoSet
##' @return character
##' @author Peter Haverty
##' @exportMethod elementLengths
##' @rdname genoset-methods
##' @aliases elementLengths,GenoSet-method
setMethod("elementLengths", "GenoSet", function(x) { return( elementLengths(locData(x)) ) } )
##' @rdname genoset-methods
##' @aliases elementLengths,GRanges-method
setMethod("elementLengths", "GRanges", function(x) {
  if ( any(duplicated(runValue(seqnames(x)))) ) {  stop("GRanges not ordered by chromosome.") }
  return( structure(runLength(seqnames(x)),names=as.character(runValue(seqnames(x)))) )
})

##' @exportMethod nrow
##' @rdname genoset-methods
##' @aliases nrow,GRanges-method
setMethod("nrow", "GRanges", function(x) { length(x) })

##' @exportMethod dim
##' @rdname genoset-methods
##' @aliases dim,GenoSet-method
setMethod("dim", "GenoSet", function(x) { c(nrow(unname(featureData(x))),nrow(unname(phenoData(x))))})

#############
# Sub-setters
#############

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
##' @rdname genoset-methods
##' @aliases [,GenoSet,ANY,ANY,ANY-method
##' @aliases [,GenoSet,ANY-method
##' @aliases [,GenoSet,RangedDataOrGRanges-method
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
              i = match(featureNames(locs),featureNames(featureData(x)))
            }
            callNextMethod(x,i,j,...,drop=drop)
          })

##' @rdname genoset-methods
##' @aliases [,GenoSet,character,ANY,ANY-method
setMethod("[", signature=signature(x="GenoSet",i="character",j="ANY"),
          function(x,i,j,...,drop=FALSE) {
            if ( ! missing(i) ) {
              indices = match(i,featureNames(x))
            }
            callNextMethod(x,indices,j,...,drop=drop)
          })

##' @rdname genoset-methods
##' @aliases [,GenoSet,RangedDataOrGRanges,ANY,ANY-method
setMethod("[", signature=signature(x="GenoSet", i="RangedDataOrGRanges", j="ANY"),
          function(x,i,j,...,drop=FALSE) {
            indices = unlist(x@locData %over% i)
            callNextMethod(x,indices,j,...,drop=drop)
          })

##' @rdname genoset-methods
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
              if (! all( sampleNames(x) == colnames(value)) || ! all( featureNames(x) == rownames(value))) {
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

##' @exportMethod show
##' @rdname genoset-methods
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
##'   test.sample.names = LETTERS[11:13]
##'   probe.names = letters[1:10]
##'   gs = GenoSet(
##'      locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18"),
##'      cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
##'      pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
##'      annotation="SNP6"
##'   )
##'   chr(gs)  # c("chr1","chr1","chr1","chr1","chr3","chr3","chrX","chrX","chrX","chrX")
##'   chr(locData(gs))  # The same
##' @author Peter Haverty
##' @export chr
##' @rdname chr-methods
setGeneric("chr", function(object) standardGeneric("chr"))
##' @rdname chr-methods
##' @aliases chr,RangedData-method
setMethod("chr", "RangedData", function(object) { return(as.character(space(object))) } )
##' @rdname chr-methods
##' @aliases chr,GenoSet-method
setMethod("chr", "GenoSet", function(object) { return(chr(slot(object,"locData"))) } )
##' @rdname chr-methods
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
##'   test.sample.names = LETTERS[11:13]
##'   probe.names = letters[1:10]
##'   gs = GenoSet(
##'      locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18"),
##'      cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
##'      pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
##'      annotation="SNP6"
##'   )
##'   pos(gs)  # 1:10
##'   pos(locData(gs))  # The same
##' @rdname pos
##' @aliases pos,RangedDataOrGenoSetOrGRanges-method
setGeneric("pos", function(object) standardGeneric("pos"))
##' @rdname pos
setMethod("pos", "RangedDataOrGenoSetOrGRanges",
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
##'   test.sample.names = LETTERS[11:13]
##'   probe.names = letters[1:10]
##'   gs = GenoSet(
##'      locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18"),
##'      cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
##'      pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
##'      annotation="SNP6"
##'   )
##'   chrNames(gs) # c("chr1","chr3","chrX")
##'   chrNames(locData(gs))  # The same
##'   chrNames(gs) = sub("^chr","",chrNames(gs))
##' @rdname chrNames
setGeneric("chrNames", function(object) standardGeneric("chrNames") )
##' @rdname chrNames
##' @aliases chrNames,GenoSet-method
setMethod("chrNames", signature(object="GenoSet"),
          function(object) {
            chrNames(locData(object))
          })
##' @rdname chrNames
##' @aliases chrNames,RangedData-method
setMethod("chrNames", signature(object="RangedData"),
          function(object) {
            names(object)
          })
##' @rdname chrNames
##' @aliases chrNames,GRanges-method
setMethod("chrNames", signature(object="GRanges"),
          function(object) {
            as.character(unique(seqnames(object)))
          })

##' @rdname chrNames
##' @export "chrNames<-"
setGeneric("chrNames<-", function(object,value) standardGeneric("chrNames<-") )
##' @rdname chrNames
##' @aliases chrNames<-,GenoSet-method
setMethod("chrNames<-", signature(object="GenoSet"),
          function(object,value) {
            chrNames(locData(object)) = value
            return(object)
          })
##' @rdname chrNames
##' @aliases chrNames<-,RangedData-method
setMethod("chrNames<-", signature(object="RangedData"),
          function(object,value) {
            names(object) = value
            return(object)
          })
##' @rdname chrNames
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
##' @rdname chrInfo
setGeneric("chrInfo", function(object) standardGeneric("chrInfo") )
##' @rdname chrInfo
##' @aliases chrInfo,RangedDataOrGenoSetOrGRanges-method
setMethod("chrInfo", signature(object="RangedDataOrGenoSetOrGRanges"),
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
##' @rdname chrIndices-methods
setGeneric("chrIndices", function(object,chr=NULL) standardGeneric("chrIndices") )
##' @rdname chrIndices-methods
##' @aliases chrIndices,RangedDataOrGenoSetOrGRanges-method
setMethod("chrIndices", signature(object="RangedDataOrGenoSetOrGRanges"),
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
##' @rdname genoPos-methods
setGeneric("genoPos", function(object) standardGeneric("genoPos") )
##' @rdname genoPos-methods
##' @aliases genoPos,RangedDataOrGenoSetOrGRanges-method
setMethod("genoPos", signature(object="RangedDataOrGenoSetOrGRanges"),
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
## Plots ##
###########

##' Plot data along the genome
##'
##' Plot location data and chromosome boundaries from a GenoSet, RangedData, or GRanges object
##' against data from a numeric or Rle. Specifying a chromosome name and optionally a 'xlim'
##' will zoom into one chromosome region. If more than one chromosome is present, the
##' chromosome boundaries will be marked. Alternatively, for a numeric x and a
##' numeric or Rle y, data in y can be plotted at genome positions x. In this case,
##' chromosome boundaries can be taken from the argument locs. If data for y-axis comes
##' from a Rle lines are plotted representing segments. X-axis tickmarks will be labeled
##' with genome positions in the most appropriate units.
##'
##' @section Methods:
##' \describe{
##' 
##' \item{\code{signature(x = "RangedDataOrGenoSetOrGRanges", y = "ANY")}}{
##' Plot feature locations and data from one sample.
##' }
##' 
##' \item{\code{signature(x = "numeric", y = "numeric")}}{
##' Plot numeric location and a vector of numeric data.
##' }
##' 
##' \item{\code{signature(x = "numeric", y = "Rle")}}{
##' Plot numeric location and a vector of Rle data. Uses lines for Rle runs.
##' }
##' }
##' 
##' @param x GenoSet (or descendant), RangedData, or GRanges
##' @param y numeric or Rle
##' @param element character, Deprecated. when x is a GenoSet, the y-th column of this assayDataElement is used for the y-axis data.
##' @param locs RangedData, like locData slot of GenoSet
##' @param chr Chromosome to plot, NULL by default for full genome
##' @param add Add plot to existing plot
##' @param xlab character, label for x-axis of plot
##' @param ylab character, label for y-axis of plot
##' @param col character, color to plot lines or points
##' @param lwd numeric, line width for segment plots from an Rle
##' @param pch character or numeric, printing character, see points
##' @param xlim integer, length two, bounds for genome positions. Used in conjunction with "chr" to subset data for plotting.
##' @param ... Additional plotting args
##' @return nothing
##' @author Peter M. Haverty
##' @export genoPlot
##' @family "genome plots"
##' @examples
##' data(genoset)
##' genoPlot( x=baf.ds,y=baf.ds[,1,"lrr"] )
##' genoPlot( genoPos(baf.ds), baf.ds[,1,"lrr"], locs=locData(baf.ds) ) # The same
##' genoPlot( 1:10, Rle(c(rep(0,5),rep(3,4),rep(1,1))) )
##' @docType methods
##' @rdname genoPlot-methods
##' @aliases genoPlot-methods
setGeneric("genoPlot", function(x,y,...) { standardGeneric("genoPlot") } )

##' @rdname genoPlot-methods
##' @aliases genoPlot,numeric,numeric-method
setMethod("genoPlot",c(x="numeric",y="numeric"),
          function(x, y, add=FALSE, xlab="", ylab="", col="black", locs=NULL, ...) {
            if (add == FALSE) {
              plot(x,y,axes=FALSE,xlab=xlab,ylab=ylab,xaxs="i",col=col,...)
              genomeAxis(locs=locs)
            } else {
              points(x,y,col=col,...)
            }
            return(invisible())
          })

##' @rdname genoPlot-methods
##' @aliases genoPlot,numeric,Rle-method
setMethod("genoPlot", c(x="numeric",y="Rle"),
          function(x, y, add=FALSE, xlab="", ylab="", col="red", locs=NULL, lwd=2, xlim=NULL, ...) {
            if (add == FALSE) {
              if (is.null(xlim)) {
                xlim=range(x,na.rm=TRUE)
              }
              plot(NA,type="n",xlim=xlim,ylim=range(y,na.rm=TRUE),xlab=xlab,ylab=ylab,xaxs="i",axes=FALSE,...)
              genomeAxis(locs=locs)
            }
            num.mark = runLength(y)
            loc.end.indices = cumsum(num.mark)
            loc.end = x[loc.end.indices]
            loc.start.indices = (loc.end.indices - num.mark) + 1
            loc.start = x[loc.start.indices]
            seg.mean = runValue(y)
            segments(loc.start, seg.mean, loc.end, seg.mean, col=col, lwd=lwd)
            return(invisible())
          })

##' @rdname genoPlot-methods
##' @aliases genoPlot,RangedDataOrGenoSetOrGRanges,ANY-method
setMethod("genoPlot", signature(x="RangedDataOrGenoSetOrGRanges",y="ANY"), function(x, y, element=NULL, chr=NULL, add=FALSE, pch=".", xlab="", ylab="", ...) {
  ## Note: zoom in by subset is much faster (10X) than xlim, so implement a zoom in with subsetting
  # Get position info, subset by chr if necessary
  if (!is.null(element)) {
    warning("The 'element' arg is deprecated. Please switch to genoPlot( RangedDataOrGenoSetOrGRanges, Rle or numeric ), e.g. genoPlot(genoset, genoset[,i,'cn']) .\n")
    if (! element %in% assayDataElementNames(x)) {
      stop("Provided assayData element, ", element, " is not a valid name of an assayData member")
    }
    y = x[ , y, element]
  }
  dot.args = list(...)
  if ( !is.null(chr) ) {
    if ( "xlim" %in% names(dot.args) ) {
      xlim = dot.args[["xlim"]]
      zoom.gr = GRanges(ranges=IRanges(start=xlim[1],end=xlim[2]),seqnames=chr)
      bounds = boundingIndicesByChr(zoom.gr, x)[1,]
      indices = bounds[1]:bounds[2]
      positions = start(x)[indices]
    } else {
      indices = chrIndices(x,chr)
      positions = start(x)[indices]
    }
    y = y[indices]
    locs = NULL
  } else {
    positions = genoPos(x)
    if (length(chrNames(x)) > 1) {
      locs = x
    } else {
      locs = NULL
    }
  }
  genoPlot(positions,y,locs=locs,add=add,xlab=xlab,ylab=ylab,pch=pch,...)
})

##' Label an axis with base positions
##'
##' Label a plot with Mb, kb, bp as appropriate, using tick locations from axTicks
##'
##' @title Label axis with base pair units
##' @param locs RangedData to be used to draw chromosome boundaries, if necessary.  Usually locData slot from a GenoSet.
##' @param side integer side of plot to put axis
##' @param log logical Is axis logged?
##' @param do.other.side logical, label non-genome side with data values at tick marks?
##' @return nothing
##' @export genomeAxis
##' @family "genome plots"
##' @examples
##'   data(genoset)
##'   genoPlot(genoPos(baf.ds), baf(baf.ds)[,1])
##'   genomeAxis( locs=locData(baf.ds) )  # Add chromosome names and boundaries to a plot assuming genome along x-axis
##'   genomeAxis( locs=locData(baf.ds), do.other.side=FALSE ) # As above, but do not label y-axis with data values at tickmarks
##'   genomeAxis()           # Add nucleotide position in sensible units assuming genome along x-axis
##' @author Peter M. Haverty
genomeAxis <- function(locs=NULL, side=1, log=FALSE, do.other.side=TRUE) {
  if (is.null(locs)) {
    label.positions = axTicks(side=side,log=log)
    if ( max(label.positions) > 1e9 ) {
      axis(side=side,at=label.positions,labels=sapply(label.positions,function(x){paste(sprintf("%.1f",x/1e9),"Gb",sep="")}))
    } else if ( max(label.positions) > 1e6 ) {
      axis(side=side,at=label.positions,labels=sapply(label.positions,function(x){paste(sprintf("%.1f",x/1e6),"Mb",sep="")}))
    } else if ( max(label.positions) > 1e3 ) {
      axis(side=side,at=label.positions,labels=sapply(label.positions,function(x){paste(sprintf("%.1f",x/1e3),"kb",sep="")}))
    } else {
      axis(side=side,at=label.positions,labels=sapply(label.positions,function(x){paste(sprintf("%.0f",x),"bp",sep="")}))
    }
  } else {
    chr.info = chrInfo(locs)
    abline(v=chr.info[-1,"start"])
    chr.labels = rownames(chr.info)
    mtext(side=rep(c(3,1),len=length(chr.labels)), text=chr.labels, line=0, at=rowMeans(chr.info[,c("start","stop"),drop=FALSE]))
  }
  box()
  if (do.other.side == TRUE) {
    if (side == 1) {
      axis(side=2)
    } else if (side == 2) {
      axis(side=1)
    }
  }
}

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
  if ("na.action" %in% attributes(ds)) {
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

##' Center continuous data on mode
##'
##' Copynumber data distributions are generally multi-modal. It is often assumed that
##' the tallest peak represents "normal" and should therefore be centered on a
##' log2ratio of zero. This function uses the density function to find the mode of
##' the dominant peak and subtracts that value from the input data.
##' 
##' @param ds numeric matrix
##' @return numeric matrix
##' @author Peter M. Haverty
##' @export
##' @examples
##'   modeCenter( matrix( rnorm(150, mean=0), ncol=3 ))
modeCenter <- function(ds) {
  if (!requireNamespace("stats",quietly=TRUE)) {
    stop("Failed to require stats package.\n")
  }
  column.modes = apply(ds,2, function(x) { 
    l2r.density = stats::density(x,na.rm=TRUE)
    density.max.index = which.max(l2r.density$y)
    return(l2r.density$x[density.max.index])
  })
  ds = sweep(ds, 2, column.modes)
  return(ds)
}

##' Convert bounding indices into a Rle
##'
##' Given a matrix of first/last indices, like from boundingIndicesByChr, and values for
##' each range, convert to a Rle.  This function takes the expected length of the Rle, n,
##' so that any portion of the full length not covered by a first/last range will be a
##' run with the value NA.  This is typical in the case where data is segmented with CBS
##' and some of the data to be segmented is NA.
##' @export 
##' @param bounds matrix, two columns, with first and last index, like from boundingIndicesByChr
##' @param values ANY, some value to be associated with each range, like segmented copy number.
##' @param n integer, the expected length of the Rle, i.e. the number of features in the
##' genome/target ranges processed by boundingIndicesByChr.
##' @return Rle
##' @family "segmented data"
##' @author Peter M. Haverty
bounds2Rle <- function( bounds, values, n ) {
  if ( length(values) != nrow(bounds) ) {
    stop("must have one value for each bound")
  }
  if (n < length(values)) {
    stop("n must be >= length(values)")
  }
  run.length = integer((2*nrow(bounds))+1)
  run.length[1] = bounds[1] - 1
  run.value = rep(NA_real_, (2*nrow(bounds))+1)
  data.indices = seq.int(from=2, by=2, length.out=nrow(bounds))
  widths = (bounds[, 2] - bounds[, 1]) + 1
  run.length[data.indices] = widths
  run.value[data.indices] = values
  run.length[data.indices+1] = diff(c(bounds[, 2], n)) - c(widths[-1], 0)
  
  if (sum(run.length) != n) {
    stop("Rle is the wrong length. Look for double counted features in your bounds table.")
  }
  return( Rle( run.value, run.length ) )
}

##' Find indices of features bounding a set of chromosome ranges/genes
##'
##' This function is similar to findOverlaps but it guarantees at least two features will be
##' covered. This is useful in the case of finding features corresponding to a set of genes.
##' Some genes will fall entirely between two features and thus would not return any ranges
##' with findOverlaps. Specifically, this function will find the indices of the features
##' (first and last) bounding the ends of a range/gene (start and stop) such that
##' first <= start <= stop <= last. Equality is necessary so that multiple conversions between
##' indices and genomic positions will not expand with each conversion. This function uses
##' findIntervals, which is for k queries and n features is O(k * log(n)) generally and
##' ~O(k) for sorted queries. Therefore will be dramatically faster for sets of query genes
##' that are sorted by start position within each chromosome.  This should give performance
##' for k genes and n features that is ~O(k) for starts and O(k * log(n)) for stops and
##' ~O(k * log(n)) overall.  Ranges/genes that are outside the range of feature positions will
##' be given the indices of the corresponding first or last index rather than 0 or n + 1 so
##' that genes can always be connected to some data.
##'
##' @param starts numeric or integer, first base position of each query range
##' @param stops numeric or integer, last base position of each query range
##' @param positions Base positions in which to search
##' @param offset integer, value to add to all returned indices. For the case where positions represents a portion of some larger array (e.g. a chr in a genome)
##' @return integer matrix of 2 columms for start and stop index of range in data
##' @family "range summaries"
##' @export boundingIndices2
##' @examples
##'   starts = seq(10,100,10)
##'   boundingIndices2( starts=starts, stops=starts+5, positions = 1:100 )
##' @author Peter M. Haverty
boundingIndices2 <- function(starts, stops, positions, offset=NULL) {
  indices = c(findInterval(starts,positions,rightmost.closed=FALSE), findInterval(stops,positions,rightmost.closed=FALSE))
  dim(indices) = c(length(starts),2)
  indices[ indices == 0L ] = 1L  # If off left end, set to 1

  right.bounds.to.expand = positions[indices[,2]] < stops
  indices[ right.bounds.to.expand,2 ] = indices[ right.bounds.to.expand,2 ] + 1L  # Right end index moves right one unless it is an exact match
  indices[ indices > length(positions) ] = length(positions)  # If off right end, set to right end
  
  if (!is.null(offset)) {  # Convert indices back to indices in full positions before subsetting by initial.bounds
    indices = indices + as.integer(offset)
  }
  return(indices)
}


##' Find indices of features bounding a set of chromosome ranges/genes
##'
##' This function is similar to findOverlaps but it guarantees at least two features will be
##' covered. This is useful in the case of finding features corresponding to a set of genes.
##' Some genes will fall entirely between two features and thus would not return any ranges
##' with findOverlaps. Specifically, this function will find the indices of the features
##' (first and last) bounding the ends of a range/gene (start and stop) such that
##' first <= start < stop <= last. Equality is necessary so that multiple conversions between
##' indices and genomic positions will not expand with each conversion. Ranges/genes that are
##' outside the range of feature positions will be given the indices of the corresponding
##' first or last index rather than 0 or n + 1 so that genes can always be connected to some data.
##'
##' This function uses some tricks from findIntervals, where is for k queries and n features it
##' is O(k * log(n)) generally and ~O(k) for sorted queries. Therefore will be dramatically
##' faster for sets of query genes that are sorted by start position within each chromosome.
##' The index of the stop position for each gene is found using the left bound from the start
##' of the gene reducing the search space for the stop position somewhat. This function has
##' important differences from boundingIndices2, which uses findInterval: boundingIndices does not
##' check for NAs or unsorted data in the subject positions. Also, the positions are
##' kept as integer, where boundingIndices2 (and findInterval) convert them to doubles.
##' These assumptions are safe for position info coming from a GenoSet, GRanges, or RangedData.
##'
##' @param starts integer vector of first base position of each query range
##' @param stops integer vector of last base position of each query range
##' @param positions Base positions in which to search
##' @param valid.indices logical, TRUE assures that the returned indices don't go off either end of the array, i.e. 0 becomes 1 and n+1 becomes n
##' @param offset integer, value to add to all returned indices. For the case where positions represents a portion of some larger array (e.g. a chr in a genome)
##' @param all.indices logical, return a list containing full sequence of indices for each query
##' @return integer matrix of 2 columms for start and stop index of range in data or a list of full sequences of indices for each query (see all.indices argument)
##' @family "range summaries"
##' @export boundingIndices
##' @examples
##'   starts = seq(10,100,10)
##'   boundingIndices( starts=starts, stops=starts+5, positions = 1:100 )
##' @author Peter M. Haverty \email{phaverty@@gene.com}
boundingIndices <- function(starts,stops,positions,valid.indices=TRUE,all.indices=FALSE, offset=0) {
  bounds = vector("integer",length(starts)*2L)
  bound.results = .C("binary_bound", as.integer(starts), as.integer(stops), as.integer(positions),
    as.integer(length(starts)), as.integer(length(positions)), bounds=bounds, as.integer(valid.indices), as.integer(offset),
    DUP=FALSE, NAOK=TRUE)
  bounds = bound.results$bounds
  dim(bounds) = c(length(starts),2)

  if (all.indices == TRUE) { # Return all covered and bounding indices
    return( apply( bounds, 1, function(x) { seq(from=x[1], to=x[2]) }) )
  } else {  # Just return left and right indices
    return(bounds)
  }
  
}

##' Find indices of features bounding a set of chromosome ranges/genes, across chromosomes
##'
##' Finds subject ranges corresponding to a set of genes (query ranges), taking chromosome
##' into account. Specifically, this function will find the indices of the features
##' (first and last) bounding the ends of a range/gene (start and stop) such that
##' first <= start < stop <= last. Equality is necessary so that multiple conversions between
##' indices and genomic positions will not expand with each conversion. Ranges/genes that are
##' outside the range of feature positions will be given the indices of the corresponding
##' first or last index on that chromosome, rather than 0 or n + 1 so that genes can always be
##' connected to some data. Checking the left and right bound for equality will tell you when
##' a query is off the end of a chromosome.
##' 
##' This function uses some tricks from findIntervals, where is for k queries and n features it
##' is O(k * log(n)) generally and ~O(k) for sorted queries. Therefore will be dramatically
##' faster for sets of query genes that are sorted by start position within each chromosome.
##' The index of the stop position for each gene is found using the left bound from the start
##' of the gene reducing the search space for the stop position somewhat.
##' 
##' This function differs from boundingIndices in that 1. it uses both start and end positions for
##' the subject, and 2. query and subject start and end positions are processed in blocks corresponding
##' to chromosomes.
##'
##' Both query and subject must be in at least weak genome order (sorted by start within chromosome blocks).
##' 
##' @param query GRanges or something coercible to GRanges
##' @param subject RangedData
##' @return integer matrix with two columns corresponding to indices on left and right bound of queries in subject
##' @export boundingIndicesByChr
##' @family "range summaries"
##' @author Peter M. Haverty \email{phaverty@@gene.com}
boundingIndicesByChr <-function(query, subject) {
  # TODO: make the whole thing work for GRanges, do coersion of GenoSet or RangedData to GRanges
  # Convert GRanges to RangedData until genome order, chr, and chrIndices ready for Granges
  if (!is(query,"GRanges")) {
    tryCatch({ query = as(query,"GRanges"); }, error=function(e) { stop("Could not convert query into GRanges.\n") })
  }

  # Subject must have features ordered by start within chromosome. Query need not really, but it's faster.  Just checking query genome order to assure data are in blocks by chromosome in a GRanges. Chromosome order doesn't matter.
  if (! isGenomeOrder(subject,strict=FALSE) ) {
    stop("subject must be in genome order.\n")
  }
  if (! isGenomeOrder(query,strict=FALSE) ) {
    stop("query must be in genome order.\n")
  }
  query.chr.indices = chrIndices(query)
  subject.chr.indices = chrIndices(subject)
  ok.chrs = intersect(rownames(subject.chr.indices),rownames(query.chr.indices))
  query.chr.indices = query.chr.indices[ok.chrs,,drop=FALSE]
  subject.chr.indices = subject.chr.indices[ok.chrs,,drop=FALSE]
  nquery = as.integer(sum(query.chr.indices[,2] - query.chr.indices[,3])) # !!!
  query.start = start(query)
  query.end = end(query)
  query.names = names(query)
  if (is.null(query.names)) { query.names = as.character(seq.int(from=1,to=nquery)) }
  subject.start = start(subject)
  subject.end = end(subject)
  return(.Call("binary_bound_by_chr", nquery, query.chr.indices, query.start, query.end, query.names, subject.chr.indices, subject.start, subject.end))
}

##' Average features in ranges per sample 
##'
##' This function takes per-feature genomic data and returns averages for each of a set of genomic ranges.
##' The most obvious application is determining the copy number of a set of genes. The features
##' corresponding to each gene are determined with boundingIndices such that all features with the bounds
##' of a gene (overlaps). The features on either side of the gene unless those positions
##' exactly match the first or last base covered by the gene.  Therefore, genes falling between two features
##' will at least cover two features. This is similar to rangeSampleMeans, but it checks the subject
##' positions for being sorted and not being NA and also treats them as doubles, not ints. Range bounding
##' performed by the boundingIndices function.
##' 
##' @param query.rd RangedData object representing genomic regions (genes) to be averaged.
##' @param subject A GenoSet object or derivative
##' @param assay.element character, name of element in assayData to use to extract data
##' @return numeric matrix of features in each range averaged by sample
##' @family "range summaries"
##' @export rangeSampleMeans
##' @examples
##'   data(genoset)
##'   my.genes = RangedData( ranges=IRanges(start=c(35e6,128e6),end=c(37e6,129e6),names=c("HER2","CMYC")), space=c("chr17","chr8"), universe="hg19")
##'   rangeSampleMeans( my.genes, baf.ds, "lrr" )
##' @author Peter M. Haverty
rangeSampleMeans <- function(query.rd, subject, assay.element) {
  ## Find feature bounds of each query in subject genoset, get feature data average for each sample
  all.indices = boundingIndicesByChr(query.rd, subject)

  # Temporary hack for DataFrame of Rle
  data.matrix = assayDataElement(subject,assay.element)

  if (class(data.matrix) == "DataFrame") {
    sample.vals = lapply( data.matrix, function(x) { rangeColMeans( all.indices, as.numeric(x)) })
    range.means = do.call(cbind,sample.vals)
  } else if (is.matrix(data.matrix)) {
    range.means = rangeColMeans( all.indices, data.matrix )
  } else {
    range.means = matrix(ncol=ncol(data.matrix),nrow=nrow(all.indices),dimnames=list(rownames(all.indices),colnames(data.matrix)))
    for (i in seq.int(length.out=ncol(data.matrix))) {
      range.means[,i] = rangeColMeans( all.indices, data.matrix[,i] )
    }
  }
  return(range.means)
}

##' Calculate column means for multiple ranges
##'
##' Essentially colMeans with a loop, all in a .Call. Designed to take a
##' 2-column matrix of row indices, bounds, for a matrix, x, and calculate
##' mean for each range in each column (or along a single vector). bounds
##' matrix need not cover all rows.
##' 
##' @param bounds A two column integer matrix of row indices
##' @param x A numeric matrix with rows corresponding to indices in bounds.
##' @return A numeric matrix or vector, matching the form of x. One row for
##' each row in bounds, one col for each col of x and appropriate dimnames.
##' If x is a vector, just a vector with names from the rownames of bounds.
##' @export
##' @family "range summaries"
##' @author Peter M. Haverty \email{phaverty@@gene.com}
rangeColMeans <- function( bounds, x ) {
  if (! is.matrix(bounds) && ncol(bounds) == 2) {
    stop("bounds must be a matrix with 2 columns\n")
  }
  if (!is.double(x)) {
    storage.mode(x) = "double"
  }
  if (!is.integer(bounds)) {
    storage.mode(bounds) = "integer"
  }
  ans = .Call("rangeColMeans", bounds, x)
  return(ans)
}

##' Load a GenoSet from a RData file
##'
##' Given a rds file or a rda file with one object (a GenoSet or related object), load it,
##' and return.
##' @param path character, path to rds or rda file
##' @return GenoSet or related object (only object in RData file)
##' @examples
##' \dontrun{ ds = readGenoSet("/path/to/genoset.RData") }
##' \dontrun{ ds = readGenoSet("/path/to/genoset.rda") }
##' \dontrun{ ds = readGenoSet("/path/to/genoset.rds") }
##' @export 
##' @author Peter M. Haverty \email{phaverty@@gene.com}
readGenoSet <- function(path) {
  header = readLines(path, 1)
  if (grepl("^RD", header)[1] == TRUE) {
    object = get(load(path)[1])    
  } else {
    object = readRDS(path)
  }
  if (!is(object,"eSet")) { stop("Loaded object is not an eSet or derived class.") }
  return( object )
}
