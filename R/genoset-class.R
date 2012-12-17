
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
##' @importMethodsFrom GenomicRanges seqnames seqlevels names "names<-" length
##' @importMethodsFrom Biobase annotation experimentData exprs fData featureNames "featureNames<-" phenoData sampleNames "sampleNames<-"
##' @importMethodsFrom IRanges as.data.frame as.list as.matrix cbind colnames "colnames<-" elementLengths end findOverlaps gsub
##' @importMethodsFrom IRanges "%in%" intersect is.unsorted lapply levels match mean na.exclude nrow order paste ranges Rle rownames
##' @importMethodsFrom IRanges "rownames<-" runLength runValue sapply space start unlist universe "universe<-"
##'
##' @importFrom Biobase assayDataElement assayDataElementNames assayDataElementReplace assayDataNew annotatedDataFrameFrom
##' @importFrom graphics abline axis axTicks box mtext plot.new plot.window points segments
##' @importFrom IRanges DataFrame IRanges RangedData
##' @importFrom GenomicRanges seqlengths GRanges
##'
##' @import methods
##' @import BiocGenerics
##' 
##' @useDynLib genoset
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
            return(unname(featureNames(locData(object))))
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
                   featureNames(slot(object,"locData")) = value
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
setMethod("locData", "GenoSet", function(object) { return(slot(object,"locData")) } )
setGeneric("locData<-", function(object,value) standardGeneric("locData<-") )
setMethod("locData<-", signature(object="GenoSet", value="RangedDataOrGRanges"),
                 function(object,value) {
                   if (! all( featureNames(value) %in% featureNames(object))) {
                       stop("Can not replace locData using rownames not in this GenoSet")
                     }
                   slot(object,"locData") = value
                   if (! all(featureNames(value) == featureNames(object))) {
                     for (adname in assayDataElementNames(object)) {
                       assayDataElement(object,adname) = assayDataElement(object,adname)[featureNames(value),]
                     }
                   }
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

##' Get chromosome names
##'
##' Get chromosome names, which are the names of the locData slot. The names method on a GenoSet is depricated. Please use chrNames.
##' @title Names for chromosome
##' @param x GenoSet
##' @return character
##' @author Peter Haverty
##' @exportMethod names
##' @rdname genoset-methods
##' @aliases names,GenoSet-method
setMethod("names", "GenoSet", function(x) {
  .Deprecated(old="names",new="chrNames",package="genoset",msg="The names method on a GenoSet is depricated. Please use chrNames.")
  return( chrNames(locData(x)) )
} )

##' Get ranges from locData slot
##'
##' Get ranges from locData slot. The ranges method on a GenoSet is depricated. Please use ranges(locData(x)).
##' @title Ranges for chromosome
##' @param x GenoSet
##' @return character
##' @author Peter Haverty
##' @exportMethod ranges
##' @rdname genoset-methods
##' @aliases ranges,GenoSet-method
setMethod("ranges", "GenoSet", function(x) {
  .Deprecated(old="ranges",package="genoset",msg="The ranges method on a GenoSet is depricated. Please use ranges(locData(x)).")
  return( ranges(locData(x)) )
})


##' Get space factor for GenoSet
##'
##' locData slot holds a RangedData, which keeps the chromosome of each
##' feature in a factor names 'space'. The ranges method on a GenoSet is depricated. Please use space(locData(x)) or seqnames(locData(x)) as appropriate for RangedData or GRanges.
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
  .Deprecated(old="space",package="genoset",msg="The ranges method on a GenoSet is depricated. Please use space(locData(x)) or seqnames(locData(x)) as appropriate for RangedData or GRanges.")
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
setMethod("dim", "GenoSet", function(x) { c(nrow(locData(x)),nrow(pData(x))) })

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
              if (is.numeric(k)) { k = assayDataElementNames(x)[k] }
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
              x@locData = x@locData[i,,drop=TRUE]
              i = match(featureNames(locData(x)),featureNames(assayData(x))) # Re-ordering of RangedData can silently disobey in order to keep its desired order of chromosomes
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
            indices = unlist(x@locData %in% i)
            callNextMethod(x,indices,j,...,drop=drop)
          })

##' @rdname genoset-methods
##' @aliases [<-,GenoSet,ANY,ANY,ANY-method
setMethod("[<-", signature=signature(x="GenoSet", i="ANY", j="ANY"),          
          function(x,i,j,k,value) {
            if ( missing(k)) {
              stop("Must specify k to replace data in the GenoSet")
            }
            if (is.numeric(k)) { k = assayDataElementNames(x)[k] }
            if (missing(i) && missing(j)) {
              if (! all( sampleNames(x) == colnames(value)) || ! all( featureNames(x) == rownames(value))) {
                stop("Dimnames for incoming assayDataElement must match this genoset.\n")
              }
              return(assayDataElementReplace(x,k,value))
            }
            if (missing(i)) {
              assayDataElement(x,k)[,j] = value
              return(x)
            }
            if (is(i,"RangedData") || is(i,"GRanges")) {
              i = unlist(locData(x) %in% i)
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
            show(slot(object,"locData"))
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
            seqlevels(object)
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
  fit = stats::lm( ds ~ gc, na.action=na.exclude )
  ds.fixed = stats::residuals(fit)
  if (retain.mean == TRUE) {
    if (is.null(dim(ds))) {
      ds.fixed = ds.fixed + mean(ds,na.rm=TRUE)
    } else {
      ds.fixed = sweep(ds.fixed,2,colMeans(ds,na.rm=TRUE),FUN="+")
    }
  }
  if (is.null(dim(ds))) {
    ds.fixed = unname(ds.fixed)
  } else {
    dimnames(ds.fixed) = dimnames(ds)
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
  cap = bounds[1,1] - 1
  tail = n - bounds[nrow(bounds),2]
  extras = (cap > 0) + (tail > 0) + sum(bounds[-1,1] - bounds[-nrow(bounds),2] > 1)
  if (extras == 0) {
    rle = Rle( values, (bounds[,2] - bounds[,1])+1 )
    if (length(rle) != n) {
      stop("Rle is the wrong length. Look for double counted features in your bounds table.")
    }
    return(rle)
  }
  extras = extras + length(values)
  # Maybe make them 2x + 1 initialized to 0 and NA, let Rle dump zeros that don't get replaced
  run.length = rep(0L,extras)
  run.value = rep(NA_real_,extras)
  if (tail > 0) {
    run.length[ length(run.length) ] = tail
  }
  widths = (bounds[,2] - bounds[,1]) + 1
  if (cap > 0) {
    run.length[1] = cap
    run.value[2] = values[1]
    run.length[2] = widths[1]
    i = 3
  } else {
    run.length[1] = widths[1]
    run.value[1] = values[1]
    i = 2
  }
  j=2

  while (i < length(run.value)) {
    pre.width = bounds[j,1] - bounds[j-1,2]
    if (pre.width > 1) {
      run.length[i] = pre.width - 1
      i = i + 1
    }
    run.length[i] = widths[j]
    run.value[i] = values[j]
    j = j + 1
    i = i + 1
  }
  if (sum(run.length) != n) {
    stop("Rle is the wrong length. Look for double counted features in your bounds table.")
  }
  return( Rle( run.value, run.length ) )
}

##' Make Rle from segments for one sample
##'
##' Take output of CBS, make Rle representing all features in 'locs' ranges. CBS output contains
##' run length and run values for genomic segmetns, which could very directly be converted into a Rle.
##' However, as NA values are often removed, especially for mBAF data, these run lengths do not
##' necessarily cover all features in every sample. Using the start and top positions of each segment 
##' and the location of each feature, we can make a Rle that represents all features.
##' 
##' @param segs data.frame of segments, formatted as output of segment function from DNAcopy package
##' @param locs RangedData, like locData slot of a GenoSet
##' @return Rle with run lengths and run values covering all features in the data set.
##' @export
##' @family "segmented data"
##' @examples
##'   data(genoset)
##'   segs = runCBS( lrr(baf.ds), locData(baf.ds), return.segs=TRUE )
##'   segs2Rle( segs[[1]], locData(baf.ds) )  # Take a data.frame of segments, say from DNAcopy's segment function, and make Rle's using probe locations in the RangedData locs
##' @author Peter M. Haverty \email{phaverty@@gene.com}
segs2Rle <- function(segs, locs) {
  if ("num.mark" %in% colnames(segs)) {
    if (sum(segs[,"num.mark"]) == nrow(locs)) {
      return(Rle( segs[,"seg.mean"], segs[,"num.mark"]))
    }
  }
  seg.gr = GRanges( ranges=IRanges(start=segs[,"loc.start"], end=segs[,"loc.end"]),
    seqnames=factor(segs[,"chrom"],levels=chrOrder(unique(as.character(segs$chrom)))), "Value"=segs[,"seg.mean"])
  seg.gr = toGenomeOrder(seg.gr)
  temp.rle = Rle(seg.gr$Value[match(locs, seg.gr)])
  #  bounds = boundingIndicesByChr( seg.gr, locs )
  #  temp.rle = bounds2Rle( bounds, values(seg.gr)$Value, nrow(locs) )  # Breaks unit test for 1st being NA
  return(temp.rle)
}

##' Given segments, make a DataFrame of Rle objects for each sample
##'
##' Take table of segments from CBS, convert DataTable of Rle objects for each sample.
##' @title CBS segments to probe matrix
##' @param seg.list list, list of data frames, one per sample, each is result from CBS
##' @param locs locData from a GenoSet object
##' @return DataFrame of Rle objects with nrows same as locs and one column for each sample
##' @export segs2RleDataFrame
##' @family "segmented data"
##' @examples
##'   data(genoset)
##'   seg.list = runCBS( lrr(baf.ds), locData(baf.ds), return.segs=TRUE )
##'   segs2RleDataFrame( seg.list, locData(baf.ds) )  # Loop segs2Rle on list of data.frames in seg.list
##' @author Peter Haverty
segs2RleDataFrame <- function(seg.list, locs) {
  rle.list = lapply(seg.list, segs2Rle, locs)
  rle.data.frame = DataFrame(rle.list, row.names=featureNames(locs))
  return(rle.data.frame)
}

##' Make a RangedData from segments
##'
##' Starting from a data.frame of segments, like from CBS and segTable, organize as a RangedData. Label data "score",
##' so it can easily be made into various genome browser formats using rtracklayer.
##' @param segs data.frame, like from segment in DNAcopy or segTable
##' @return RangedData
##' @family "segmented data"
##' @export 
##' @author Peter M. Haverty \email{phaverty@@gene.com}
segs2RangedData <- function(segs) {
  rd = RangedData(ranges=IRanges(start=segs$loc.start,end=segs$loc.end),space=segs$chrom,score=segs$seg.mean,num.mark=segs$num.mark)
  return(rd)
}

##' Convert Rle objects to tables of segments
##'
##' Like the inverse of segs2Rle and segs2RleDataFrame. Takes a
##' Rle or a DataFrame with Rle
##' columns and the locData RangedData both from a GenoSet object
##' and makes a list of data.frames each like the result of CBS's
##' segment.  Note the loc.start and loc.stop will correspond
##' exactly to probe locations in locData and the input to
##' segs2RleDataFrame are not necessarily so. For a DataFrame, the
##' argument \code{stack} combines all of the individual data.frames
##' into one large data.frame and adds a "Sample" column of sample ids.
##'
##' For a Rle, the user can provide \code{locs} or \code{chr.ind},
##' \code{start} and \code{stop}.  The latter is surprisingly much faster
##' and this is used in the DataFrame version.
##'
##' @param object Rle or list/DataFrame of Rle vectors
##' @param locs RangedData with rows corresponding to rows of df
##' @param chr.ind matrix, like from chrIndices method
##' @param start integer, vector of feature start positions
##' @param end integer, vector of feature end positions
##' @param factor.chr scalar logical, make 'chrom' column a factor?
##' @return one or a list of data.frames with columns chrom, loc.start, loc.end, num.mark, seg.mean
##' @export segTable
##' @family "segmented data"
##' @examples
##'   data(genoset)
##'   seg.list = runCBS( lrr(baf.ds), locData(baf.ds), return.segs=TRUE )
##'   df = segs2RleDataFrame( seg.list, locData(baf.ds) )  # Loop segs2Rle on list of data.frames in seg.list
##'   assayDataElement( baf.ds, "lrr.segs" ) = df
##'   segTable( df, locData(baf.ds) )
##'   segTable( assayDataElement(baf.ds,"lrr.segs"), locData(baf.ds) )
##'   segTable( assayDataElement(baf.ds,"lrr.segs")[,1], locData(baf.ds), sampleNames(baf.ds)[1] )
##' @author Peter M. Haverty
##' @docType methods
##' @rdname segTable-methods
setGeneric("segTable", function(object,...) standardGeneric("segTable"))

##' @rdname segTable-methods
##' @aliases segTable,Rle-method
setMethod("segTable", signature(object="Rle"), function(object,locs=NULL,chr.ind=NULL,start=NULL,end=NULL,factor.chr=TRUE) {

  if (!is.null(locs)) {
    chr.ind = chrIndices(locs)
    start = start(locs)
    end = end(locs)
  } else {
    if (is.null(chr.ind) || is.null(start) || is.null(end)) {
      stop("If locs arg is not provided then chr.ind, start, and end must be provided.")
    }
  }

  # Get union of all breakpoints in Rle and chromosomes
  object.ends = cumsum(runLength(object))
  
  all.ends = sort(unique(c(chr.ind[,2],object.ends)))
  all.starts = c(1L,all.ends[-length(all.ends)]+1L)
  num.mark = (all.ends - all.starts) + 1L

  # Look up runValue with binary search on cumsum runValue starts. Starts rather than ends because findInterval is < rather than <=.
  object.starts = c(1L,object.ends[-length(object.ends)]+1L)
  object.vals = runValue(object)[ findInterval( all.ends, object.starts ) ]

  # Assign chrom,start,stop to each segment
  if (factor.chr == TRUE) {
    chrom = factor(rownames(chr.ind)[ findInterval(all.starts,chr.ind[,1]) ],levels=rownames(chr.ind))
  } else {
    chrom = rownames(chr.ind)[ findInterval(all.starts,chr.ind[,1]) ]
  }
  loc.end = end[all.ends]
  loc.start = start[all.starts]

  sample.seg = data.frame(chrom = chrom, loc.start = loc.start, loc.end = loc.end, num.mark = num.mark, seg.mean = object.vals, row.names=NULL, stringsAsFactors=FALSE, check.names=FALSE, check.rows=FALSE)
  return(sample.seg)
})

##' @rdname segTable-methods
##' @aliases segTable,DataFrame-method
##' @param stack logical, rbind list of segment tables for each sample and add "Sample" column?
setMethod("segTable", signature(object="DataFrame"), function(object,locs,factor.chr=TRUE,stack=FALSE) {
  internal.factor.chr = ifelse(factor.chr == TRUE && stack == FALSE,TRUE,FALSE)
  chr.ind = chrIndices(locs)
  start = start(locs)
  end = end(locs)
  
  segs = lapply( object,
    function(x) {
      return(segTable(x,chr.ind=chr.ind, start=start, end=end,factor.chr=internal.factor.chr))
    })
  if (stack == FALSE) {
    return(segs)
  } else {
    segs.df = do.call(rbind,segs)
    segs.df = cbind(data.frame(Sample = rep(names(segs),sapply(segs,nrow)),stringsAsFactors=FALSE),segs.df,row.names=NULL)
    if (factor.chr == TRUE) {
      chr.names = chrNames(locs)
      segs.df$chrom = factor(segs.df$chrom,levels=chr.names)
    }
    return(segs.df)
  }
})


##' Convert Rle objects to tables of segments
##'
##' Like segTable, but for two Rle objects. Takes a
##' pair of Rle or DataFrames with Rle
##' columns and makes one or more data.frames with bounds of each new 
##' segment.  Rle objects are broken up so that each resulting segment 
##' has one value from each Rle. For a DataFrame, the
##' argument \code{stack} combines all of the individual data.frames
##' into one large data.frame and adds a "Sample" column of sample ids.
##'
##' For a Rle, the user can provide \code{locs} or \code{chr.ind},
##' \code{start} and \code{stop}.  The latter is surprisingly much faster
##' and this is used in the DataFrame version.
##'
##' @param x Rle or list/DataFrame of Rle vectors
##' @param y Rle or list/DataFrame of Rle vectors
##' @param locs RangedData with rows corresponding to rows of df
##' @param chr.ind matrix, like from chrIndices method
##' @param start integer, vector of feature start positions
##' @param end integer, vector of feature end positions
##' @param factor.chr scalar logical, make 'chrom' column a factor?
##' @return one or a list of data.frames with columns chrom, loc.start, loc.end, num.mark, seg.mean
##' @export segPairTable
##' @family "segmented data"
##' @examples
##'   cn = Rle(c(3,4,5,6),rep(3,4))
##'   loh = Rle(c(2,4,6,8,10,12),rep(2,6))
##'   start = c(9:11,4:9,15:17)
##'   end = start
##'   locs = RangedData(IRanges(start=start,end=end),space=c(rep("chr1",3),rep("chr2",6),rep("chr3",3)))
##'   segPairTable(cn,loh,locs)
##' @author Peter M. Haverty
##' @docType methods
##' @rdname segPairTable-methods
setGeneric("segPairTable", function(x,y,...) standardGeneric("segPairTable"))

##' @rdname segPairTable-methods
##' @aliases segPairTable,Rle,Rle-method
setMethod("segPairTable", signature(x="Rle",y="Rle"), function(x,y,locs=NULL,chr.ind=NULL,start=NULL,end=NULL,factor.chr=TRUE) {
  # Fill in missing args if locs given
  # Maybe use ... rather than x and y and get names from that to use in colnames
  if (!is.null(locs)) {
    chr.ind = chrIndices(locs)
    start = start(locs)
    end = end(locs)
  } else {
    if (is.null(chr.ind) || is.null(start) || is.null(end)) {
      stop("If locs arg is not provided then chr.ind, start, and end must be provided.")
    }
  }

  # Get union of all breakpoints in two Rles and chromosomes
  x.ends = cumsum(runLength(x))
  y.ends = cumsum(runLength(y))
  
  all.ends = sort(unique(c(chr.ind[,2],x.ends,y.ends)))
  all.starts = c(1L,all.ends[-length(all.ends)]+1L)
  num.mark = (all.ends - all.starts) + 1L

  # Look up runValue with binary search on cumsum runValue starts. Starts rather than ends because findInterval is < rather than <=.
  x.starts = c(1L,x.ends[-length(x.ends)]+1L)
  x.vals = runValue(x)[ findInterval( all.ends, x.starts ) ]
  y.starts = c(1L,y.ends[-length(y.ends)]+1L)
  y.vals = runValue(y)[ findInterval( all.ends, y.starts ) ]

  # Assign chrom,start,stop to each segment
  if (factor.chr == TRUE) {
    chrom = factor(rownames(chr.ind)[ findInterval(all.starts,chr.ind[,1]) ],levels=rownames(chr.ind))
  } else {
    chrom = rownames(chr.ind)[ findInterval(all.starts,chr.ind[,1]) ]
  }
  loc.end = end[all.ends]
  loc.start = start[all.starts]

  # Make output object
  sample.seg = data.frame(chrom=chrom,loc.start = loc.start, loc.end = loc.end, num.mark = num.mark, x=x.vals, y=y.vals, row.names=NULL, stringsAsFactors=FALSE, check.names=FALSE, check.rows=FALSE)
  #  sample.seg = GRanges(ranges=IRanges(start=loc.start, end=loc.end), seqnames=chrom,  num.mark = num.mark, x=x.vals, y=y.vals)
  # Later, when incoming pinfo is a GRanges, will want to pass on chrlengths in new GRanges
  return(sample.seg)
})

##' @rdname segPairTable-methods
##' @aliases segPairTable,DataFrame,DataFrame-method
##' @param stack logical, rbind list of segment tables for each sample and add "Sample" column?
setMethod("segPairTable", signature(x="DataFrame",y="DataFrame"), function(x,y,locs,stack=FALSE,factor.chr=TRUE) {
  internal.factor.chr = ifelse(factor.chr == TRUE && stack == FALSE,TRUE,FALSE)
  chr.ind = chrIndices(locs)
  start = start(locs)
  end = end(locs)
  segs = mapply(
    function(one,two) {
      return(segPairTable(one,two,chr.ind=chr.ind, start=start, end=end,factor.chr=internal.factor.chr))
    },
    x,y,
    SIMPLIFY=FALSE
    )
  if (stack == FALSE) {
    return(segs)
  } else {
    segs.df = do.call(rbind,segs)
    segs.df = cbind(Sample = rep(names(segs),sapply(segs,nrow)),segs.df,stringsAsFactors=FALSE,row.names=NULL,check.rows=FALSE)
    if (factor.chr == TRUE) {
      chr.names = chrNames(locs)
      segs.df$chrom = factor(segs.df$chrom,levels=chr.names)
    }
    return(segs.df)
  }
})

##' Fix NA runs in a Rle
##'
##' Fix NA runs in a Rle when the adjacent runs have equal values
##' @param x Rle to be fixed
##' @param max.na.run integer, longest run of NAs that will be fixed
##' @return Rle
##' @export 
##' @author Peter M. Haverty
fixSegNAs <- function(x,max.na.run=3) {
  if (is.na(runValue(x)[1]) & runLength(x)[1] <= max.na.run) {
    runValue(x)[1] = runValue(x)[2]
  }
  if (is.na(runValue(x)[nrun(x)]) & runLength(x)[nrun(x)] <= max.na.run) {
    runValue(x)[nrun(x)] = runValue(x)[nrun(x)-1]
  }
  bad = which(is.na(runValue(x)) & runLength(x) <= max.na.run)
  bad = bad[ runValue(x)[bad-1] == runValue(x)[bad+1] ]
  runValue(x)[bad] = runValue(x)[bad+1]
  return(x)
}

##' Utility function to run CBS's three functions on one or more samples
##' 
##' Takes care of running CBS segmentation on one or more samples. Makes appropriate
##' input, smooths outliers, and segment
##' 
##' @title Run CBS Segmentation
##' @aliases runCBS
##' @param data numeric matrix with continuous data in one or more columns
##' @param locs RangeData, like locData slot of GenoSet
##' @param return.segs logical, if true list of segment data.frames return, otherwise a DataFrame of Rle vectors. One Rle per sample.
##' @param n.cores numeric, number of cores to ask mclapply to use
##' @param smooth.region number of positions to left and right of individual positions to consider when smoothing single point outliers
##' @param outlier.SD.scale number of SD single points must exceed smooth.region to be considered an outlier
##' @param smooth.SD.scale floor used to reset single point outliers
##' @param trim fraction of sample to smooth
##' @param alpha pvalue cutoff for calling a breakpoint
##' @return data frame of segments from CBS
##' @family "segmented data"
##' @export runCBS
##' @examples
##'     sample.names = paste("a",1:2,sep="")
##'     probe.names =  paste("p",1:30,sep="")
##'     ds = matrix(c(c(rep(5,20),rep(3,10)),c(rep(2,10),rep(7,10),rep(9,10))),ncol=2,dimnames=list(probe.names,sample.names))
##'     locs = RangedData(ranges=IRanges(start=c(1:20,1:10),width=1,names=probe.names),space=paste("chr",c(rep(1,20),rep(2,10)),sep=""))
##'   
##'     seg.rle.result = DataFrame( a1 = Rle(c(rep(5,20),rep(3,10))), a2 = Rle(c(rep(2,10),rep(7,10),rep(9,10))), row.names=probe.names )
##'     seg.list.result = list(
##'       a1 = data.frame( ID=rep("a1",2), chrom=factor(c("chr1","chr2")), loc.start=c(1,1), loc.end=c(20,10), num.mark=c(20,10), seg.mean=c(5,3), stringsAsFactors=FALSE),
##'       a2 = data.frame( ID=rep("a2",3), chrom=factor(c("chr1","chr1","chr2")), loc.start=c(1,11,1), loc.end=c(10,20,10), num.mark=c(10,10,10), seg.mean=c(2,7,9), stringsAsFactors=FALSE)
##'       )
##'
##'     runCBS(ds,locs)  # Should give seg.rle.result
##'     runCBS(ds,locs,return.segs=TRUE) # Should give seg.list.result
##' @author Peter M. Haverty
runCBS <- function(data, locs, return.segs=FALSE, n.cores=1, smooth.region=2, outlier.SD.scale=4, smooth.SD.scale=2, trim=0.025, alpha=0.001) {
  if (!requireNamespace("DNAcopy",quietly=TRUE)) {
    stop("Failed to require DNAcopy package.\n")
  }
  sample.name.list = colnames(data)
  names(sample.name.list) = sample.name.list
  loc.start = as.numeric(start(locs))
  loc.chr = chr(locs)
  presorted = isGenomeOrder(locs,strict=TRUE)
  
  # mclapply over samples. cbs can loop over the columns of data, but want to use multiple forks
  if (n.cores > 1 && is.loaded("mc_fork", PACKAGE="parallel")) {
    mcLapply <- get('mclapply', envir=getNamespace('parallel'))
    loopFunc = function(...) { mcLapply(...,mc.cores=n.cores, mc.preschedule=FALSE) }
    cat("Using mclapply for segmentation ...\n")
  } else {
    loopFunc = lapply
  }
  segs = loopFunc(sample.name.list,
    function(sample.name) {
      writeLines(paste("Working on segmentation for sample number",match(sample.name,sample.name.list),":",sample.name))
      temp.data = as.numeric(data[,sample.name,drop=TRUE])
      ok.indices = !is.na(temp.data)
      CNA.object <- DNAcopy::CNA(temp.data[ok.indices], loc.chr[ok.indices], loc.start[ok.indices], data.type = "logratio", sampleid = sample.name, presorted=presorted)
      smoothed.CNA.object <- DNAcopy::smooth.CNA(CNA.object, smooth.region=smooth.region, outlier.SD.scale=outlier.SD.scale, smooth.SD.scale=smooth.SD.scale, trim=trim)
      segment.smoothed.CNA.object <- DNAcopy::segment(smoothed.CNA.object, verbose=0, alpha=alpha)$output
      if (return.segs == TRUE) {
        segment.smoothed.CNA.object$chrom = factor(as.character(segment.smoothed.CNA.object$chrom),levels=chrNames(locs))
        return(segment.smoothed.CNA.object)
      } else {
        return(segs2Rle(segment.smoothed.CNA.object,locs))
      }
    })

  if (return.segs == TRUE) {
    return(segs)
  } else {
    return( DataFrame(segs, row.names=featureNames(locs) ) )
  }
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

##' Order chromosome names in proper genome order
##'
##' Chromosomes make the most sense orded by number, then by letter.
##' 
##' @param chr.names character, vector of unique chromosome names
##' @return character vector of chromosome names in proper order
##' @export chrOrder
##' @family "genome ordering"
##' @examples
##'    chrOrder(c("chr5","chrX","chr3","chr7","chrY"))  #  c("chr3","chr5","chr7","chrX","chrY")
##' @author Peter M. Haverty
chrOrder <- function(chr.names) {
  simple.names = gsub("^chr","",chr.names)
  name.is.numeric = grepl("^[0-9]+$",simple.names,perl=T)
  numeric.names = chr.names[name.is.numeric][ order(as.numeric(simple.names[name.is.numeric])) ]
  non.numeric.names = chr.names[! name.is.numeric][ order(chr.names[ !name.is.numeric]) ]
  all.names = c(numeric.names,non.numeric.names)
  return(all.names)
}

##' Check if a GRanges, GenoSet or RangedData is in genome order
##'
##' Checks that rows in each chr are ordered by start.  If strict=TRUE, then chromosomes
##' must be in order specified by chrOrder. isGenomeOrder for GRanges differs from order
##' in that it orders by chromsome and start position only,
##' rather than chromsome, strand, start, and width.
##' 
##' @param ds GenoSet, GRanges, or RangedData
##' @param strict logical, should space/chromosome order be identical to that from chrOrder?
##' @return logical
##' @export isGenomeOrder
##' @family "genome ordering"
##' @examples
##'   data(genoset)
##'   isGenomeOrder( locData(genoset.ds) )
##' @author Peter M. Haverty
##' @rdname isGenomeOrder-methods
setGeneric("isGenomeOrder", function(ds,...) standardGeneric("isGenomeOrder"))

##' @aliases isGenomeOrder,RangedDataOrGenoSet-method
##' @rdname isGenomeOrder-methods
setMethod("isGenomeOrder",signature=signature(ds="RangedDataOrGenoSet"),
          function(ds, strict=TRUE) {
            if (strict) {
              if ( ! all( chrNames(ds) == chrOrder( chrNames(ds) ) ) ) {
                return(FALSE)
              }
            }
            # Check each chr for ordered start
            chr.ind = chrIndices(ds)
            return(!any(aggregate( start(ds), start=chr.ind[,1], end=chr.ind[,2], FUN=is.unsorted)))
          })

##' @aliases isGenomeOrder,GRanges-method
##' @rdname isGenomeOrder-methods
setMethod("isGenomeOrder",signature=signature(ds="GRanges"),
          function(ds, strict=TRUE) {
            if ( any(duplicated(runValue(seqnames(ds)))) ) { stop("GRanges not in blocks by chromosome.") }
            if (strict == TRUE) {
              if (!isTRUE(all.equal(chrOrder(seqlevels(ds)), seqlevels(ds)))) {
                return(FALSE)
              }
            }
            chr.ind = chrIndices(ds)
            return(!any(aggregate( start(ds), start=chr.ind[,1], end=chr.ind[,2], FUN=is.unsorted)))
          })

##' Set a GRanges, GenoSet, or RangedData to genome order
##'
##' Returns a re-ordered object sorted by chromosome and start position. If strict=TRUE, then
##' chromosomes must be in order specified by chrOrder.
##' If ds is already ordered, no re-ordering is done. Therefore, checking order with isGenomeOrder,
##' is unnecessary if order will be corrected if isGenomeOrder is FALSE.
##'
##' toGenomeOrder for GRanges differs from sort in that it orders by chromsome and start position only,
##' rather than chromsome, strand, start, and width.
##' 
##' @param ds GenoSet, GRanges, or RangedData
##' @param strict logical, should chromosomes be in order specified by chrOrder?
##' @return re-ordered ds
##' @export toGenomeOrder
##' @examples
##'   data(genoset)
##'   toGenomeOrder( baf.ds, strict=TRUE )
##'   toGenomeOrder( baf.ds, strict=FALSE )
##'   toGenomeOrder( locData(baf.ds) )
##' @author Peter M. Haverty
##' @docType methods
##' @family "genome ordering"
##' @rdname toGenomeOrder-methods
setGeneric("toGenomeOrder", function(ds,...) standardGeneric("toGenomeOrder"))

##' @rdname toGenomeOrder-methods
##' @aliases toGenomeOrder,RangedData-method
setMethod("toGenomeOrder",signature=signature(ds="RangedData"),
          function(ds, strict=TRUE) {
            if (strict == TRUE) {
              if (!isTRUE(all.equal(chrOrder(chrNames(ds)), chrNames(ds)))) {
                ds = ds[ chrOrder(chrNames(ds)) ]
              }
            }
            row.order = order(as.integer(space(ds)),start(ds))
            if (is.unsorted(row.order)) {
              return( ds[row.order,,drop=FALSE] )
            } else {
              return( ds )
            }
          })

##' @rdname toGenomeOrder-methods
##' @aliases toGenomeOrder,GRanges-method
setMethod("toGenomeOrder",signature=signature(ds="GRanges"),
          function(ds, strict=TRUE) {
            if (strict == TRUE) {
              if (!isTRUE(all.equal(chrOrder(seqlevels(ds)), seqlevels(ds)))) {
                seqlevels(ds) = chrOrder(seqlevels(ds))
              }
            }
            row.order = order(as.integer(seqnames(ds)),start(ds))
            if (is.unsorted(row.order)) {
              ds = ds[row.order,,drop=FALSE]
            }
            return(ds)
          })

##' @rdname toGenomeOrder-methods
##' @aliases toGenomeOrder,GenoSet-method
setMethod("toGenomeOrder", signature=signature(ds="GenoSet"),
          function(ds,strict=TRUE) {
            locData(ds) = toGenomeOrder(locData(ds),strict=strict) # locData<- fixes row ordering in ds
            return(ds)
          })

##' Load a GenoSet from a RData file
##'
##' Given a RData file with one object (a GenoSet or related object), load it,
##' and return.
##' @param path character, path to RData file
##' @return GenoSet or related object (only object in RData file)
##' @examples
##' \dontrun{ ds = readGenoSet("/path/to/genoset.RData") }
##' @export 
##' @author Peter M. Haverty \email{phaverty@@gene.com}
readGenoSet <- function(path) {
  object = get(load(path)[1])
  if (!is(object,"eSet")) { stop("Loaded object is not an eSet or derived class.") }
  return( object )
}
