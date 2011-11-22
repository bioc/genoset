######  Class definition for GenoSet, which will extend eSet
######   GenoSet will provide a locData slot containing a RangedData object from the IRanges
######   package to hold genome locations of the features and allow for easy subsetting
######   by location.
######   Intended to be subset by other classes to add one or more data matrices to
######   the assayData slot.

##' @importClassesFrom Biobase AnnotatedDataFrame AssayData eSet ExpressionSet MIAME Versioned VersionedBiobase
##' @importClassesFrom IRanges DataFrame RangedData RangesList Rle
##' @importClassesFrom methods ANY character matrix numeric
##' @importClassesFrom BSgenome BSgenome
##' @importFrom Biostrings getSeq letterFrequency
##' @importMethodsFrom Biobase annotation experimentData exprs fData featureNames "featureNames<-" phenoData sampleNames "sampleNames<-"
##' @importMethodsFrom IRanges as.data.frame as.list as.matrix cbind colnames "colnames<-" elementLengths end findOverlaps gsub
##' @importMethodsFrom IRanges "%in%" intersect is.unsorted lapply levels match mean na.exclude nrow order paste ranges Rle rownames
##' @importMethodsFrom IRanges "rownames<-" runLength runValue sapply space start universe "universe<-" unlist
##' @importMethodsFrom methods coerce show
##'
##' @importFrom Biobase assayDataElement assayDataElementNames assayDataElementReplace assayDataNew annotatedDataFrameFrom
##'
##' @importFrom DNAcopy CNA segment smooth.CNA
##'
##' @importFrom graphics abline axis axTicks box mtext plot plot.new plot.window points segments
##'
##' @importFrom IRanges DataFrame IRanges RangedData
##'
##' @importFrom methods slot "slot<-" callNextMethod is new
##'
##' @importFrom stats density lm residuals
##'
##' @importFrom GenomicRanges seqlengths
##'
##' @importFrom bigmemory as.big.matrix attach.big.matrix is.nil is.big.matrix
##' 
##' @include DataFrame-methods.R
##' @useDynLib genoset


###############
# Class GenoSet
###############

##' @exportClass GenoSet
setClass("GenoSet", contains=c("eSet"), representation=representation(locData="RangedData"))

setValidity("GenoSet", function(object) {
  return( all( rownames(locData(object)) == featureNames(object) ) )
})

# Create class union of GenoSet and RangedData so method signatures can be set for either
setClassUnion("RangedDataOrGenoSet",c("RangedData","GenoSet"))

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
##' @param locData A RangedData object specifying feature chromosome
##' locations. Rownames are required to match featureNames.
##' @param pData A data frame with rownames matching all data matrices
##' @param annotation character, string to specify chip/platform type
##' @param universe character, a string to specify the genome universe
##' for locData
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
initGenoSet <- function(type, locData, pData=NULL, annotation="", universe=NULL, ...) {
  # Function to clean up items for slots and call new for GenoSet and its children
  # ... will be the matrices that end up in assayData
  # all dimnames "fixed" with make names because eSet is inconsistent about that

  if (is.null(universe)) {
    if( is.null(universe(locData)) ) {
      stop("Arg universe must be provided")
    } else {
      # Do nothing, universe is set
    }
  } else {
    universe(locData) = universe
  }
  
  # Check/set genome order of locData
  clean.loc.rownames = make.names(rownames(locData),unique=TRUE)
  if ( ! all(rownames(locData) == clean.loc.rownames) ) {
    rownames(locData) = clean.loc.rownames
  }
  if ( ! isGenomeOrder(locData, strict=TRUE) ) {
    locData = toGenomeOrder(locData, strict=TRUE )
  }

 # Create assayData
  ad = assayDataNew(storage.mode="environment",...)
  clean.featureNames = make.names(featureNames(ad),unique=TRUE)
  if ( ! all(featureNames(ad) == clean.featureNames) ) {
    featureNames(ad) = clean.featureNames
  }

  # Check colnames of all data matrices identical and set to same order if necessary
  first.name = assayDataElementNames(ad)[1]
  for (mat.name in assayDataElementNames(ad)[-1]) {
    if (! setequal(colnames(ad[[mat.name]]), colnames(ad[[first.name]]) ) ) {
      stop(paste("Mismatch between rownames of first data matrix and", mat.name))
    }
    if ( any( colnames(ad[[mat.name]]) != colnames(ad[[first.name]])) ) {
      ad[[mat.name]] == ad[[mat.name]][,colnames(ad[[first.name]])]
    }
  }

  # Check sampleNames are same as check.names
  clean.sampleNames = make.names(sampleNames(ad),unique=TRUE)
  if ( ! all(sampleNames(ad) == clean.sampleNames) ) {
    sampleNames(ad) = clean.sampleNames
  }
  
  # Set row order to match locData
  for (  ad.name in assayDataElementNames(ad) ) {
    if (! setequal(rownames(ad[[ad.name]]), rownames(locData)) ) {
      stop(paste("Mismatch between rownames of data matrix", ad.name, "and probe location info 'locData'"))
    }
    if ( any( rownames(ad[[ad.name]]) != rownames(locData) ) ) {
      ad[[ad.name]] = ad[[ad.name]][rownames(locData),]
    }
  }

  # Done editing assayData members, lock
  lockEnvironment(ad, bindings=TRUE)
  
  # Create or check phenoData
  if (is.null(pData)) {
    pData = data.frame(Sample=sampleNames(ad),row.names=sampleNames(ad))
  } else {
    rownames(pData) = make.names(rownames(pData),unique=TRUE)
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
GenoSet <- function(locData, pData=NULL, annotation="", universe=NULL, ...) {
  object = initGenoSet(type="GenoSet", locData=locData, pData=pData, annotation=annotation, universe=universe, ...)
  return(object)
}

#########
# Methods
#########

#####################
# Getters and Setters
#####################

setMethod("sampleNames<-", signature(object="GenoSet",value="ANY"),
          function(object,value) {
            value = make.names(value,unique=TRUE)
            object = callNextMethod(object,value)
            return(object)
          })

##' Set featureNames
##'
##' Set featureNames including rownames of position info
##' @title Set featureNames
##' @param object GenoSet 
##' @param value ANY
##' @return A new object of the class of supplied object
##' @exportMethod "featureNames<-"
##' @author Peter M. Haverty
setMethod("featureNames<-",
                 signature=signature(object="GenoSet", value="ANY"),
                 function(object, value) {
                   value = make.names(value,unique=TRUE)
                   object = callNextMethod(object,value)
                   rownames(slot(object,"locData")) = value
                   return(object)
                 })

##' Access the feature genome position info
##'
##' The position information for each probe/feature is stored as an IRanges RangedData object.
##' The locData functions allow this data to be accessed or re-set.
##'
##' @title Get and set probe set info
##' @param object GenoSet
##' @export locData
##' @author Peter M. Haverty
##' @param object A GenoSet object
##' @rdname locData
##' @examples
##'   data(genoset)
##'   rd = locData(genoset.ds)
##'   locData(genoset.ds) = rd
setGeneric("locData", function(object) standardGeneric("locData"))
##' @rdname locData
setMethod("locData", "GenoSet", function(object) { return(slot(object,"locData")) } )

##' @rdname locData
setGeneric("locData<-", function(object,value) standardGeneric("locData<-") )

##' Set locData
##'
##' Set locData
##' @title Set position info
##' @param object GenoSet
##' @param value RangedData describing features
##' @return A GenoSet object
##' @author Peter Haverty
##' @export "locData<-"
##' @rdname locData
setMethod("locData<-", signature(object="GenoSet", value="RangedData"),
                 function(object,value) {
                   if (! all( rownames(value) %in% featureNames(object))) {
                       stop("Can not replace locData using rownames not in this GenoSet")
                     }
                   slot(object,"locData") = value
                   if (! all(rownames(value) == featureNames(object))) {
                     for (adname in assayDataElementNames(object)) {
                       assayDataElement(object,adname) = assayDataElement(object,adname)[rownames(value),]
                     }
                   }
                   return(object)
                   })

##' Genome universe for locData
##'
##' The genome positions of the features in locData. The UCSC notation (e.g. hg18, hg19, etc.) should be used.
##'
##' @title Get and set the genome universe annotation.
##' @param x GenoSet
##' @return character, e.g. hg19
##' @author Peter M. Haverty
##' @exportMethod universe
##' @rdname universe
##' @examples
##'   data(genoset)
##'   universe(genoset.ds)
##'   universe(genoset.ds) = "hg19"
setMethod("universe", "GenoSet", function(x) { return(universe(x@locData)) } )

##' Set genome universe
##'
##' Set genome universe
##' 
##' @param x GenoSet
##' @param value character, new universe string, e.g. hg19
##' @return A GenoSet object
##' @author Peter Haverty
##' @exportMethod "universe<-"
##' @rdname universe
setMethod("universe<-", signature(x="GenoSet"),
                 function(x,value) {
                   universe(x@locData) = value
                   return(x)
                   })

###########################################
# Shared API between GenoSet and RangedData
###########################################

##' Get space factor for GenoSet
##'
##' locData slot holds a RangedData, which keeps the chromsome of each
##' feature in a factor names 'space'.
##' @param x GenoSet
##' @return factor
##' @author Peter M. Haverty
##' @rdname genoset-methods
##' @examples
##' data(genoset)
##' space(genoset.ds)
##' start(genoset.ds)
##' end(genoset.ds)
##' names(genoset.ds)
##' ranges(genoset.ds) # Returns a RangesList
##' elementLengths(genoset.ds) # Returns the number of probes per chromosome
setMethod("space", "GenoSet", function(x) { return(space(locData(x))) } )

##' Get start of location for each feature
##'
##' locData slot holds a RangedData.
##' @param x GenoSet
##' @return integer
##' @author Peter M. Haverty
##' @rdname genoset-methods
setMethod("start", "GenoSet", function(x) { return(start(locData(x))) } )

##' Get space factor for GenoSet
##'
##' locData slot holds a RangedData.
##' @param x GenoSet
##' @return integer
##' @author Peter M. Haverty
##' @rdname genoset-methods
setMethod("end", "GenoSet", function(x) { return(end(locData(x))) } )

##' Get chromosome names
##'
##' Get chromosome names, which are the names of the locData slot.
##' @title Names for chromosome
##' @param x GenoSet
##' @return character
##' @author Peter Haverty
##' @exportMethod names
##' @rdname genoset-methods
setMethod("names", "GenoSet", function(x) { return( names(locData(x)) ) } )

##' Get ranges from locData slot
##'
##' Get ranges from locData slot
##' @title Ranges for chromosome
##' @param x GenoSet
##' @return character
##' @author Peter Haverty
##' @exportMethod ranges
##' @rdname genoset-methods
setMethod("ranges", "GenoSet", function(x) { return( ranges(locData(x)) ) } )

##' Get elementLengths from locData slot
##'
##' Get elementLengths from locData slot
##' @title ElementLengths for chromosome
##' @param x GenoSet
##' @return character
##' @author Peter Haverty
##' @exportMethod elementLengths
##' @rdname genoset-methods
setMethod("elementLengths", "GenoSet", function(x) { return( elementLengths(locData(x)) ) } )

#############
# Sub-setters
#############

##' @exportMethod "["
##' @param x GenoSet
##' @param i character, RangedData, RangesList, logical, integer
##' @param j character, RangedData, RangesList, logical, integer
##' @param k chracter or integer
##' @param drop logical drop levels of space factor?
##' @param ... additional subsetting args
##' @examples
##'   data(genoset)
##'   genoset.ds[1:5,2:3]  # first five probes and samples 2 and 3
##'   genoset.ds[ , "K"]  # Sample called K
##'   rd = RangedData(ranges=IRanges(start=seq(from=15e6,by=1e6,length=7),width=1),names=letters[8:14],space=rep("chr17",7))
##'   genoset.ds[ rd, "K" ]  # sample K and probes overlapping those in rd, which overlap specifed ranges on chr17
##'  @rdname genoset-methods
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
              i = match(rownames(x@locData),featureNames(x)) # Re-ordering of RangedData can silently disobey in order to keep its desired order of chromosomes
            }
            callNextMethod(x,i,j,...,drop=drop)
          })
##' @rdname genoset-methods
setMethod("[", signature=signature(x="GenoSet",i="character",j="ANY"),
          function(x,i,j,...,drop=FALSE) {
            if ( ! missing(i) ) {
              indices = match(i,featureNames(x))
            }
            callNextMethod(x,indices,j,...,drop=drop)
          })

##' @rdname genoset-methods
setMethod("[", signature=signature(x="GenoSet", i="RangedData", j="ANY"),
          function(x,i,j,...,drop=FALSE) {
            indices = unlist(x@locData %in% i)
            callNextMethod(x,indices,j,...,drop=drop)
          })

##' @rdname genoset-methods
setMethod("[", signature=signature(x="GenoSet", i="RangesList", j="ANY"),
          function(x,i,j,...,drop=FALSE) {
            indices = unlist(x@locData %in% i)
            callNextMethod(x,indices,j,...,drop=drop)
          })

#######
# Other
#######

##' @exportMethod show
setMethod("show","GenoSet",
          function(object) {
            callNextMethod(object)
            cat("Feature Locations:\n")
            show(slot(object,"locData"))
            cat("Universe: ",universe(object),"\n")
          })

########################
# Get genome information
########################

##' Chromsome name for each feature
##'
##' Get chromosome name for each feature.  Returns character, not the factor 'space'.
##' @title Look up chromosome for each feature 
##' @param object RangedData or GenoSet
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
setMethod("chr", "RangedData", function(object) { return(as.character(space(object))) } )
##' @rdname chr-methods
setMethod("chr", "GenoSet", function(object) { return(as.character(space(slot(object,"locData")))) } )

##' Chromosome position of features
##'
##' Get chromsome position of features/ranges. Defined as floor of mean of start and end.
##' @title Positions for features
##' @param object RangedData or GenoSet
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
setGeneric("pos", function(object) standardGeneric("pos"))
##' @rdname pos
setMethod("pos", "RangedData", function(object) { return( (start(object) + end(object)) %/% 2L ) } )
##' @rdname pos
setMethod("pos", "GenoSet",    function(object) { return( (start(locData(object)) + end(locData(object))) %/% 2L ) } )

##' Get list of unique chromosome names
##'
##' Get list of unique chromosome names. A synonym for names().
##' 
##' @param object RangedData or GenoSet
##' @return character vector with names of chromosomes
##' @author Peter M. Haverty
##' @export uniqueChrs
##' @examples
##'   test.sample.names = LETTERS[11:13]
##'   probe.names = letters[1:10]
##'   gs = GenoSet(
##'      locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18"),
##'      cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
##'      pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
##'      annotation="SNP6"
##'   )
##'   uniqueChrs(gs) # c("chr1","chr3","chrX")
##'   uniqueChrs(locData(gs))  # The same
##' @rdname uniqueChrs
setGeneric("uniqueChrs", function(object) standardGeneric("uniqueChrs") )
##' @rdname uniqueChrs
setMethod("uniqueChrs", signature(object="RangedDataOrGenoSet"),
          function(object) {
            names(object)
          })


##' Get chromosome names in genome order
##' 
##' Get chromosome names from locData data in a GenoSet.  Order numerically, for
##' numeric chromosomes, then lexically for the rest.
##' 
##' @param object GenoSet or RangedData
##' @return character vector with chrs in genome order
##' @examples
##'   test.sample.names = LETTERS[11:13]
##'   probe.names = letters[1:10]
##'   gs = GenoSet(
##'      locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chrX",2),rep("chr3",4)),universe="hg18"),
##'      cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
##'      pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
##'      annotation="SNP6"
##'   )
##'   orderedChrs(gs) # c("chr1","chr3","chrX")
##'   orderedChrs(locData(gs))  # The same       
##' @author Peter M. Haverty
##' @export orderedChrs
##' @rdname orderedChrs
setGeneric("orderedChrs", function(object) standardGeneric("orderedChrs") )
##' @rdname orderedChrs
setMethod("orderedChrs", signature(object="RangedDataOrGenoSet"),
          function(object) {
            chr.names = names(object)
            chr.names = chrOrder(chr.names)
            return(chr.names)
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
setMethod("chrInfo", signature(object="RangedDataOrGenoSet"),
          function(object) {
            # Get max end value for each chr
            max.val = as.list(max(end(ranges(object))))
            
            # Alternatively, get from R library storing chr info for this genome
            #library(org.Hs.eg.db)
            #max.val = as.list(org.Hs.egCHRLENGTHS)
            
            max.val = max.val[ orderedChrs(object) ]
            
            chr.info = matrix(ncol=3,nrow=length(max.val), dimnames=list(names(max.val),c("start","stop","offset")))
            chr.info[,"stop"]    = cumsum(max.val)
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
##' @param object GenoSet or RangedData
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
setMethod("chrIndices", signature(object="RangedDataOrGenoSet"),
          function(object,chr=NULL) {
            chr.names = names(object)
            chr.info = matrix(ncol=3,nrow=length(chr.names), dimnames=list(chr.names,c("first","last","offset")))
            chr.info[,"last"] = cumsum( elementLengths(object) )
            chr.info[,"first"] = c(1,chr.info[- nrow(chr.info),"last"] + 1)
            if (!is.null(chr)) {
              if (! chr %in% names(object)) { stop("Must specify a valid chromosome name in chrIndices.\n") }
              return( seq.int( chr.info[chr,"first"], chr.info[chr,"last"]) )
            }
            chr.info[,"offset"] = chr.info[,"first"] -1
            return(chr.info)
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
##'   genoPos(genoset.ds)
##'   genoPos(locData(genoset.ds))  # The same
##' @export genoPos
##' @rdname genoPos-methods
setGeneric("genoPos", function(object) standardGeneric("genoPos") )
##' @rdname genoPos-methods
setMethod("genoPos", signature(object="RangedDataOrGenoSet"),
          function(object) {

            # For single chr objects, just return pos
            if ( length(names(object)) == 1 ) {
              return(pos(object))
            }
            
            ### Add offset to pos by chr
            offset = chrInfo(object)[,"offset"]
            genopos = pos(object) + unlist(offset[chr(object)])
            
            return(genopos)
          })

#######
# Plots
#######

##' Plot data along the genome
##'
##' For a GenoSet object, data for a specified sample in a specified assayDataElement
##' can be plotted along the genome.  One chromosome can be specified if desired. If
##' more than one chromosome is present, the chromosome boundaries will be marked.
##' Alternatively, for a numeric x and a
##' numeric or Rle y, data in y can be plotted at genome positions y. In this case,
##' chromosome boundaries can be taken from the argument locs. If data for y-axis comes
##' from a Rle, either specified directly or coming from the specified assayData
##' element and sample, lines are plotted representing segments.
##' 
##' @param sample A index or sampleName to plot
##' @param element character, name of element in assayData to plot
##' @param x GenoSet (or descendant) or numeric with chromosome or genome positions
##' @param y numeric or Rle, values to be used for y-dimension, run start and stop indices or numeric with all values mapped to values in x for x-dimension or index of sample to be plotted if x is a GenoSet.
##' @param element character, when x is a GenoSet, the name of the assayDataElement to plot from.
##' @param locs RangedData, like locData slot of GenoSet
##' @param chr Chromosome to plot, NULL by default for full genome
##' @param add Add plot to existing plot
##' @param xlab character, label for x-axis of plot
##' @param ylab character, label for y-axis of plot
##' @param col character, color to plot lines or points
##' @param lwd numeric, line width for segment plots from an Rle
##' @param pch character or numeric, printing charactater, see points
##' @param ... Additional plotting args
##' @return nothing
##' @author Peter M. Haverty
##' @export genoPlot
##' @examples
##'   data(genoset)
##'   genoPlot( baf.ds,1,element="lrr")
##'   genoPlot( genoPos(baf.ds), assayDataElement(baf.ds,"lrr")[,1], locs=locData(baf.ds) ) # The same
##'   genoPlot( 1:10, Rle(c(rep(0,5),rep(3,4),rep(1,1))) )
##' @rdname genoPlot
setGeneric("genoPlot", function(x,y,...) { standardGeneric("genoPlot") } )
##' @rdname genoPlot
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

##' @rdname genoPlot
setMethod("genoPlot", c(x="numeric",y="Rle"),
          function(x, y, add=FALSE, xlab="", ylab="", col="red", locs=NULL, lwd=2, ...) {
            if (add == FALSE) {
              plot.new()
              plot.window(range(x,na.rm=TRUE),range(y,na.rm=TRUE),xlab=xlab,ylab=ylab,xaxs="i",...)
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

##' @rdname genoPlot
setMethod("genoPlot", signature(x="GenoSet",y="ANY"), function(x, y, element, chr=NULL, add=FALSE, pch=".", xlab="", ylab="", ...) {

  # Get position info, subset by chr if necessary
  if (! element %in% assayDataElementNames(x)) {
    stop("Provided assayData element, ", element, " is not a valid name of an assayData member")
  }
  if ( !is.null(chr) ) {
    indices = chrIndices(x,chr)
    element.values = assayDataElement(x,element)[indices,y]
    positions = start(x)[indices]
    locs = NULL
  } else {
    element.values = assayDataElement(x,element)[,y]
    positions = genoPos(x)
    if (length(uniqueChrs(x)) > 1) {
      locs = locData(x)
    } else {
      locs = NULL
    }
  }
  genoPlot(positions,element.values,locs=locs,add=add,xlab=xlab,ylab=ylab,pch=pch,...)
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
##' @examples
##'   data(genoset)
##'   genoPlot(genoPos(baf.ds), baf(baf.ds)[,1])
##'   genomeAxis( locs=locData(baf.ds) )  # Add chromsome names and boundaries to a plot assuming genome along x-axis
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

##' Load local GC percentage around features
##'
##' Local GC content  can be used to remove GC artifacts from copynumber data
##' see Diskin, 2008). GC% column will be added to the feature data.  The dataset
##' may be truncated to remove positions without GC information.  GC data are
##' accessible with locData(). Uses a cool BSgenome trick from Michael Lawrence.
##' This takes 5.6 hours for 2Mb windows on 2.5M probes, so look for some custom C
##' in future releases.
##' 
##' @param object A GenoSet object or derivative
##' @param expand numeric, expand each feature location by this many bases on each side
##' @param bsgenome, sequence db object from BSgenome (e.g. Hsapiens)
##' @return An updated object, with GC percentage information added to the locData slot.
##' @export loadGC
##' @rdname genoset-methods
##' @author Peter M. Haverty
setGeneric("loadGC", function(object,expand,bsgenome) standardGeneric("loadGC"))
##' @rdname genoset-methods
setMethod("loadGC", signature=signature(object="RangedData",expand="numeric",bsgenome="BSgenome"),
          function(object,expand=1e6,bsgenome) {
            expanded.ranges = ranges(object) + expand # Zoom
            # Check chr name matches
            if ( ! all( grepl("^chr",names(expanded.ranges)))) {
              names(expanded.ranges) = gsub("^chr","",names(expanded.ranges)) # Get rid of any chr prefixes that may exist
              names(expanded.ranges) = paste("chr",names(expanded.ranges),sep="") # Add chr prefix to all
            }
            # Check and fix zooming off of the chr ends
            expanded.ranges = restrict(expanded.ranges, start=1L,
              end=seqlengths(bsgenome)[ as.character(space(expanded.ranges)) ],
              keep.all.ranges=TRUE)
            # Get seqs and get GC content
            allSeqs = getSeq(bsgenome, expanded.ranges, as.character=FALSE)
            object$gc = letterFrequency(allSeqs,letters=c("GC"),as.prob=TRUE)[,1]
            return(object)
          })

##' @rdname genoset-methods
setMethod("loadGC", signature=signature(object="GenoSet",expand="numeric",bsgenome="BSgenome"), function(object,expand=1e6,bsgenome) {
  # Load gc into locData of GenoSet object
  ds = loadGC( locData(object),expand,bsgenome )
  locData(object) = ds
  return(object)
})


##' Correct copy number for GC content
##'
##' Copy number estimates from various platforms show "Genomic Waves" (Diskin et al.,
##' Nucleic Acids Research, 2008) where copy number trends with local GC content.
##' This function regresses copy number on GC percentage and removes the effect
##' (returns residuals). GC content should be smoothed along the genome in wide
##' windows >= 100kb.
##' 
##' @title cgCorrect
##' @param ds numeric matrix of copynumber or log2ratio values, samples in columns
##' @param gc numeric vector, GC percentage for each row of ds, must not have NAs
##' @param retain.mean logical, center on zero or keep same mean?
##' @return numeric matrix, residuals of ds regressed on gc
##' @export gcCorrect
##' @examples
##'   gc = runif(n=100, min=1, max=100)
##'   ds = rnorm(100) + (0.1 * gc)
##'   gcCorrect(ds, gc)
##' @author Peter M. Haverty
gcCorrect <- function(ds, gc, retain.mean=TRUE) {
  fit = lm( ds ~ gc, na.action=na.exclude )
  ds.fixed = residuals(fit)
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
  column.modes = apply(ds,2, function(x) { 
    l2r.density = density(x,na.rm=TRUE)
    density.max.index = which.max(l2r.density$y)
    return(l2r.density$x[density.max.index])
  })
  ds = sweep(ds, 2, column.modes)
  return(ds)
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
##' @examples
##'   data(genoset)
##'   segs = runCBS( lrr(baf.ds), locData(baf.ds), return.segs=TRUE )
##'   segs2Rle( segs[[1]], locData(baf.ds) )  # Take a data.frame of segments, say from DNAcopy's segment function, and make Rle's using probe locations in the RangedData locs
##' @author Peter M. Haverty \email{phaverty@@gene.com}
segs2Rle <- function(segs, locs) {
#  if (sum(segs[,"num.mark"],na.rm=TRUE) == nrow(locs)) {
#    return(Rle( segs[,"seg.mean"], segs[,"num.mark"]))
#  } else {
    temp.rle = Rle(as.numeric(NA),nrow(locs))
    seg.rd = RangedData( ranges=IRanges(start=segs[,"loc.start"], end=segs[,"loc.end"]),
      space=segs[,"chrom"], "Value"=segs[,"seg.mean"])
    seg.overlap = as.matrix( findOverlaps(seg.rd, locs ) )
    temp.rle[ seg.overlap[,2], drop=FALSE ] = seg.rd[ seg.overlap[,1], ]$Value
#  }
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
##' @examples
##'   data(genoset)
##'   seg.list = runCBS( lrr(baf.ds), locData(baf.ds), return.segs=TRUE )
##'   segs2RleDataFrame( seg.list, locData(baf.ds) )  # Loop segs2Rle on list of data.frames in seg.list
##' @author Peter Haverty
segs2RleDataFrame <- function(seg.list, locs) {
  rle.list = lapply(seg.list, segs2Rle, locs)
  rle.data.frame = DataFrame(rle.list, row.names=rownames(locs))
  return(rle.data.frame)
}

##' Make a RangedData from segments
##'
##' Starting from a data.frame of segments, like from CBS and segTable, organize as a RangedData. Label data "score",
##' so it can easily be made into various genome browser formats using rtracklayer.
##' @param segs data.frame, like from segment in DNAcopy or segTable
##' @return RangedData
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
##' and make a list of data.frames each like the result of CBS's
##' segment.  Note the loc.start and loc.stop will correspond
##' exactly to probe locations in locData and the input to
##' segs2RleDataFrame are not necessarily so.
##'
##' @param object Rle or list/DataFrame of Rle vectors
##' @param locs RangedData with rows corresponding to rows of df
##' @param sample.name character for single Rle optionally include "ID" column with this sample name
##' @return one or a list of data.frames with columns ID, chrom, loc.start, loc.end, num.mark, seg.mean
##' @export segTable
##' @examples
##'   data(genoset)
##'   seg.list = runCBS( lrr(baf.ds), locData(baf.ds), return.segs=TRUE )
##'   df = segs2RleDataFrame( seg.list, locData(baf.ds) )  # Loop segs2Rle on list of data.frames in seg.list
##'   assayDataElement( baf.ds, "lrr.segs" ) = df
##'   segTable( df, locData(baf.ds) )
##'   segTable( assayDataElement(baf.ds,"lrr.segs"), locData(baf.ds) )
##'   segTable( assayDataElement(baf.ds,"lrr.segs")[,1], locData(baf.ds), sampleNames(baf.ds)[1] )
##' @author Peter M. Haverty
setGeneric("segTable", function(object,...) standardGeneric("segTable"))
setMethod("segTable", signature(object="Rle"), function(object,locs,sample.name=NULL) {
  # All the time goes into start, end, and space, particularly the unlist.  Maybe faster by chr then cbind?
  chr.ind = chrIndices(locs)
  num.mark = unlist(aggregate(object, FUN=runLength, start=chr.ind[,"first"], end=chr.ind[,"last"]))
  seg.mean = unlist(aggregate(object, FUN=runValue, start=chr.ind[,"first"], end=chr.ind[,"last"]))

  loc.end.indices = cumsum(num.mark)
  loc.end = end(locs)[loc.end.indices]  # unlist here is the big time waster
  loc.start.indices = (loc.end.indices - num.mark) + 1L
  loc.start = start(locs)[loc.start.indices] # unlist here is the big time waster
  chrom = space(locs)[loc.start.indices]
  if (is.null(sample.name)) {
    sample.seg = data.frame(chrom = chrom, loc.start = loc.start, loc.end = loc.end, num.mark = num.mark, seg.mean = seg.mean, row.names=NULL)
  } else {
    sample.seg = data.frame(ID = sample.name, chrom = chrom, loc.start = loc.start, loc.end = loc.end, num.mark = num.mark, seg.mean = seg.mean, row.names=NULL)
  }
  return(sample.seg)
})

setMethod("segTable", signature(object="DataFrame"), function(object,locs) {
  segs = sapply( names(object),
    function(x) {
      temp.rle = object[[x]]
      sample.seg = segTable(temp.rle,locs,sample.name=x)
      return(sample.seg)
    },simplify=FALSE, USE.NAMES=TRUE)
  return(segs)
})

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
##' @param n.cores numeric, number of cores to ask multicore to use
##' @param smooth.region number of positions to left and right of individual positions to consider when smoothing single point outliers
##' @param outlier.SD.scale number of SD single points must exceed smooth.region to be considered an outlier
##' @param smooth.SD.scale floor used to reset single point outliers
##' @param trim fraction of sample to smooth
##' @param alpha pvalue cutoff for calling a breakpoint
##' @return data frame of segments from CBS
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
runCBS <- function(data, locs, return.segs=FALSE, n.cores=getOption("cores"), smooth.region=2, outlier.SD.scale=4, smooth.SD.scale=2, trim=0.025, alpha=0.001) {
  sample.name.list = colnames(data)
  names(sample.name.list) = sample.name.list
  loc.pos = as.numeric(pos(locs))
  loc.chr = chr(locs)
  
  # mclapply over samples. cbs can loop over the columns of data, but want to use multiple forks
  if (n.cores > 1 && is.loaded("mc_fork", PACKAGE="multicore")) {
    mcLapply <- get('mclapply', envir=getNamespace('multicore'))
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
      CNA.object <- CNA(temp.data[ok.indices], loc.chr[ok.indices], loc.pos[ok.indices], data.type = "logratio", sampleid = sample.name)
      smoothed.CNA.object <- smooth.CNA(CNA.object, smooth.region=smooth.region, outlier.SD.scale=outlier.SD.scale, smooth.SD.scale=smooth.SD.scale, trim=trim)
      segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=0, alpha=alpha)
      if (return.segs == TRUE) {
        return(segment.smoothed.CNA.object$output)
      } else {
        return(segs2Rle(segment.smoothed.CNA.object$output,locs))
      }
    })

  if (return.segs == TRUE) {
    return(segs)
  } else {
    return( DataFrame(segs, row.names=rownames(locs) ) )
  }
}

##' Find indices of features bounding a set of chromsome ranges/genes
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


##' Find indices of features bounding a set of chromsome ranges/genes
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
##' important differences from intervalBound, which uses findInterval: boundingIndices does not
##' check for NAs or unsorted data in the subject positions. Also, the positions are
##' kept as integer, where intervalBound (and findInterval) convert them to doubles. These
##' three once-per-call differences account for much of the speed improvement in boundingIndices.
##' These three differences are meant for position info coming from GenoSet objects
##' and boundingIndices2 is safer for general use. boundingIndices works on integer postions and
##' does not check that the positions are ordered. The starts and stops need not be sorted, but
##' it will be much faster if they are.
##'
##' @param starts integer vector of first base position of each query range
##' @param stops integer vector of last base position of each query range
##' @param positions Base positions in which to search
##' @param valid.indices logical, TRUE assures that the returned indices don't go off either end of the array, i.e. 0 becomes 1 and n+1 becomes n
##' @param offset integer, value to add to all returned indices. For the case where positions represents a portion of some larger array (e.g. a chr in a genome)
##' @param all.indices logical, return a list containing full sequence of indices for each query
##' @return integer matrix of 2 columms for start and stop index of range in data or a list of full sequences of indices for each query (see all.indices argument)
##' @seealso boundingIndices2
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
##' @seealso boundingIndices intervalBound
##' @export rangeSampleMeans
##' @examples
##'   data(genoset)
##'   my.genes = RangedData( ranges=IRanges(start=c(35e6,128e6),end=c(37e6,129e6),names=c("HER2","CMYC")), space=c("chr17","chr8"), universe="hg19")
##'   rangeSampleMeans( my.genes, baf.ds, "lrr" )
##' @author Peter M. Haverty
rangeSampleMeans <- function(query.rd, subject, assay.element) {
  ## Find feature bounds of each query in subject genoset, get feature data average for each sample

  if (! isGenomeOrder(locData(subject),strict=FALSE) ) {
    cat("Setting subject to genome order.")
    subject = toGenomeOrder(subject)
  }
  if (! isGenomeOrder(query.rd,strict=FALSE) ) {
    cat("Setting query to genome order.")
    query.rd = toGenomeOrder(query.rd)
  }

  chr.names = intersect( names(query.rd), names(locData(subject)) )  # All chrs in both sets
  chr.indices = chrIndices(subject)  # Bounds of each chr

  # Foreach chr
  subject.ranges = ranges(subject)
  subject.starts = start(subject.ranges)
  subject.stops = end(subject.ranges)
  query.ranges = ranges(query.rd)
  query.starts = start(query.ranges)
  query.stops = end(query.ranges)
  
  ranges.by.chr = lapply( chr.names, function(x) {  # lapply here so I can switch to mclapply later
    indices = boundingIndices( starts=query.starts[[x]], stops=query.stops[[x]], positions=subject.starts[[x]], offset=chr.indices[x,"offset"] ) # Find bounds of query in this chr
    rownames(indices) = names(query.ranges[[x]])
    return(indices)
  })
  all.indices = do.call(rbind,ranges.by.chr)

  # Temporary hack for DataFrame of Rle
  data.matrix = assayDataElement(subject,assay.element)

  if (class(data.matrix) == "DataFrame") {
    sample.vals = sapply( colnames(data.matrix), function(x) { rangeColMeans( all.indices, as.numeric(data.matrix[,x] )) }, USE.NAMES=TRUE, simplify=FALSE)
    range.means = do.call(cbind,sample.vals)
  } else {
    range.means = rangeColMeans( all.indices, data.matrix )
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

##' Check if a RangedData or GenoSet is in genome order
##'
##' Checks that rows in each chr are ordered by start.  If strict=TRUE, then chromsomes
##' must be in order specified by chrOrder.
##' 
##' @param ds RangedData or GenoSet
##' @param strict logical, should space/chromosome order be identical to that from chrOrder?
##' @return logical
##' @export isGenomeOrder
##' @examples
##'   data(genoset)
##'   isGenomeOrder( locData(genoset.ds) )
##' @author Peter M. Haverty
isGenomeOrder <- function(ds, strict=FALSE) {
  if (strict) {
    if ( ! all( names(ds) == chrOrder(names(ds) ) ) ) {
      return(FALSE)
    }
  }
  for (chr.name in names(ds)) { # Check each chr for ordered start
    if ( is.unsorted(start(ranges(ds))[[chr.name]]) ) {
      return(FALSE)
    }
  }
  return(TRUE)
}

##' @rdname genomeorder
setGeneric("toGenomeOrder", function(ds,...) standardGeneric("toGenomeOrder"))

##' Get indices to set a RangedData or GenoSet to genome order
##'
##' Returns a vector of idices to use in re-ordering a RangedData or
##' GenoSet to genome order. If strict=TRUE, then chromsomes must be in order specified by chrOrder.
##' 
##' @param ds RangedData or GenoSet
##' @param strict logical, should chromosomes be in order specified by chrOrder?
##' @return numeric vector of indices for re-ordering
##' @export toGenomeOrder
##' @examples
##'   data(genoset)
##'   toGenomeOrder( baf.ds, strict=TRUE )
##'   toGenomeOrder( baf.ds )
##'   toGenomeOrder( locData(baf.ds) )
##' @author Peter M. Haverty
##' @rdname genomeorder
setMethod("toGenomeOrder",signature=signature(ds="RangedData"),
          function(ds, strict=FALSE) {
            if (strict == TRUE) {
              ds = ds[ chrOrder(names(ds)) ]
            }
            return( ds[order(as.integer(space(ds)),start(ds)),,drop=FALSE] )
          }
        )
##' @rdname genomeorder
setMethod("toGenomeOrder", signature=signature(ds="GenoSet"),
          function(ds,strict=TRUE) {
            locData(ds) = toGenomeOrder(locData(ds),strict=strict) # locData<- fixes row ordering in ds
            return(ds)
          })

#############################
### Working with data on disk
#############################

##' Attach on-disk matrices into assayData
##'
##' GenoSet objects can hold big.matrix objects in their assayData slot environment.
##' After re-loading the GenoSet from disk, these objects will each need to be re-attached
##' to their on-disk component using their resource locators stored in their "desc"
##' attributes. This function checks each assayDataElement to see if it is an un-attached
##' big.matrix object, re-attaching if necessary. All other assayDataElements are left
##' untouched. In later releases this function will also handle other on-disk types,
##' like HDF5-based matrices.
##' @param object GenoSet
##' @return GenoSet
##' @export 
##' @author Peter M. Haverty \email{phaverty@@gene.com}
attachAssayDataElements <- function(object) {
  for( ad.name in assayDataElementNames(object)) {
    if ( is.big.matrix( assayDataElement(object,ad.name) ) && is.nil( assayDataElement(object,ad.name)@address ) ) {
      if (is.null(attr(assayDataElement(object,ad.name),"desc"))) {
        stop("Failed to attach assayDataElement",ad.name,". No 'desc' attribute.")
      } else {
        assayDataElement(object,ad.name)@address = attach.big.matrix( attr(assayDataElement(object,ad.name),"desc") )@address
      }
    }
  }
  return(object)
}

##' Load a GenoSet from a RData file
##'
##' Given a RData file with one object (a GenoSet or related object), load it, attach bigmatrix
##' objects as necessary, and return.
##' @param path character, path to RData file
##' @return GenoSet or related object (only object in RData file)
##' @examples
##' \dontrun{ ds = readGenoSet("/path/to/genoset.RData") }
##' @export 
##' @author Peter M. Haverty \email{phaverty@@gene.com}
readGenoSet <- function(path) {
  object = get(load(path)[1])
  if (!is(object,"eSet")) { stop("Loaded object is not an eSet or derived class.") }
  return( attachAssayDataElements(object) )
}

##' Make standard matrices in a GenoSet filebacked bigmatrix objects
##'
##' Make standard matrices in a GenoSet filebacked bigmatrix objects. Something like a
##' factor can be obtained using integer assayDataElements with a "levels" attribute. The
##' levels attribute will be maintained. Such objects will be stored as char on disk if
##' there are < 128 levels, and integer otherwise. "nlevels" and "levels" will work on
##' these objects as they only require the levels attribute. The "as.character" functionality
##' of a factor can be obtained like this: levels(assayDataElement(ds,"geno"))[ ds[1:5,1:5,"geno"] ]
##' for a GenoSet called "ds" with a factor-like element called "geno".
##' @param object GenoSet
##' @param prefix character, prefix for all bigmatrix related files
##' @param path character, directory to be created for all bigmatrix files, can be pre-existing.
##' @return GenoSet or related, updated copy of "object"
##' @examples
##' \dontrun{ ds = convertToBigMatrix(ds) }
##' @export 
##' @author Peter M. Haverty \email{phaverty@@gene.com}
convertToBigMatrix <- function(object,prefix="bigmat",path="bigmat") {
  orig.umask = Sys.umask()
  Sys.umask("002")
  dir.create(path,showWarnings=FALSE)
  path = normalizePath(path)
  for( ad.name in assayDataElementNames(object)) {

    if (is.matrix( assayDataElement(object,ad.name) ) ) {
      back.file = paste(prefix,ad.name,"bin",sep=".")
      desc.file = paste(prefix,ad.name,"desc",sep=".")

      if (is.integer( assayDataElement(object,ad.name))) {
        if (nlevels( assayDataElement(object,ad.name)) > 0) {  # Incoming "factor" from asFactorMatrix
          options(bigmemory.typecast.warning=FALSE)

          if ( nlevels(assayDataElement(object,ad.name)) > 127 ) {
            mat.type = "integer"
          } else {
            mat.type = "char"
          }
        } else {
          mat.type = "integer"
        }
      } else if (is.double( assayDataElement(object,ad.name))) {
        mat.type = "double"
      } else {
        # Other types not handled yet
        ### character matrices probably too big to unique, need to handle conversion to factor with asFactorMatrix
        next;
      }

      cat("Converting", ad.name, "to big.matrix ...\n")
      new.matrix = as.big.matrix(assayDataElement(object,ad.name),
        type=mat.type,
        backingfile=back.file,
        descriptorfile=desc.file,
        backingpath=path)
      attr(new.matrix,"desc") = file.path(path,desc.file)
      if (nlevels(assayDataElement(object,ad.name)) > 0) {
        attr( new.matrix, "levels") = levels(assayDataElement(object,ad.name))
        options(bigmemory.typecast.warning=TRUE)
      }
      assayDataElement(object,ad.name) = new.matrix
    }
  }
  Sys.umask(orig.umask)
  return(object)
}

##' Make factor matrix from character matrix
##'
##' Make factor matrix from character matrix for use with convertToBigMatrix.
##' Makes an integer matrix with levels since as.big.matrix would make a
##' factor matrix into a 1D object for some reason. Character matrices should
##' be converted to factors with explicit levels as huge matrices are likely
##' too big to unique.
##'
##' Caution: use asFactorMatrix on matrices already in an eSet.  The eSet constructor will
##' apparently wipe out the levels.
##'
##' @param object matrix of characters
##' @param levels character
##' @return factor with dimensions matching object
##' @export 
##' @author Peter M. Haverty \email{phaverty@@gene.com}
asFactorMatrix <- function(object, levels) {
  big.factor = as.integer(factor(object,levels=levels))
  attributes(big.factor) = attributes(object)
  levels(big.factor) = levels
  return(big.factor)
}

#geno.levels = paste( c(rep("A",4),rep("C",4),rep("G",4),rep("T",4)), rep(c("A","C","G","T"),4), sep="")
#ab.levels = c("AA","AB","BB")

##' Update "desc" attributes for big.matrix assayDataElement to new location
##'
##' Update "desc" attributes for big.matrix assayDataElement to new location. Assumes files have already
##' been moved on the filesystem. Assumes names of description and data files are the same.
##' 
##' @param ds eSet
##' @param new.bigmat.dir character, path to directory holding desc and data files
##' @return eSet
##' @export 
##' @author Peter M. Haverty \email{phaverty@@gene.com}
relocateAssayData <- function(ds, new.bigmat.dir) {
  for (ad.name in assayDataElementNames(ds)) {
    if (is.big.matrix( assayDataElement(ds,ad.name) )) {
      attr( assayDataElement(ds,ad.name), "desc" ) = file.path( new.bigmat.dir, basename( attr(assayDataElement(ds,ad.name),"desc") ) )
    }
  }
  return(ds)
}
