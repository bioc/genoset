#############################
### Working with data on disk
#############################

##' @include genoset-class.R
{}

# Allow eSet constructor to make featureNames from a big.matrix as if it were a matrix
setMethod("annotatedDataFrameFrom",
          signature(object="big.matrix"),
          Biobase:::annotatedDataFrameFromMatrix)

##' Get assayDataElement, attaching on-disk resource if necessary
##'
##' Get assayDataElement, attaching on-disk resource if necessary
##' @param object eSet
##' @param elt character
##' @return assayDataElement, matrix, DataFrame, or the like
##' @export 
##' @author Peter M. Haverty \email{phaverty@@gene.com}
assayDataElement <- function(object, elt) {
  if ( is.big.matrix( assayData(object)[[elt]] ) && is.nil( assayData(object)[[elt]]@address ) ) {
    cat("Attaching big.matrix class assayData elements to their on-disk components ...\n")
    attachAssayDataElements(assayData(object))
  }
  return(assayData(object)[[elt]])
}

##' Set assayDataElement, attaching on-disk resource if necessary
##'
##' Set assayDataElement, attaching on-disk resource if necessary
##' @param object eSet
##' @param elt character, assayDataElement name
##' @param value input data to assayDataElement
##' @return eSet
##' @export 'assayDataElement<-'
##' @usage assayDataElement(object, elt) <- value
##' @author Peter M. Haverty \email{phaverty@@gene.com}
'assayDataElement<-' <- function(object, elt, value) {
  if ( is.big.matrix( assayData(object)[[elt]] ) ) {
    if ( file.access(attr(object[,,elt],"desc"),2) < 0 ) {
      stop("You do not have write permission on the 'desc' file for this big.matrix class assayDataElement")
    }
  }
  return(assayDataElementReplace(object, elt, value))
}

##' Attach on-disk matrices into assayData
##'
##' GenoSet objects can hold big.matrix objects in their assayData slot environment.
##' After re-loading the GenoSet from disk, these objects will each need to be re-attached
##' to their on-disk component using their resource locators stored in their "desc"
##' attributes. This function checks each assayDataElement to see if it is an un-attached
##' big.matrix object, re-attaching if necessary. All other assayDataElements are left
##' untouched. In later releases this function will also handle other on-disk types,
##' like HDF5-based matrices.
##'
##' *** Intentional side-effects *** Environment type assayData objects, even
##' "lockedEnvironment" objects, will be updated in place (same pointer). This allows
##' for functions trying to access assayDataElements to attach before access, rather
##' than crashing R.
##' 
##' @param object eSet
##' @return assayData in storage mode of input assayData, invisibly.  Re-assignment back
##' original eSet only necessary if using a list type assayData.
##' @export 
##' @author Peter M. Haverty \email{phaverty@@gene.com}
attachAssayDataElements <- function(object) {
  # Most of the time goes into "dget", which reads and parses the description from file.  Could cache those ...
  aData = assayData(object)
  storage.mode <- storageMode(aData)
  
  for( ad.name in assayDataElementNames(aData)) {
    if ( is.big.matrix( aData[[ad.name]] ) && is.nil( aData[[ad.name]]@address ) ) {
      if (is.null(attr(aData[[ad.name]],"desc"))) {
        stop("Failed to attach assayDataElement",ad.name,". No 'desc' attribute.")
      } else {
        if (storage.mode == "lockedEnvironment") {
          unlockBinding(ad.name,aData)
        }
        aData[[ad.name]]@address = attach.big.matrix( attr(aData[[ad.name]],"desc") )@address
      }
    }
  }

  if (storage.mode == "lockedEnvironment") { Biobase:::assayDataEnvLock(aData) }
  return(invisible(aData))
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
  attachAssayDataElements(object)
  return( object )
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
