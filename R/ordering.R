#############################################################################
#####  Functions related to the order of chromosomes and range features #####
#############################################################################


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
            if ( any(duplicated(runValue(seqnames(ds)))) ) { return(FALSE) }
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