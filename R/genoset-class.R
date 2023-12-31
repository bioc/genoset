###### Class definition for GenoSet, which will extend RangedSummarizedExperiment

##' GenoSet: An eSet for data with genome locations
##'
##' Load, manipulate, and plot copynumber and BAF data.
##'
##' @docType package
##' @name genoset-package
##' @aliases genoset genoset-package
##' @seealso genoset-datasets
##'
##' @importClassesFrom GenomicRanges GRanges GenomicRanges DelegatingGenomicRanges GNCList GPos
##' @importClassesFrom SummarizedExperiment SummarizedExperiment RangedSummarizedExperiment
##'
##' @importMethodsFrom GenomicRanges names 'names<-' length width
##' @importMethodsFrom IRanges as.data.frame as.list as.matrix cbind colnames 'colnames<-' end findOverlaps gsub
##' @importMethodsFrom IRanges intersect lapply mean nrow order ranges rownames
##'
##' @importFrom graphics abline axis axTicks box mtext plot.new plot.window points segments
##' @importFrom IRanges IRanges '%over%' Views RleList PartitioningByEnd
##' @importFrom GenomicRanges GRanges
##'
##' @import methods
##' @import BiocGenerics
##' @import S4Vectors
##' @import GenomeInfoDb
##' @import SummarizedExperiment
##' @useDynLib genoset, .registration=TRUE
NULL

############### Class GenoSet

##' @exportClass GenoSet
setClass("GenoSet", contains = "RangedSummarizedExperiment")

##' @exportClass GenoSetOrGenomicRanges
setClassUnion("GenoSetOrGenomicRanges", c("GenoSet", "GenomicRanges"))

##' Create a GenoSet object
##'
##' This function is the preferred method for creating a new GenoSet object. Currently,
##' a GenoSet is simply a RangedSummarizedExperiment with some API changes and extra
##' methods. Therefore, a GenoSet must always have a rowRanges.
##'
##' locations. Rownames are required to match featureNames.
##' @param rowRanges GenomicRanges, not a GenomicRangesList
##' @param assays list, SimpleList or matrix-like object
##' @param colData a data.frame or DataFrame of sample metadata with rownames matching the colnames of the matrices in assays
##' @param metadata a list of any other data you want to attach to the GenoSet object
##' @param x A GenoSet
##' @return A GenoSet object
##' @examples
##' test.sample.names = LETTERS[11:13]
##' probe.names = letters[1:10]
##' assays=list(matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)))
##' rowRanges=GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep('chr1',4),rep('chr3',2),rep('chrX',4)))
##' colData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5])))
##' rse=SummarizedExperiment(rowRanges=rowRanges,assays=assays,colData=colData,metadata=metadata)
##' gs = GenoSet(rowRanges, assays, colData)
##' @export GenoSet
##' @family GenoSet
##' @rdname genoset-methods
GenoSet <- function(rowRanges, assays, colData, metadata = list()) {
    if (!is(rowRanges, "GenomicRanges")) {
        stop("'rowRanges' must be a 'GenomicRanges' or one of its subclasses.")
    }
    if (!is(colData, "DataFrame")) {
        colData = as(colData, "DataFrame")
    }
    if (length(assays) > 0) {
        an = colnames(assays[[1]])
        stopifnot(!is.null(an) && identical(an,rownames(colData)))
    }

    rse = SummarizedExperiment(assays = assays, rowRanges = rowRanges, colData = colData,
        metadata = metadata)
    new("GenoSet", rse)
}

############# Sub-setters

##' Subset a GenoSet
##'
##' Subset a GenoSet
##' @exportMethod '['
##' @exportMethod '[<-'
##' @param x GenoSet
##' @param i character, GRanges, logical, integer
##' @param j character, logical, integer
##' @param k character or integer
##' @param withDimnames scalar logical, put dimnames on returned assay?
##' @param drop logical drop levels of space factor?
##' @param ... additional subsetting args
##' @examples
##'   data(genoset,package='genoset')
##'   genoset.ds[1:5,2:3]  # first five probes and samples 2 and 3
##'   genoset.ds[ , 'K']  # Sample called K
##'   gr = GRanges(ranges=IRanges(start=seq(from=15e6,by=1e6,length=7),width=1,names=letters[8:14]),seqnames=rep('chr17',7))
##'   genoset.ds[ gr, 'K' ]  # sample K and probes overlapping those in rd, which overlap specifed ranges on chr17
##' @rdname genoset-subset
setMethod("[", signature = signature(x = "GenoSet", i = "ANY"), function(x, i, j,
    k, ..., withDimnames = TRUE, drop = FALSE) {
    if (!missing(i) && is(i, "GenomicRanges")) {
        i = unlist(rowRanges(x) %over% i)
    }
    if (missing(k)) {
        rval = callNextMethod(x, i, j, ..., drop = drop)
        return(rval)
    } else {
        if (missing(i) && missing(j)) {
            return(assay(x, k))
        } else {
            if (missing(i)) {
                return(assay(x, k)[, j])
            } else if (missing(j)) {
                return(assay(x, k)[i, ])
            } else {
                return(assay(x, k)[i, j])
            }
        }
    }
})

##' @param value incoming data for assay 'k', rows 'i' and cols 'j'
##' @rdname genoset-subset
setMethod("[<-", signature = signature(x = "GenoSet", i = "ANY"), function(x, i,
    j, k, value) {
    if (missing(k)) {
        stop("Must specify k to replace data in the GenoSet")
    }
    if (!missing(i) && is(i, "GenomicRanges")) {
        i = unlist(rowRanges(x) %over% i)
    }
    if (missing(i) && missing(j)) {
        if (!all(colnames(x) == colnames(value)) || !all(rownames(x) == rownames(value))) {
            stop("Dimnames for incoming assay must match this genoset.\n")
        }
        assay(x, k) = value
    } else {
        if (missing(i)) {

            assay(x, k)[, j] = value
        } else if (missing(j)) {
            assay(x, k)[i, ] = value
        } else {
            assay(x, k)[i, j] = value
        }
    }
    return(x)
})

######################## Get genome information

##' Chromosome name for each feature
##'
##' Get chromosome name for each feature.  Returns character.
##' @param object GRanges GenoSet
##' @return character vector of chromosome positions for each feature
##' @examples
##'   data(genoset,package='genoset')
##'   chr(genoset.ds)  # c('chr1','chr1','chr1','chr1','chr3','chr3','chrX','chrX','chrX','chrX')
##'   chr(rowRanges(genoset.ds))  # The same
##' @export chr
##' @rdname chr-methods
setGeneric("chr", function(object) standardGeneric("chr"))
##' @rdname chr-methods
setMethod("chr", "GenoSet", function(object) {
    return(chr(slot(object, "rowRanges")))
})
##' @rdname chr-methods
setMethod("chr", "GenomicRanges", function(object) {
    return(as.character(seqnames(object)))
})

##' Chromosome position of features
##'
##' Get chromosome position of features/ranges. Defined as floor of mean of start and end.
##' @param x GRanges GenoSet
##' @return numeric vector of feature positions within a chromosome
##' @export pos
##' @examples
##'   data(genoset,package='genoset')
##'   pos(genoset.ds)  # 1:10
##'   pos(rowRanges(genoset.ds))  # The same
##' @rdname pos-methods
##' @importFrom GenomicRanges pos
setMethod("pos", "GenoSetOrGenomicRanges", function(x) {
    return(start(x) + (width(x) - 1L)%/%2L)
})

##' Get list of unique chromosome names
##'
##' Get list of unique chromosome names
##'
##' @param object GenomicRanges or GenoSet
##' @param value return value of chrNames
##' @return character vector with names of chromosomes
##' @examples
##'   data(genoset,package='genoset')
##'   chrNames(genoset.ds) # c('chr1','chr3','chrX')
##'   chrNames(rowRanges(genoset.ds))  # The same
##'   chrNames(genoset.ds) = sub('^chr','',chrNames(genoset.ds))
##' @rdname chrNames-methods
##' @export 'chrNames'
##' @export 'chrNames<-'
setGeneric("chrNames", function(object) standardGeneric("chrNames"))

##' @rdname chrNames-methods
setMethod("chrNames", signature(object = "GenoSet"), function(object) {
    chrNames(rowRanges(object))
})

##' @rdname chrNames-methods
setMethod("chrNames", signature(object = "GenomicRanges"), function(object) {
    as.character(unique(seqnames(object)))
})

##' @rdname chrNames-methods
##' @export 'chrNames<-'
setGeneric("chrNames<-", function(object, value) standardGeneric("chrNames<-"))

##' @rdname chrNames-methods
setMethod("chrNames<-", signature(object = "GenoSet"), function(object, value) {
    chrNames(rowRanges(object)) = value
    return(object)
})

##' @rdname chrNames-methods
setMethod("chrNames<-", signature(object = "GenomicRanges"), function(object, value) {
    seqlevels(object) = value
    return(object)
})

##' Get chromosome start and stop positions
##'
##' Provides a matrix of start, stop and offset, in base numbers for each chromosome.
##' @param object A GenoSet object or similar
##' @return list with start and stop position, by ordered chr
##' @export chrInfo
##' @examples
##'   data(genoset,package='genoset')
##'   chrInfo(genoset.ds)
##'   chrInfo(rowRanges(genoset.ds))  # The same
##' @rdname chrInfo-methods
setGeneric("chrInfo", function(object) standardGeneric("chrInfo"))

##' @rdname chrInfo-methods
setMethod("chrInfo", signature(object = "GenoSetOrGenomicRanges"), function(object) {
    # Get max end value for each chr
    if (is(object, "GenomicRanges") && !any(is.na(seqlengths(object)))) {
        max.val = seqlengths(object)
    } else {
        max.val = max(relist(end(object), chrPartitioning(object)))
    }
    if (length(max.val) == 1) {
        names(max.val) = chrNames(object)
    } else {
        max.val = max.val[chrOrder(chrNames(object))]
    }

    chr.info = matrix(ncol = 3, nrow = length(max.val), dimnames = list(names(max.val),
        c("start", "stop", "offset")))
    chr.info[, "stop"] = cumsum(as.numeric(max.val))
    chr.info[, "offset"] = c(0, chr.info[-nrow(chr.info), "stop"])
    chr.info[, "start"] = chr.info[, "offset"] + 1

    return(chr.info)
})

##' Partitioning by Chromosome
##'
##' Get indices of first and last element in each chromosome.
##' @param object GenoSet or GenomicRanges
##' @return PartitioningByEnd
##' @export
chrPartitioning <- function(object) {
    rle = Rle(seqnames(object))  # Redundant in some cases
    ends = structure(cumsum(runLength(rle)), names = as.character(runValue(rle)))
    return(PartitioningByEnd(ends))
}

##' Get a matrix of first and last index of features in each chromosome
##'
##' Sometimes it is handy to know the first and last index for each chr.
##' This is like chrInfo but for feature indices rather than chromosome
##' locations. If chr is specified, the function will return a sequence
##' of integers representing the row indices of features on that chromosome.
##'
##' @param object GenoSet or GRanges
##' @param chr character, specific chromosome name
##' @return data.frame with 'first' and 'last' columns
##' @export chrIndices
##' @examples
##'   data(genoset,package='genoset')
##'   chrIndices(genoset.ds)
##'   chrIndices(rowRanges(genoset.ds))  # The same
##' @rdname chrIndices-methods
setGeneric("chrIndices", function(object, chr = NULL) standardGeneric("chrIndices"))

##' @rdname chrIndices-methods
setMethod("chrIndices", signature(object = "GenoSetOrGenomicRanges"), function(object,
    chr = NULL) {
    partitions = chrPartitioning(object)
    chr.first = start(partitions)
    chr.last = end(partitions)
    chr.info = matrix(c(chr.first, chr.last, chr.first - 1), ncol = 3, nrow = length(chr.first),
        dimnames = list(names(partitions), c("first", "last", "offset")))
    if (!is.null(chr)) {
        if (!chr %in% rownames(chr.info)) {
            stop("Must specify a valid chromosome name in chrIndices.\n")
        }
        return(seq.int(chr.info[chr, "first"], chr.info[chr, "last"]))
    } else {
        return(chr.info)
    }
})

##' Get base positions of features in genome-scale units
##'
##' Get base positions of array features in bases counting from the start of the
##' genome. Chromosomes are ordered numerically, when possible, then lexically.
##' @param object A GenoSet object or a GenomicRanges object
##' @return numeric position of each feature in whole genome units, in original order
##' @examples
##'   data(genoset,package='genoset')
##'   head(genoPos(genoset.ds))
##'   head(genoPos(rowRanges(genoset.ds)))  # The same
##' @export genoPos
##' @rdname genoPos-methods
setGeneric("genoPos", function(object) standardGeneric("genoPos"))

##' @rdname genoPos-methods
setMethod("genoPos", signature(object = "GenoSetOrGenomicRanges"), function(object) {

    # For single chr objects, just return pos
    if (length(chrNames(object)) == 1) {
        return(pos(object))
    }

    ### Add offset to pos by chr
    offset = chrInfo(object)[, "offset"]
    genopos = pos(object) + unlist(offset[chr(object)])

    return(genopos)
})

##' @rdname genoset-methods
##' @export lengths
setMethod("lengths", signature(x = "GenoSet"), function(x) {
    lengths(rowRanges(x))
})

##' GenomicRanges API Additions
##'
##' I have extended the API for GenomicRanges a bit so that genoset
##' and GenomicRanges can have the same API, at least as far as
##' genome location based features go.
##' @param x A GenomicRanges
##' @rdname genomicranges-methods
##' @export nrow
setMethod("nrow", signature(x = "GenomicRanges"), function(x) {
    length(x)
})
