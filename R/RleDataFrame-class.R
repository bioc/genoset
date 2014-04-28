##### Additional methods for DataFrame class to allow use of DataTable of Rle vectors as assayData element

##' @export
setClass("RleDataFrame",
         representation(
           rownames = "characterORNULL",
           nrows = "integer"
           ),
         prototype(rownames = NULL,
                   nrows = 0L,
                   listData = structure(list(),  names = character())),
         contains = c("SimpleRleList", "DataFrame")
)

##' @export
RleDataFrame <- function(..., row.names=NULL, check.names=TRUE) {
  x = DataFrame(..., row.names=row.names, check.names=check.names)
  return( new("RleDataFrame", listData=x@listData, rownames=rownames(x), nrows=nrow(x)) )
}

##' @export
setMethod("show", "RleDataFrame",
          function(object) {
            message(sprintf("RleDataFrame with %i rows and %i columns\n", nrow(object), ncol(object)))
            if (is.null(rownames(object))) {
              message("rownames: NULL\n")
            } else {
              message(sprintf("rownames: %s\n", paste(head(rownames(object)), collapse=", ")))
            }
            show(object@listData)
          })
