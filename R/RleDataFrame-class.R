##### Additional methods for DataFrame class to allow use of DataTable of Rle vectors as assayData element

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

RleDataFrame <- function(..., row.names=NULL, check.names=TRUE) {
  x = DataFrame(..., row.names=row.names, check.names=check.names)
  return( new("RleDataFrame", listData=x@listData, rownames=rownames(x), nrows=nrow(x)) )
}

setMethod("show", "RleDataFrame",
          function(object) {
            message(sprintf("RleDataFrame with %i rows and %i columns\n", nrow(object), ncol(object)))
            show(object@listData)
          })
