\name{RleDataFrame-class}
\Rdversion{1.1}
\docType{class}
\alias{RleDataFrame-class}
\alias{colMeans,RleDataFrame-method}
\alias{colSums,RleDataFrame-method}
\alias{rowMeans,RleDataFrame-method}
\alias{rowSums,RleDataFrame-method}
\alias{show,RleDataFrame-method}

\title{Class \code{"RleDataFrame"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("RleDataFrame", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{rownames}:}{Object of class \code{"characterORNULL"} ~~ }
    \item{\code{nrows}:}{Object of class \code{"integer"} ~~ }
    \item{\code{elementType}:}{Object of class \code{"character"} ~~ }
    \item{\code{elementMetadata}:}{Object of class \code{"DataTableORNULL"} ~~ }
    \item{\code{metadata}:}{Object of class \code{"list"} ~~ }
    \item{\code{listData}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{SimpleRleList}"}, directly.
Class \code{"\linkS4class{DataFrame}"}, directly.
Class \code{"\linkS4class{RleList}"}, by class "SimpleRleList", distance 2.
Class \code{"\linkS4class{SimpleAtomicList}"}, by class "SimpleRleList", distance 2.
Class \code{"\linkS4class{DataTable}"}, by class "DataFrame", distance 2.
Class \code{"\linkS4class{AtomicList}"}, by class "SimpleRleList", distance 3.
Class \code{"\linkS4class{SimpleList}"}, by class "SimpleRleList", distance 3.
Class \code{"\linkS4class{DataTableORNULL}"}, by class "DataFrame", distance 3.
Class \code{"\linkS4class{List}"}, by class "DataFrame", distance 3.
Class \code{"\linkS4class{Vector}"}, by class "DataFrame", distance 4.
Class \code{"\linkS4class{Annotated}"}, by class "DataFrame", distance 5.
}
\section{Methods}{
  \describe{
    \item{colMeans}{\code{signature(x = "RleDataFrame")}: ... }
    \item{colSums}{\code{signature(x = "RleDataFrame")}: ... }
    \item{rowMeans}{\code{signature(x = "RleDataFrame")}: ... }
    \item{rowSums}{\code{signature(x = "RleDataFrame")}: ... }
    \item{show}{\code{signature(object = "RleDataFrame")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("RleDataFrame")
}
\keyword{classes}
