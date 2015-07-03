\name{designCount}
\alias{designCount}
\title{ Internal function }
\description{
This is an internal function. The user should not use the function directly.
}
\usage{
designCount(array, design)
}
\arguments{
  \item{array}{ \code{Data matrix} }
  \item{design}{ \code{Design matrix} }
}
\author{ Alberto Cassese, Marta Blangiardo }
\examples{

## The function is currently defined as
function(array,design) {

sum1intersects <- function(c1,c2) return(all(c1 == c2))
res <- vector(mode="numeric",length=(nrow(design)-1))

for(i in 2:nrow(design)){
res[i-1] <- sum(apply(array,1,sum1intersects,design[i,]))
  }
return(res)
  }
}
