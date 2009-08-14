\name{designMatrix}
\alias{designMatrix}
\title{ Internal function}
\description{
This is an internal function. The user should not use the function directly.
}
\usage{
designMatrix(lists)
}
\arguments{
  \item{lists}{ \code{The number of lists to be compared} }
}
\author{ Marta Blangiardo, Alberto Cassese }
\examples{

## The function is currently defined as
function(lists){
rows = 2^(lists)
ncycles = rows
x = matrix(0,rows,lists)
for (k in 1:lists){   
settings = c(0,1)   
ncycles = ncycles/2   
nreps = rows/(2*ncycles)   
settings = matrix(rep(settings,nreps),nreps,length(settings),byrow=TRUE)
settings = as.vector(settings) #impila in un vettore settings, una colonna sotto l'altra
settings = matrix(rep(settings,ncycles),length(settings),ncycles)
x[,lists-k+1] = as.vector(settings)
  }
return(x)
  }
}
