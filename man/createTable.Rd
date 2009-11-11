\name{createTable}
\Rdversion{1.3}
\alias{createTable}
\title{Function to create an output table}
\description{
This function reports the results from the Frequentist and Bayesian model for qmax and for q2. It also creates an output table with the results for all the thresholds in a csv format, so the user can select additional thresholds of interest.
}
\usage{
createTable(output.ratio, output.bay, dir = getwd(),q=NULL)
}
\arguments{
  \item{output.ratio}{ \code{The output object from the Frequentist model (ratio function)} }
  \item{output.bay}{ 	\code{The output object from the Bayesian model (baymod function)} }
  \item{dir}{ \code{Directory for storing the table} }
  \item{q}{\code{Additional thresholds in the form of a vector}}
}
\details{
To select a list of interesting features from the Bayesian model we suggest two decision rules in the paper: 1) the maximum of Median(R(q)) only for the subset of credibility intervals which do not include 1; 2) the largest threshold q for which the ratio R(q) il bigger than 2.

The first one is pointing out the strongest deviation from independence, whilst the second is the largest threshold where the number of features called in common at least doubles the number of features in common under independence. 
}
\value{
  \item{max }{The results of the R(qmax) statistic}
  \item{rule2 }{The results using the rule R(q) larger than 2 (see details)}
  \item{ruleq}{The results using additional thresholds}
}
\references{ 1. M.Blangiardo and S.Richardson (2007) Statistical tools for synthesizing lists of differentially expressed features in related experiments , Genome Biology, 8, R54  }
\author{ Alberto Cassese, Marta Blangiardo }
\examples{
data = simulation(n=500,GammaA=1,GammaB=1,r1=0.5,r2=0.8,
DEfirst=300,DEsecond=200,DEcommon=100)
Tq<- ratio(data=data$Pval)
Rq<- baymod(iter=100,output.ratio=Tq)
output.table <- createTable(output.ratio=Tq,output.bay=Rq)

}
