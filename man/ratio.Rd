\name{ratio}
\alias{ratio}
\title{Ratio Th between the observed features in common and the expected ones}
\description{
The function for each experiment calculates the ratio T(h) for each threshold h, using the list of p-values or the other measure used in the experiment to rank the features, (e.g. posterior probability, correlation).
}
\usage{
ratio(data, pvalue = TRUE, interval = 0.01,
name ="Distribution of T(h)",dir = getwd(),
dataname = "dataratio")
}
\arguments{
  \item{data}{ \code{Lists of pvalues to be compared} }
  \item{pvalue}{ \code{Indicate if the data are pvalues (TRUE) or posterior probability (FALSE). If they are posterior probability they are transformed in pvalues. If correlation is used, select p-value=FALSE} }
  \item{interval}{ \code{The interval between two threshold} }
  \item{name}{ \code{The name to be used in the plots} }
  \item{dir}{ \code{Directory for storing the plots} }
  \item{dataname}{ \code{The name of the file containing the data (Pvalue)} }
}
\details{
This function calculates the ratio T(h) of observed number of features in common between the lists vs the expected number under the hypothesis of independence for each threshold h. The expected numbers are calculated as the product among the marginals divided by \eqn{(numbers of features)^(number of lists-1)}. T(hmax) identifies the maximum of the statistic T(h) and it is shown on the plot. 
}
\value{
Returns a plot with the distribution of T(h) showing where T(hmax) and hmax are located. The same plot is also saved in the directory specified by the user. It returns also an object of class list with the ratio, the thresholds and other attributes. In particular:
  \item{h }{Threshold corresponding to T(h) values}
  \item{DE }{Differentially expressed features in each experiment}
  \item{ratios }{Vector or T(h) values for each threshold}
  \item{Common }{Features in common corresponding to the T(h) values}
  \item{interval }{Interval on the p-value scale defined by the user (default is 0.01)}
  \item{name }{Names to be used in the plots (defined by the user)}
  \item{pvalue }{Logical: TRUE if the measures used for the analysis are p-value, FALSE if they are posterior probabilities}
  \item{dataname }{The name of the file where the data has been saved}
}
\references{ Stone et al.(1988), Investigations of excess environmental risks around putative sources: statistical problems and a proposed test,Statistics in Medicine, 7, 649-660.

M.Blangiardo and S.Richardson (2007) Statistical tools for synthesizing lists of differentially expressed features in related experiments, Genome Biology, 8, R54. }
\author{ Alberto Cassese, Marta Blangiardo }

\examples{
data = simulation(n=500,GammaA=1,GammaB=1,r1=0.5,r2=0.8,
DEfirst=300,DEsecond=200,DEcommon=100)
Th<- ratio(data=data$Pval)

}
