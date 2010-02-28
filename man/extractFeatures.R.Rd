\name{extractFeatures.R}
\Rdversion{1.4}
\alias{extractFeatures.R}
\title{Extracting the lists of features of interest}
\description{
The function returns the list of features in common using the two suggested rules hmax and h2 (Bayesian model) and additional ones defined by the user.
}
\usage{
extractFeatures.R(output.ratio, output.bay, feat.names, h = NULL)
}
\arguments{
  \item{output.ratio}{ \code{The output object from the Frequentist model (ratio function)} }
  \item{output.bay}{ \code{The output object from the Bayesian model (baymod function)} }
  \item{feat.names}{ \code{Names of the features (e.g Affy ID for genes)} }
  \item{h}{ \code{Additional thresholds in the form of a vector to select a list of features in common. If it is NULL only the Rmax and rule2 are used to select the lists of features of interest} }
}
\details{
To select a list of interesting features from the Bayesian model we suggest two decision rules in the paper: 1) the maximum of Median(R(h)) only for the subset of credibility intervals which do not include 1; 2) the largest threshold h for which the ratio R(h) is bigger than 2.

The first one is pointing out the strongest deviation from independence, whilst the second is the largest threshold where the number of features called in common at least doubles the number of features in common under independence. The user can define additional thresholds of interest and obtain the list of associated features.
}
\value{
The function returns an object of the class list. Each element is a matrix where the first column contains the name of the features while the other columns contain the p-values* from the experiments. It also saves a .csv file with the same information.

* instead of the p-values any other measure used to rank the features in the experiments can be used
  \item{max }{The list of features of interest selected on the basis of the threshold associated to R(hmax)}
  \item{rule2 }{The list of features of interest selected on the basis of the threshold associated to R(h2)}
  \item{User}{The list of features of interest selected on the basis of the additional thresholds selected by the user}
}
\references{ 1. M.Blangiardo and S.Richardson (2007) Statistical tools for synthesizing lists of differentially expressed features in related experiments, Genome Biology, 8, R54 }
\author{ Alberto Cassese, Marta Blangiardo  }
\examples{
data = simulation(n=500,GammaA=1,GammaB=1,r1=0.5,r2=0.8,
DEfirst=300,DEsecond=200,DEcommon=100)
Th<- ratio(data=data$Pval)
Rh<- baymod(iter=100,output.ratio=Th)
feat.names = data$names
feat.lists <- extractFeatures.R(output.ratio=Th,output.bay=Rh,
feat.names=feat.names,h=NULL)

}
