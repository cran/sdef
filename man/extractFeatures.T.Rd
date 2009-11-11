\name{extractFeatures.T}
\Rdversion{1.3}
\alias{extractFeatures.T}
\title{Extracting the lists of features of interest}
\description{
The function returns the list of features in common using the qmax rule (Frequentist model).
}
\usage{
extractFeatures.T(output.ratio, feat.names)
}
\arguments{
  \item{output.ratio}{ \code{The output object from the Frequentist model (ratio function)} }
  \item{feat.names}{ \code{names of the features (e.g Affy ID for genes)} }
}
\details{
To select a list of interesting features from the frequentist model we suggest two decision rules in the paper: 1) the maximum of Median(T(q); 2) the largest threshold q for which the ratio T(q) is bigger than 2.

The first one is pointing out the strongest deviation from independence, whilst the second is the largest threshold where the number of features called in common at least doubles the number of features in common under independence.
}
\value{
The function returns an object of the class list. Each element is a matrix where the first column contains the name of the features while the other columns contain the p-values* from the experiments. It also saves a .csv file with the same information.

*instead of the p-values any other measure used to rank the features in the experiments can be used.
  \item{max }{The list of features in common selected on the basis of the threshold associated to T(qmax)}
}
\references{ 1. M.Blangiardo and S.Richardson (2007) Statistical tools for synthesizing lists of differentially expressed features in related experiments, Genome Biology, 8, R54 }
\author{ Alberto Cassese, Marta Blangiardo  }
\examples{
data = simulation(n=500,GammaA=1,GammaB=1,r1=0.5,r2=0.8,
DEfirst=300,DEsecond=200,DEcommon=100)
Tq<- ratio(data=data$Pval)
feat.names = data$names
feat.lists.T <- extractFeatures.T(output.ratio=Tq,
feat.names=feat.names)

}

