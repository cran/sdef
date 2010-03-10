\name{Tmc}
\Rdversion{1.5}
\alias{Tmc}
\title{Empirical null distribution of max T(h)}
\description{
The function uses Monte Carlo permutations to calculate the empirical distribution of max T(h)=T(hmax) under the null hypothesis of independence among the experiments. An empirical p-value is calculated to evaluate where T(hmax) is located under the null distribution.
}
\usage{
Tmc(iter = 1000, output.ratio)
}
\arguments{
  \item{iter}{ \code{Number of iteration to be performed} }
  \item{output.ratio}{ \code{The output object from the ratio function} }
}
\details{
This function uses Monte Carlo permutations to calculate the empirical distribution of the maximum of T(h) (i.e. T(hmax)) under the null hypothesis of independence among the experiments. While the p-values* for the first list are fixed, the ones for the other lists are independently permutate B times. In this way, any relationship among the lists is destroyed. At each permutation b (b varies from 1 to B) a Tb(h) is calculated for each h and a maximum statistic Tb(hmax) is returned; its distribution represents the null distribution of T(hmax) under the condition of independence. The relative frequency of Tb(hmax) larger than T(hmax) identifies the p-value: it returns the proportion of Tb(hmax) from permuted dataset greater than the observed one (so indicates where the observed T(hmax) is located under the null distribution).

* instead of the p-values any other measure used to rank the features in the experiments can be used 
}
\value{
Returns the empirical pvalue from testing T(hmax) and a plot of the Tb(hmax) distribution. The same plot is also saved in the directory specified by the user.
}
\references{ Stone et al.(1988), Investigations of excess environmental risks around putative sources: statistical problems and a proposed test,Statistics in Medicine, 7, 649-660.

M.Blangiardo and S.Richardson (2007) Statistical tools for synthesizing lists of differentially expressed features in related experiments, Genome Biology, 8, R54. }
\author{ Alberto Cassese, Marta Blangiardo }
\examples{
data = simulation(n=500,GammaA=1,GammaB=1,
r1=0.5,r2=0.8,DEfirst=300,DEsecond=200,
DEcommon=100)
Th<- ratio(data=data$Pval)
MC<- Tmc(iter=100,output.ratio=Th)

}
