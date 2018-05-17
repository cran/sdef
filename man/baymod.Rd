\name{baymod}
\alias{baymod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian model for the ratio of observed to expected probability of features to be in common
}
\description{
The function specifies a Bayesian model for the ratio of observed to expected probability of features to be in common. A multinomial distribution is specified on the probabilities of being significant in any combination of the experiments (e.g. if two experiments are considered, the probability of being significant in none, one and two experiments is specified) and a prior distribution is put on their parameters. The quantity of interest is the ratio of the probability that a feature is in common, to the probability that a feature is in common by chance, called R(h).
}
\usage{
baymod(output.ratio, iter = 1000, dir = getwd(), conf = 95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output.ratio}{\code{The output object from the ratio function}
%%     ~~Describe \code{output.ratio} here~~
}
  \item{iter}{\code{Number of iterations to be performed}
%%     ~~Describe \code{iter} here~~
}
  \item{dir}{\code{Directory for storing the plots}
%%     ~~Describe \code{dir} here~~
}
  \item{conf}{\code{Size of Confidence Interval}
%%     ~~Describe \code{conf} here~~
}
}
\details{
It returns an object of class list with the ratio R(h) for each threshold and its quantiles specified by conf. R(h) is significant if its CI does not include 1. We consider two rules for selecting the list of genes of interest: 1) hmax is the maximum of Median(R(h)) only for the subset of credibility intervals which do not include 1; 2) h2 is the largest threshold where the number of features called in common at least doubles the number of features in common under independence (where R(h) larger than 2).

The function returns also a plot of the credibility interval for each threshold. The same plot is also saved in the
directory specified by the user.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A matrix and a plot with the quantiles of R(h) identified by \emph{conf} for each p-value threshold.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
 M.Blangiardo and S.Richardson (2007) Statistical tools for synthesizing lists of differentially expressed features in related experiments, Genome Biology, 8, R54.
 %% ~put references to the literature/web site here ~
}
\author{
Alberto Cassese, Marta Blangiardo
%%  ~~who you are~~
}
\examples{
data = simulation(n=500,GammaA=1,GammaB=1,r1=0.5,r2=0.8,DEfirst=300,
DEsecond=200,DEcommon=100)
Th<- ratio(data=data$Pval)
Rh<- baymod(iter=100,output.ratio=Th)
}

