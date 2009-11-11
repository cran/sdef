\name{simulation.indep}
\Rdversion{1.3}
\alias{simulation.indep}
\title{ Simulate p-values for two indipendent experiments }
\description{
The function simulate two vectors of p-values using the procedure described in Hwang et al. for independent experiments
}
\usage{
simulation.indep(n, GammaA = 2, GammaB = 2, epsilonM = 0,
epsilonSD = 1, r1, r2, DEfirst, DEsecond)
}
\arguments{
  \item{n}{ \code{Number of features to be simulated} }
  \item{GammaA}{ \code{Parameter of the Gamma distribution} }
  \item{GammaB}{ \code{Parameter of the Gamma distribution} }
  \item{epsilonM}{ \code{Parameter of the Gaussian noise} }
  \item{epsilonSD}{ \code{Parameter of the Gaussian noise} }
  \item{r1}{ \code{Additional experiment-specific noise} }
  \item{r2}{ \code{Additional experiment-specific noise} }
  \item{DEfirst}{ \code{Number of DE features in the first experiment} }
  \item{DEsecond}{ \code{Number of DE features in in the second experiment} }
}
\details{
Considering two experiments (k=1,2), each of them with two classes, and n genes, for each gene we simulate a true difference between the classes delta(g), drawn from a Gamma distribution with random sign. The true difference delta(g) is 0 if the gene is not differentially expressed. We then add two normal random noise components, r[k] that act as experiment specific components and epsilon(gk), that is the gene-experiment components. The former is assigned deterministically, whilst the latter is drawn from a standard Gaussian distribution. So the log fold change (FC(gk)) is the sum of all these components for each gene and experiment. We divide the n genes in three groups: genes differentially expressed only in the first experiment, genes differentially expressed only in the second experiment and genes differentially expressed in neither experiment. There are not true positive genes (i.e. truly DE in both experiments), so we should find no genes in common using our method.

Then, as described in Hwang et al., a two tails T-test is performed for each FC(gk) and a p-value is generated as: P(gk) = 2 Normal cdf(-absolute value (FC(gk)/r(k))) where FC(gk) is the t statistic that evaluates the differential expression between the two classes for the g gene and k experiment. 
}
\value{
  \item{names }{Which group each simulated gene expression value belongs to}
  \item{FC1}{T statistic for the first experiment}
  \item{FC2}{T statistic for the second experiment}
  \item{Pval}{p-values for the experiment to be compared}
}
\references{ Hwang D, Rust A, Ramsey S, Smith J, Leslie D, Weston A, de Atauri P, Aitchison J, Hood L, Siegel A, Bolouri H (2005): A data integration methodology for system biology. PNAS 2005.

M.Blangiardo and S.Richardson (2007) Statistical tools for synthesizing lists of differentially expressed features in related experiments, Genome Biology, 8, R54.  }
\author{ Alberto Cassese, Marta Blangiardo }
\examples{
data.indep = simulation.indep(n=500,GammaA=1,
GammaB=1,r1=0.5,r2=0.8,DEfirst=300,DEsecond=200)

}
