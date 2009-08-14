\name{simulation}
\Rdversion{1.1}
\alias{simulation}
\title{ Simulate p-values for two related experiments}
\description{
The function simulates two vectors of p-values using the procedure described in Hwang et al.
}
\usage{
simulation(n, GammaA, GammaB, epsilonM = 0, epsilonSD = 1, r1, r2, DEfirst, DEsecond, DEcommon)
}
\arguments{
  \item{n}{ \code{Number of features to simulate} }
  \item{GammaA}{ \code{Parameter of the Gamma distribution} }
  \item{GammaB}{ \code{Parameter of the Gamma distribution} }
  \item{epsilonM}{ \code{Parameter of the Gaussian noise specific to the genes and experiment} }
  \item{epsilonSD}{ \code{Parameter of the Gaussian noise specific to the genes and experiment} }
  \item{r1}{ \code{Additional experiment-specific noise} }
  \item{r2}{ \code{Additional experiment-specific noise} }
  \item{DEfirst}{ \code{Number of DE features in each experiment} }
  \item{DEsecond}{ \code{Number of DE features in each experiment} }
  \item{DEcommon}{ \code{Number of DE features in common between the two experiments} }
}
\details{
Considering two experiments (k=1,2), each of them with two classes, and n genes, for each gene we simulate a true difference between the classes delta(g), drawn from a Gamma distribution with random sign. The true difference delta(g) is 0 if the gene is not differentially expressed. We then add two normal random noise components, r(k) that act as experiment specific components and epsilon(gk), that are the gene-experiment components. The former is assigned deterministically, whilst the latter is drawn from a standard Gaussian distribution. The log fold change (FC(gk)) is the sum of all these components for each gene and experiment. We assign the n genes to four groups: genes differentially expressed (DE) in both experiments, genes differentially expressed only in the first experiment, genes differentially expressed only in the second experiment and genes differentially expressed in neither experiment. When the genes are differentially expressed in both experiments, they share the same delta(g) and the only difference between them is given by the random components: FC(g1) = delta(g) + r(1) times epsilon(g1) FC(g2) = delta(g) + r(2) times epsilon(g2) This group represents the true positive genes (i.e. truly DE in both experiments) that we are interested in finding using our method. The two groups of genes differentially expressed only in one of the two experiments act like additional noise and make the simulation more realistic.

Then, as described in Hwang et al., a two tails T-test is performed for each FC(gk) and a p-value is generated as: P(gk) = 2 Normal cdf(-absolute value (FC(gk)/r(k))). 
}
\value{
  \item{names}{Which group each simulated gene expression value belongs to}
  \item{FC1}{T statistic for the first experiment}
  \item{FC2}{T statistic for the second experiment}
  \item{Pval}{p-value for the experiments to be compared}
}
\references{ Hwang D, Rust A, Ramsey S, Smith J, Leslie D, Weston A, de Atauri P, Aitchison J, Hood L, Siegel A, Bolouri H (2005): A data integration methodology for system biology. PNAS 2005.

M.Blangiardo and S.Richardson (2007) Statistical tools for synthesizing lists of differentially expressed features in related experiments, Genome Biology, 8, R54. }
\author{ Marta Blangiardo, Alberto Cassese }

\examples{
data = simulation(n=500,GammaA=1,GammaB=1,r1=0.5,r2=0.8,DEfirst=300,DEsecond=200,DEcommon=100)

## The function is currently defined as
function(n,GammaA,GammaB,epsilonM=0,epsilonSD=1,r1,r2,DEfirst,DEsecond,DEcommon){

FC1=c()
FC2=c()
delta=rgamma(n,GammaA,1/GammaB)
epsilon1=rnorm(n,epsilonM,epsilonSD)
epsilon2=rnorm(n,epsilonM,epsilonSD)
names=c()
#Group 1 : DE in common
for(i in 1: DEcommon){
x=rbinom(1,1,0.5)
if(x==1) {
FC1[i] <- delta[i] + epsilon1[i]*r1;
FC2[i] <- delta[i] + epsilon2[i]*r2
  }
if(x==0) {
FC1[i] <- -delta[i] - epsilon1[i]*r1;
FC2[i] <- -delta[i] - epsilon2[i]*r2
  }
names[i] <- "DEcommon"
  }

#Group 2 : DE in the first experiment
for(i in (DEcommon+1):(DEfirst)){
x=rbinom(1,1,0.5)
if(x==1){
FC1[i] <- delta[i] + epsilon1[i]*r1
  }
if(x==0){
FC1[i] <- -delta[i] - epsilon1[i]*r1
  }
FC2[i] <- epsilon2[i]*r2
names[i] <- "DEfirst"
  }

#Group 3 : DE in the second experiment
for(i in (DEfirst+1):(DEfirst+DEsecond-DEcommon)){
x=rbinom(1,1,0.5)
FC1[i] <- epsilon1[i]*r1
if(x==1){
FC2[i] <- delta[i] + epsilon2[i]*r2
  }
if(x==0){
FC2[i] <- -delta[i] - epsilon2[i]*r2
  }
names[i] <- "DEsecond"

  }

#Group 4 : Not DE in Both experiments
for(i in (DEfirst+DEsecond-DEcommon+1):(n)){
FC1[i] <- epsilon1[i]*r1
FC2[i] <- epsilon2[i]*r2
names[i] <- "Null"
  }

##############################################
#Assign the Pvalues

Pval1 = c()
Pval2 = c()

for(i in 1:n){
Pval1[i] <- 2*pnorm(-abs(FC1[i]/r1))
Pval2[i] <- 2*pnorm(-abs(FC2[i]/r2))
  }

##############################################
return(list(names=names,FC1=FC1,FC2=FC2,Pval=cbind(Pval1,Pval2)))
  }
}
