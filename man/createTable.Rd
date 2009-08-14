\name{createTable}
\Rdversion{1.1}
\alias{createTable}
\title{Function to create an output table}
\description{
This function creates an output table with the results from the Bayesian model in a csv format.
}
\usage{
createTable(output.ratio, output.bay, dir = getwd())
}
\arguments{
  \item{output.ratio}{ \code{The output object from the Frequentist model (ratio function)} }
  \item{output.bay}{ 	\code{The output object from the Bayesian model (baymod function)} }
  \item{dir}{ \code{Directory for storing the table} }
}
\details{
To select a list of interesting features from the Bayesian model we suggest two decision rules in the paper: 1) the maximum of Median(R(q)) only for the subset of credibility intervals which do not include 1; 2) the largest threshold q for which the ratio R(q) il bigger than 2.

The first one is pointing out the strongest deviation from independence, whilst the second is the largest threshold where the number of genes called in common at least doubles the number of genes in common under independence. 
}
\value{
  \item{max }{The results of the Rmax statistic}
  \item{rule2 }{The results using the rule R larger than 2 (see details)}
}
\references{ 1. M.Blangiardo and S.Richardson (2007) Statistical tools for synthesizing lists of differentially expressed features in related experiments , Genome Biology, 8, R54  }
\author{ Marta Blangiardo, Alberto Cassese }
\examples{
data = simulation(n=500,GammaA=1,GammaB=1,r1=0.5,r2=0.8,DEfirst=300,DEsecond=200,DEcommon=100)
Tq<- ratio(data=data$Pval)
Rq<- baymod(iter=100,output.ratio=Tq)
output.table <- createTable(output.ratio=Tq,output.bay=Rq)

## The function is currently defined as
function(output.ratio,output.bay,dir=getwd()){

if(output.ratio$pvalue==TRUE){
matrix.results =  cbind(output.ratio$q,output.ratio$ratio,round(output.bay,3),output.ratio$Common,output.ratio$DE)
}
if(output.ratio$pvalue==FALSE){
matrix.results =  cbind(1-output.ratio$q,output.ratio$ratio,round(output.bay,3),output.ratio$Common,output.ratio$DE)
}

lists = dim(output.ratio$DE)[2]
namesDE = paste("O",seq(1,lists),rep("+",lists))
names.matrix = c("q","T",colnames(output.bay)[1],colnames(output.bay)[2],colnames(output.bay)[3],"O11",namesDE)
dimnames(matrix.results)[[2]]<-names.matrix

setwd(dir)
write.csv(matrix.results,paste("Output",output.ratio$name,".csv"),row.names=FALSE)

#Decision rules:
#1) Maximum for CI not including 1
if(length(output.bay[round(output.bay[,1],2)>1,2])==0){
cat("WARNING: the requested contrast is under-represented in the data (Rmax<1)\n")
}
if(length(output.bay[round(output.bay[,1],2)>1,2])>0){
max.R = max(output.bay[round(output.bay[,1],2)>1,2])
maximum1 = matrix.results[matrix.results[,4]==round(max.R,3),]
maximum=c(maximum1[1],maximum1[2],maximum1[3],maximum1[4],maximum1[5],maximum1[6])
for(i in 1:lists){
maximum=c(maximum,maximum1[6+i])
}
maximum=matrix(maximum,nrow=1,ncol=6+lists)
dimnames(maximum)[[2]]<-names.matrix

if(length(output.bay[output.bay[round(output.bay[,1],2)>1,2]>=2,1])>0){
#2) Rule 2
R2 = max(matrix.results[round(output.bay[,2],2)>=2 & round(output.bay[,1],2)>1 ,1])
rule2_temp = matrix.results[matrix.results[,1]==R2,]
rule2=c(rule2_temp[1],rule2_temp[2],rule2_temp[3],rule2_temp[4],rule2_temp[5],rule2_temp[6])
for(i in 1:lists){
rule2=c(rule2,rule2_temp[6+i])
}
rule2=matrix(rule2,nrow=1,ncol=6+lists)
dimnames(rule2)[[2]]<-names.matrix
return(list(maximum=maximum,rule2=rule2))
}

if(length(output.bay[output.bay[round(output.bay[,1],2)>1,2]>=2,1])==0){
#2) Rule 2
return(list(maximum=maximum))
}

}
}
}
