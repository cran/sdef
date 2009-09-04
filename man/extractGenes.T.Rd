\name{extractGenes.T}
\Rdversion{1.1}
\alias{extractGenes.T}
\title{Extracting the lists of genes of interest}
\description{
The function returns the list of genes in common using the two suggested rules qmax and q2 (Frequentist model).
}
\usage{
extractGenes.T(output.ratio, gene.names)
}
\arguments{
  \item{output.ratio}{ \code{The output object from the Frequentist model (ratio function)} }
  \item{gene.names}{ \code{ID of the genes (e.g Affy ID)} }
}
\details{
To select a list of interesting features from the frequentist model we suggest two decision rules in the paper: 1) the maximum of Median(T(q)); 2) the largest threshold q for which the ratio T(q) il bigger than 2.

The first one is pointing out the strongest deviation from independence, whilst the second is the largest threshold where the number of genes called in common at least doubles the number of genes in common under independence.
}
\value{
The function returns an object of the class list. Each element is a matrix where the first column contains the name of the genes while the other columns contain the p-values from the experiments: 
  \item{max }{The list of genes in common selected on the basis of the threshold associated to Tmax}
  \item{rule2 }{The list of genes in common selected on the basis of the threshold associated to T2}
}
\references{ 1. M.Blangiardo and S.Richardson (2007) Statistical tools for synthesizing lists of differentially expressed features in related experiments, Genome Biology, 8, R54 }
\author{ Alberto Cassese, Marta Blangiardo  }
\examples{
data = simulation(n=500,GammaA=1,GammaB=1,r1=0.5,r2=0.8,
DEfirst=300,DEsecond=200,DEcommon=100)
Tq<- ratio(data=data$Pval)
gene.names = data$names
gene.lists.T <- extractGenes.T(output.ratio=Tq,
gene.names=gene.names)

## The function is currently defined as
function(output.ratio,gene.names){
load(paste(output.ratio$dataname,".Rdata"))
if(output.ratio$pvalue==FALSE){
data = 1 - data
  }

dim1=dim(data)[1]
lists = dim(data)[2]

#Decision rules:
max.T = max(output.ratio$ratios)
threshold.max = 
output.ratio$q[output.ratio$ratios==max.T]

#function Table
table=function(threshold){
temp=matrix(0,dim1,lists)
for(i in 1:dim1){
for(j in 1:lists){
if(data[i,j]<= threshold){temp[i,j]<-1}
}
}
return(temp)
}

#Table
temp<-table(threshold.max)

name="Names"
for(i in 1:ncol(data)){
name=c(name,paste("List",as.character(i)))
  }

table.max <- data[apply(temp,1,sum)==lists,]
names.max <- gene.names[apply(temp,1,sum)==lists]

if(output.ratio$pvalue==FALSE){
table.max <- 1-table.max
}

if(is.matrix(table.max)==FALSE){
table.max=as.matrix(table.max)
table.max=t(table.max)
}

table.max <- data.frame(Names=names.max,
RankingStat = table.max)
colnames(table.max)<-name

if(length(output.ratio$q
[output.ratio$ratios>=2])>0){

#2) Rule 2
threshold.2 = max(output.ratio$q
[output.ratio$ratios>=2])

#Table
temp<-table(threshold.2)

table.2 <- data[apply(temp,1,sum)==lists,]
names.2 <- gene.names[apply(temp,1,sum)==lists]

if(output.ratio$pvalue==FALSE){
table.2 <- 1-table.2
}

if(is.matrix(table.2)==FALSE){
table.2=as.matrix(table.2)
table.2=t(table.2)
}

table.2 <- data.frame(Names=names.2,
RankingStat = table.2)
colnames(table.2)<-name

return(list(max = table.max,rule2 = table.2))
}
if(length(output.ratio$q
[output.ratio$ratios>=2])==0){

return(list(max = table.max))
}
}
}
