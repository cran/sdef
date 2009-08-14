\name{extractGenes.R}
\Rdversion{1.1}
\alias{extractGenes.R}
\title{Extracting the lists of genes of interest}
\description{
The function returns the list of genes in common using the two suggested rules qmax and q2 (Bayesian model) and additional ones defined by the user.
}
\usage{
extractGenes.R(output.ratio, output.bay, gene.names, q = NULL)
}
\arguments{
  \item{output.ratio}{ \code{The output object from the Frequentist model (ratio function)} }
  \item{output.bay}{ \code{The output object from the Bayesian model (baymod function)} }
  \item{gene.names}{ \code{ID of the genes (e.g Affy ID)} }
  \item{q}{ \code{Additional thresholds in the form of a vector to select a list of genes in common. If it is NULL only the Rmax and rule2 are used to select the lists of genes of interest} }
}
\details{
To select a list of interesting features from the Bayesian model we suggest two decision rules in the paper: 1) the maximum of Median(R(q)) only for the subset of credibility intervals which do not include 1; 2) the largest threshold q for which the ratio R(q) il bigger than 2.

The first one is pointing out the strongest deviation from independence, whilst the second is the largest threshold where the number of genes called in common at least doubles the number of genes in common under independence. The user can define an additional threshold q of interest and obtain the list of genes associated with it.
}
\value{
The function returns an object of the class list. Each element is a matrix where the first column contains the name of the genes while the other columns contain the p-values from the experiments: 
  \item{max }{The list of genes of interest selected on the basis of the threshold associated to Rmax}
  \item{rule2 }{The list of genes of interest selected on the basis of the threshold associated to R2}
  \item{table.q}{The list of genes of interest selected on the basis of the additional thresholds selected by the user}
}
\references{ 1. M.Blangiardo and S.Richardson (2007) Statistical tools for synthesizing lists of differentially expressed features in related experiments, Genome Biology, 8, R54 }
\author{ Marta Blangiardo, Alberto Cassese  }
\examples{
data = simulation(n=500,GammaA=1,GammaB=1,r1=0.5,r2=0.8,DEfirst=300,DEsecond=200,DEcommon=100)
Tq<- ratio(data=data$Pval)
Rq<- baymod(iter=100,output.ratio=Tq)
gene.names = data$names
gene.lists <- extractGenes.R(output.ratio=Tq,output.bay=Rq,gene.names=gene.names,q=NULL)

## The function is currently defined as
function(output.ratio,output.bay,gene.names,q=NULL){
load(paste(output.ratio$dataname,".Rdata"))
if(output.ratio$pvalue==FALSE){
data = 1 - data
  }

dim1=dim(data)[1]
lists = dim(data)[2]


#Decision rules:
#1) Maximum for CI not including 1
if(length(output.bay[round(output.bay[,1],2)>1,2])==0){
cat("WARNING: the requested contrast is under-represented in the data (Rmax<1)\n")
if(!is.null(q)){
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


l = length(q)
table.q = list()
        for(r in 1:l){
temp <- table(q[r])
table.q[[r]] <- data[apply(temp,1,sum)==lists,]
names.q <- gene.names[apply(temp,1,sum)==lists]
if(output.ratio$pvalue==FALSE){
        table.q[[r]] <- 1-table.q[[r]]
        }

table.q[[r]] <- data.frame(Names=names.q,RankingStat = table.q[[r]])

names(table.q)[[r]] <- paste("q=",q[r]) 
  }
return(User=table.q)
  }
  }
if(length(output.bay[round(output.bay[,1],2)>1,2])>0){
max.R = max(output.bay[round(output.bay[,1],2)>1,2])
threshold.max = output.ratio$q[output.bay[,2]==max.R]

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

table.max <- data.frame(Names=names.max,RankingStat = table.max)
colnames(table.max)<-name

if(length(output.ratio$q[output.bay[round(output.bay[,1],2)>1,2]>=2])>0){

#2) Rule 2
threshold.2 = max(output.ratio$q[round(round(output.bay[,2],2),3)>=2 & round(output.bay[,1],2)>1])

#Table
temp<-table(threshold.2)

table.2 <- data[apply(temp,1,sum)==lists,]
names.2 <- gene.names[apply(temp,1,sum)==lists]

if(output.ratio$pvalue==FALSE){
table.2 <- 1-table.2
  }


table.2 <- data.frame(Names=names.2,RankingStat = table.2)
colnames(table.2)<-name

if(is.null(q)){return(list(max = table.max,rule2 = table.2))}

if(!is.null(q)){
l = length(q)
table.q = list()
        for(r in 1:l){

temp <- table(q[r])

table.q[[r]] <- data[apply(temp,1,sum)==lists,]
names.q <- gene.names[apply(temp,1,sum)==lists]
if(output.ratio$pvalue==FALSE){
        table.q[[r]] <- 1-table.q[[r]]
        }

table.q[[r]] <- data.frame(Names=names.q,RankingStat = table.q[[r]])

names(table.q)[[r]] <- paste("q=",q[r]) 
  }
  }
return(list(max = table.max,rule2 = table.2, User = table.q))
        }

if(length(output.ratio$q[output.bay[round(output.bay[,1],2)>1,2]>=2])==0){

if(is.null(q)){return(list(max = table.max))}

if(!is.null(q)){
l = length(q)
table.q = list()
for(r in 1:l){

temp <- table(q[r])

table.q[[r]] <- data[apply(temp,1,sum)==lists,]
names.q <- gene.names[apply(temp,1,sum)==lists]
if(output.ratio$pvalue==FALSE){
table.q[[r]] <- data.frame(Names=names.q,RankingStat = 1-table.q[[r]])
  }

if(output.ratio$pvalue==TRUE){
table.q[[r]] <- data.frame(Names=names.q,RankingStat = table.q[[r]])
  }
names(table.q)[[r]] <- paste("q=",q[r]) 
  }
  }
return(list(max = table.max,User = table.q))
  }

  }
  }
}
