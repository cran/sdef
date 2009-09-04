\name{Tmc}
\Rdversion{1.1}
\alias{Tmc}
\title{Empirical null distribution of max(T)}
\description{
The function uses bootstrap for calculating the empirical distribution of max(T)=T(q*) under the null hypothesis of independence among the experiments. An empirical p-value is calculated to evaluate how the data are far from the hypothesis of independence.
}
\usage{
Tmc(iter = 1000, output.ratio)
}
\arguments{
  \item{iter}{ \code{Number of iteration to be performed} }
  \item{output.ratio}{ \code{The output object from the ratio function} }
}
\details{
This function uses bootstrap for calculating the empirical distribution of the maximum of T (i.e. T(q*)) under the null hypothesis of independence among the experiments.  While the p-values for the first list are fixed, the ones for the other lists are independently permutate B times. In this way, any relationship among the lists is destroyed. At each permutation b (b varies from 1 to B) a Tb(q) is calculated for each q and a maximum statistic Tb(q*) is returned; its distribution represents the null distribution of T(q*) under the condition of independence. The relative frequency of Tb(q*) larger than T(q*) identifies the p-value: it returns the proportion of Tb(q*) from permuted dataset greater than the observed one (so indicates where the observed T(q*) is located on the null distribution). 
}
\value{
Returns the empirical pvalue from testing T(q*).
}
\references{ Stone et al.(1988), Investigations of excess environmental risks around putative sources: statistical problems and a proposed test,Statistics in Medicine, 7, 649-660.

M.Blangiardo and S.Richardson (2007) Statistical tools for synthesizing lists of differentially expressed features in related experiments, Genome Biology, 8, R54. }
\author{ Alberto Cassese, Marta Blangiardo }
\examples{
data = simulation(n=500,GammaA=1,GammaB=1,
r1=0.5,r2=0.8,DEfirst=300,DEsecond=200,
DEcommon=100)
Tq<- ratio(data=data$Pval)
bootstrap<- Tmc(iter=100,output.ratio=Tq)

## The function is currently defined as
function(iter=1000,output.ratio){
load(paste(output.ratio$dataname,".Rdata"))
if(output.ratio$pvalue==FALSE){
data=1-data
}

dim1=dim(data)[1]
lists = ncol(data)
l=length(output.ratio$Common)
Tmax = max(output.ratio$ratios,na.rm=TRUE)
Tmax.null = rep(NA,iter)
ratios.null = matrix(NA,l,iter)
sample = matrix(NA,dim1,lists)
    
sample[,1] <- data[,1]
for(k in 1:iter){
int = c()
L=matrix(0,l,lists)
data1 = matrix(NA,dim1,lists)
data1[,1] <- data[,1]

for(j in 2:lists){
sample[,j] = sample(data[,j])
data1[,j] = sample[,j]
}
    
threshold = output.ratio$q
for(i in 1:l){
temp = data1<=threshold[i]
for(j in 1:lists){
L[i,j] <- sum(temp[,j])
temp[temp[,j]==FALSE,j]<-0
temp[temp[,j]==TRUE,j]<-1
}
int[i] <- sum(apply(temp,1,sum)==lists)
}


expected = apply(L,1,prod)/(dim1)^(lists-1)
observed = int
ratios = matrix(0,l,1)

for(i in 1:l){
ratios[i,1] <- observed[i]/expected[i]
}
ratios.null[,k] <- ratios
ratios <- ratios[threshold>0]
Tmax.null[k] = max(ratios)
if(k\%\%10==0)cat(k,"of",iter,"completed","\n")
}

ID=seq(1,iter)
p=length(ID[Tmax.null>=Tmax])
pvalue<- p/iter

main="("

for(i in 1:(length(lists)-1)){
main=paste(main,"1,")
  }
main=paste(main,"1)")

postscript(paste("Pvalue",output.ratio$name,
".ps"))
hist(Tmax.null,main=main,xlab="T",
ylab="",xaxt="n",cex.main=0.9,
xlim=c(min(Tmax.null),
max(c(Tmax,max(Tmax.null)))),yaxt="n",
cex.axis=0.9)
axis(1,at = seq(min(Tmax.null),
max(c(Tmax,max(Tmax.null))),
length.out=10),
labels = round(seq(min(Tmax.null),
max(c(Tmax,max(Tmax.null))),
length.out=10),2))
if(pvalue>0){
legend(x=max(c(Tmax,max(Tmax.null)))-
min(Tmax.null)/2,y=dim1/100,
legend=paste("P value =",pvalue),
bty="n",cex=0.9)
abline(v=Tmax,lty=2)
dev.off()
return(list(pvalue=pvalue))
}
if(pvalue==0){
legend(x=max(c(Tmax,max(Tmax.null)))
-min(Tmax.null)/2,y=dim1/100,
legend=paste("P value <",1/iter),bty="n",
cex=0.9)
abline(v=Tmax,lty=2)
dev.off()
return(noquote(paste("pvalue < ",1/iter)))
}
}
}
