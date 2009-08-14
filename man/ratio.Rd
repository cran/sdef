\name{ratio}
\Rdversion{1.1}
\alias{ratio}
\title{Ratio T between the observed genes in common and the expected ones}
\description{
The function for each experiment calculates the ratio T(q) for each threshold q, using the list of p-values.
}
\usage{
ratio(data, pvalue = TRUE, interval = 0.01, name = NULL, dir = getwd(), dataname = "dataratio")
}
\arguments{
  \item{data}{ \code{Lists of pvalues to be compared} }
  \item{pvalue}{ \code{Indicate if the data are pvalues (TRUE) or posterior probability (FALSE). If they are posterior probability they are transformed in pvalues} }
  \item{interval}{ \code{The interval between two threshold} }
  \item{name}{ \code{The name to be used in the plots} }
  \item{dir}{ \code{Directory for storing the plots} }
  \item{dataname}{ \code{The name of the file containing the data (Pvalue)} }
}
\details{
This function calculates the ratio T(q) of observed number of genes in common between  the lists vs the expected number under the hypothesis of independence for each threshold q. The expected numbers are calculated as the product among the marginals divided by $(numbers of genes)^(number of lists-1)$. T(qmax) identifies the maximum of the statistic T(q) and it is shown on the plot.
}
\value{
Returns a plot with the distribution of T(q) showing where T(qmax) and qmax are located. It returns also an object of class list with the ratio, the thresholds and other attributes. In particular: 
  \item{ratios }{Vector or T values for each threshold}
  \item{q }{Threshold corresponding to T values}
  \item{DE }{Differentially expressed features in each experiment}
  \item{Common }{Features in common corresponding to the T values}
  \item{interval }{Interval on the p-value scale defined by the user (default is 0.01)}
  \item{name }{Names to be used in the plots (defined by the user)}
  \item{pvalue }{Logical: TRUE if the measures used for the analysis are p-value, FALSE if they are posterior probabilities}
  \item{dataname }{The name of the file where the data has been saved}
}
\references{ Stone et al.(1988), Investigations of excess environmental risks around putative sources: statistical problems and a proposed test,Statistics in Medicine, 7, 649-660.

M.Blangiardo and S.Richardson (2007) Statistical tools for synthesizing lists of differentially expressed features in related experiments, Genome Biology, 8, R54. }
\author{ Marta Blangiardo, Alberto Cassese }

\examples{
data = simulation(n=500,GammaA=1,GammaB=1,r1=0.5,r2=0.8,DEfirst=300,DEsecond=200,DEcommon=100)
Tq<- ratio(data=data$Pval)

## The function is currently defined as
function(data,pvalue=TRUE,interval=0.01,name=NULL,dir=getwd(),dataname="dataratio"){

#Define how many lists for the comparison
lists = ncol(data)
dim1=dim(data)[1]

if(pvalue==FALSE){
data=1-data;
  }

main="("
for(i in 1:(length(lists)-1)){
main=paste(main,"1,")
  }
main=paste(main,"1)")


if(length(name)==0){
name=main
  }

ID=seq(1,dim1)
threshold = seq(0,1,interval)
l=length(threshold)
int=c()
L=matrix(0,l,lists)

for(i in 1:l){
temp = data<=threshold[i]
for(j in 1:lists){
L[i,j] <- sum(temp[,j])
temp[temp[,j]==FALSE,j]<-0
temp[temp[,j]==TRUE,j]<-1
        }
        int[i] <- sum(apply(temp,1,sum)==lists)
  }

#Calculate the ratio for the number of observed genes/number of expected ones
expected = apply(L,1,prod)/(dim1)^(lists-1)
ratios = matrix(0,l,1)

for(i in 1:l){
ratios[i,1] <- int[i]/expected[i]
  }

#Plot

ratios=ratios[int>0]
thresh.ratios=threshold[int>0]
L = L[int>0,] 
int=int[int>0]

Tmax = max(ratios)
qmax = thresh.ratios[ratios==Tmax]
if(length(qmax)>1){qmax=qmax[1]}
if(Tmax<1){cat("WARNING: the requested contrast is under-represented in the data (Tmax<1)\n")}
ps.options(paper="a4",horizontal=TRUE)
setwd(dir)
ps.options(horizontal=FALSE)
postscript(paste("Ratio",name,".ps"))

if(length(thresh.ratios[ratios>=2])>0){
q2 = max(thresh.ratios[ratios>=2])
T2 = ratios[thresh.ratios==q2]
if(q2==qmax){
plot(thresh.ratios,ratios,type="l",
ylab= "T",xlab="P-value",main=main,yaxt="n",xaxt="n",cex.main=0.7,cex.axis=1.2,ylim=c(0,(max(ratios,na.rm=TRUE)+sd(ratios,na.rm=TRUE))))
if(Tmax<1.1){
axis(2, at = c(0,0.5,Tmax), labels = c(0,0.5,expression(T[max])), tick = TRUE,cex.axis=0.9)
  }
if(Tmax>1.1){
axis(2, at = c(0,0.5,1,Tmax), labels = c(0,0.5,1,expression(T[max])), tick = TRUE,cex.axis=0.9)
  }
if(qmax>0.1 & qmax<0.9){
axis(1, at = c(0,qmax,1), labels = c(0,expression(q[max]),1), tick=TRUE,cex=0.9)
  }
if(qmax<0.1){
axis(1, at = c(qmax,1), labels = c(expression(q[max]),1), tick=TRUE,cex=0.9)
  }
if(qmax>0.9 & qmax != 1){
axis(1, at = c(0,qmax), labels = c(0,expression(q[max])), tick=TRUE,cex=0.9)
  }
if(qmax==1){
axis(1, at = c(0,1), labels = c(0,paste(expression(q[max]),"=1")), tick=TRUE,cex=0.9)
  }
axis(4, at = c(1,Tmax),labels = c(dim1,int[thresh.ratios==qmax]),tick=TRUE,cex=0.9)
dev.off()
  }
if(q2!=qmax){
plot(thresh.ratios,ratios,type="l",
ylab= "T",xlab="P-value",main=main,yaxt="n",xaxt="n",cex.main=0.7,cex.axis=1.2,ylim=c(0,(max(ratios,na.rm=TRUE)+sd(ratios,na.rm=TRUE))))
if(Tmax<1.1){
axis(2, at = c(0,0.5,Tmax,T2), labels = c(0,0.5,expression(T[max]),expression(T[2])), tick = TRUE,cex.axis=0.9)
  }
if(Tmax>1.1){
axis(2, at = c(0,0.5,1,Tmax,T2), labels = c(0,0.5,1,expression(T[max]),expression(T[2])), tick = TRUE,cex.axis=0.9)
  }

if(qmax>0.1 & qmax<0.9){
axis(1, at = c(0,qmax,q2,1), labels = c(0,expression(q[max]),expression(q[2]),1), tick=TRUE,cex=0.9)
  }
if(qmax<0.1){
axis(1, at = c(qmax,q2,1), labels = c(expression(q[max]),expression(q[2]),1), tick=TRUE,cex=0.9)
  }
if(qmax>0.9 & qmax != 1){
axis(1, at = c(0,qmax,q2), labels = c(0,expression(q[max]),expression(q[2])), tick=TRUE,cex=0.9)
  }
if(qmax==1){
axis(1, at = c(0,q2,1), labels = c(0,expression(q[2]),paste(expression(q[max]),"=1")), tick=TRUE,cex=0.9)
  }

axis(4, at = c(1,T2,Tmax),labels = c(dim1,int[thresh.ratios==q2],int[thresh.ratios==qmax]),tick=TRUE,cex=0.9)
dev.off()
  }
  }
    
if(length(thresh.ratios[ratios>=2])==0){
plot(thresh.ratios,ratios,type="l",
ylab= "T",xlab="P-value",main=main,yaxt="n",xaxt="n",cex.main=0.7,cex.axis=1.2,ylim=c(0,(max(ratios,na.rm=TRUE)+sd(ratios,na.rm=TRUE))))

if(Tmax<1.1){
axis(2, at = c(0,0.5,Tmax), labels = c(0,0.5,expression(T[max])), tick = TRUE,cex.axis=0.9)
  }
if(Tmax>1.1){
axis(2, at = c(0,0.5,1,Tmax), labels = c(0,0.5,1,expression(T[max])), tick = TRUE,cex.axis=0.9)
  }
if(qmax>0.1 & qmax<0.9){
axis(1, at = c(0,qmax,1), labels = c(0,expression(q[max]),1), tick=TRUE,cex=0.9)
  }
if(qmax<0.1){
axis(1, at = c(qmax,1), labels = c(expression(q[max]),1), tick=TRUE,cex=0.9)
  }
if(qmax>0.9 & qmax != 1){
axis(1, at = c(0,qmax), labels = c(0,expression(q[max])), tick=TRUE,cex=0.9)
  }
if(qmax==1){
axis(1, at = c(0,1), labels = c(0,paste(expression(q[max]),"=1")), tick=TRUE,cex=0.9)
  }

axis(4, at = c(1,Tmax),labels = c(dim1,int[thresh.ratios==qmax]),tick=TRUE,cex=0.9)
dev.off() 
  }
save(data,file = paste(dataname,".Rdata"))
return(list(DE = L, ratios=ratios,q=thresh.ratios,Common=int,interval=interval,name=name,pvalue=pvalue,dataname=dataname))

  }
}
