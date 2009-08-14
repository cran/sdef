ratio <-
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

