createTable <-
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

