extractGenes.T <-
function(output.ratio,gene.names){
load(paste(output.ratio$dataname,".Rdata"))
if(output.ratio$pvalue==FALSE){
data = 1 - data
}

dim1=dim(data)[1]
lists = dim(data)[2]


#Decision rules:
#1) Maximum for CI not including 1
max.T = max(output.ratio$ratios)
threshold.max = output.ratio$q[max(output.ratio$ratios)]

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

if(length(output.ratio$q[output.ratio$ratios>=2])>0){

#2) Rule 2
threshold.2 = max(output.ratio$q[output.ratio$ratios>=2])

#Table
temp<-table(threshold.2)

table.2 <- data[apply(temp,1,sum)==lists,]
names.2 <- gene.names[apply(temp,1,sum)==lists]

if(output.ratio$pvalue==FALSE){
table.2 <- 1-table.2
}


table.2 <- data.frame(Names=names.2,RankingStat = table.2)
colnames(table.2)<-name

return(list(max = table.max,rule2 = table.2))
}
if(length(output.ratio$q[output.ratio$ratios>=2])==0){

return(list(max = table.max))
}
}

