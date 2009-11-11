extractFeatures.R <-
function(output.ratio,output.bay,feat.names,q=NULL){
load(paste(output.ratio$dataname,".Rdata"))
if(output.ratio$pvalue==FALSE){
data = 1 - data
}

dim1=dim(data)[1]
lists = dim(data)[2]


#Decision rules:
#1) Maximum for CI not including 1
if(length(output.bay[round(output.bay[,1],3)>1,2])==0){
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
names.q <- feat.names[apply(temp,1,sum)==lists]
if(output.ratio$pvalue==FALSE){
        table.q[[r]] <- 1-table.q[[r]]
        }

table.q[[r]] <- data.frame(Names=names.q,RankingStat = table.q[[r]])

names(table.q)[[r]] <- paste("q=",q[r]) 
}
return(User=table.q)
}
}
if(length(output.bay[round(output.bay[,1],3)>1,2])>0){
max.R = max(output.bay[round(output.bay[,1],3)>1,2])
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
names.max <- feat.names[apply(temp,1,sum)==lists]

if(output.ratio$pvalue==FALSE){
table.max <- 1-table.max
}

if(is.vector(table.max)==TRUE){
table.max=as.matrix(table.max)
table.max=t(table.max)
}


table.max <- data.frame(Names=names.max,RankingStat = table.max)
colnames(table.max)<-name

if(length(output.ratio$q[output.bay[round(output.bay[,1],3)>1,2]>=2])>0){

#2) Rule 2
threshold.2 = max(output.ratio$q[round(round(output.bay[,2],3),3)>=2 & round(output.bay[,1],3)>1])

#Table
temp<-table(threshold.2)

table.2 <- data[apply(temp,1,sum)==lists,]
names.2 <- feat.names[apply(temp,1,sum)==lists]

if(output.ratio$pvalue==FALSE){
table.2 <- 1-table.2
}


if(is.vector(table.2)==TRUE){
table.2=as.matrix(table.2)
table.2=t(table.2)
}


table.2 <- data.frame(Names=names.2,RankingStat = table.2)
colnames(table.2)<-name

if(is.null(q)){
write.csv(table.max,"featuresRmax.csv")
write.csv(table.2,"featuresR2.csv")
return(list(max = table.max,rule2 = table.2))
}

if(!is.null(q)){
l = length(q)
table.q = list()
        for(r in 1:l){

temp <- table(q[r])

table.q[[r]] <- data[apply(temp,1,sum)==lists,]
names.q <- feat.names[apply(temp,1,sum)==lists]
if(output.ratio$pvalue==FALSE){
        table.q[[r]] <- 1-table.q[[r]]
        }

table.q[[r]] <- data.frame(Names=names.q,RankingStat = table.q[[r]])

names(table.q)[[r]] <- paste("q=",q[r]) 
}
}
write.csv(table.max,"featuresRmax.csv")
write.csv(table.2,"featuresR2.csv")
for (i in 1:length(q)){
write.csv(table.q[[i]],paste("features",q[i],".csv"))
}
return(list(max = table.max,rule2 = table.2, User = table.q))
        }

if(length(output.ratio$q[output.bay[round(output.bay[,1],3)>1,2]>=2])==0){

if(is.null(q)){
write.csv(table.max,"featuresRmax.csv")
return(list(max = table.max))
}

if(!is.null(q)){
l = length(q)
table.q = list()
for(r in 1:l){

temp <- table(q[r])

table.q[[r]] <- data[apply(temp,1,sum)==lists,]
names.q <- feat.names[apply(temp,1,sum)==lists]
if(output.ratio$pvalue==FALSE){
table.q[[r]] <- data.frame(Names=names.q,RankingStat = 1-table.q[[r]])
}

if(output.ratio$pvalue==TRUE){
table.q[[r]] <- data.frame(Names=names.q,RankingStat = table.q[[r]])
}
names(table.q)[[r]] <- paste("q=",q[r]) 
}
}
write.csv(table.max,"featuresRmax.csv")
for (i in 1:length(q)){
write.csv(table.q[[i]],paste("features",q[i],".csv"))
}
return(list(max = table.max,User = table.q))
}

}
}

