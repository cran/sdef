\name{Example3Lists}
\Rdversion{1.5}
\alias{Example3Lists}
\docType{data}
\title{
Molecular Differences between Mammalian Sexes.
}
\description{
This dataset contains three lists of p-values obtained from a publicly available experiment to evaluate differential expression between mammalian sexes in three tissues (hypothalamus, kidney and liver).
}
\usage{data(Example3Lists)}
\format{
The format is: a matrix with 6477 rows and 3 columns. For each gene (row) it reports the p-values of being differentially expressed between male and female mice for the three tissues.
}
\source{
\url{http://www.ncbi.nlm.nih.gov/geo}, accession number GSE1147-GSE1148
}
\references{
Rinn J, Rozowsky J, Laurenzi I, Petersen P, Zou ZW K, Gerstein M, Snyder1 M: Major
Molecular Differences between Mammalian Sexes Are Involved in Drug Metabolism and
Renal Function.
Developmental Cell 2004, 6:791-800
}
\examples{
#data(Example3Lists)
#Th<- ratio(data=Example3Lists)

#Rh<- baymod(iter=100,output.ratio=Th)

#MC<- Tmc(iter=100,output.ratio=Th)

#The gene names can be obtained using the command dimnames:
#feat.names = dimnames(Example3Lists)[[1]]
#feat.lists <- extractFeatures.R(output.ratio=Th,output.bay=Rh,feat.names=feat.names,h=NULL)
#feat.lists.T <- extractFeatures.T(output.ratio=Th,feat.names=feat.names)

#output.table <- createTable(output.ratio=Th,output.bay=Rh)

}
