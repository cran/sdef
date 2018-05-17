\name{Liver.Muscle}
\alias{Liver.Muscle}
\docType{data}
\title{
Diabetes susceptibility in liver and skeletal muscle of mice.
}
\description{
This dataset contains two lists of p-values obtained from a publicly available experiment to evaluate differential expression between diabetes susceptibility in liver and skeletal muscle of obese and normal mice.
}
\usage{data(Liver.Muscle)}
\format{
The format is: a matrix with 2912 rows and two columns. For each gene (row) it reports the p-values of being differentially expressed between obese and normal mice for the two tissues (columns).
}

\source{
\url{http://www.ncbi.nlm.nih.gov/geo}, accession number GDS1443
}
\references{
Lan H, Rabaglia ME, Stoehr JP, Nadler ST et al. Gene expression profiles of non diabetic and diabetic obese mice suggest a role of hepatic lipogenic capacity in diabetes susceptibility. Diabetes 2003 Mar;52(3):688-700.
}
\examples{
#data(Liver.Muscle)
#Th<- ratio(data=Liver.Muscle)

#Rh<- baymod(iter=100,output.ratio=Th)

#MC<- Tmc(iter=100,output.ratio=Th)

#The gene names can be obtained using the command dimnames:
#feat.names = dimnames(Liver.Muscle)[[1]]
#feat.lists <- extractFeatures.R(output.ratio=Th,output.bay=Rh,feat.names=feat.names,h=NULL)
#feat.lists.T <- extractFeatures.T(output.ratio=Th,feat.names=feat.names)

#output.table <- createTable(output.ratio=Th,output.bay=Rh)

}
