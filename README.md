fNet
====

an R program, short of Filter by highly ranked gene for Support Vector Machine, for microarray classification. 
R program and Supplment meterial for "Yupeng Cun, Holger Fröhlich (2012) Integrating Prior Knowledge Into Prognostic Biomarker Discovery Based on Network Structure. arXiv:1212.3214 "

FrSVM is the short of Filter by highly ranked gene for Support Vector Mchine. Download source codes here.  FrSVM is a feature selection algorithm which integrates protein-protein interaction network information into gene selection for prognostic biomarker discovery. 

How to run FrSVM:

1. Geting gene expression profiles (GEP), PPi Network.

##############################################
# Geing GEP
#----------------------------------------------------------------------------------
library(GEOquery) 
a = getGEO("GSExxxxx", destdir="/home/YOURPATH/")
## Normalized the GEP by limma
x= t(normalizeBetweenArrays(exprs(a), method="quantile") )
## defien your classes labes, y, as a factor
y= facotr("Two Class")

##############################################
# mapping probest IDs to Entrez IDs
# take hgu133a paltform as example
#---------------------------------------------------------------------------------
library('hgu133a.db')
mapped.probes<-mappedkeys(hgu133aENTREZID)
refseq<-as.list(hgu133aENTREZID[mapped.probes])
times<-sapply(refseq, length)
mapping <- data.frame(probesetID=rep(names(refseq),times=times), graphID=unlist(refseq),row.names=NULL, stringsAsFactors=FALSE)
mapping<- unique(mapping)

##############################################
# Summarize probests to genes of x by limma
# ad.ppi: Adjacencen matrix of PPI network
#---------------------------------------------------------------------------------
Gsub=ad.ppi
mapping <- mapping[mapping[,'probesetID'] %in% colnames(x),]
int <- intersect(rownames(Gsub), mapping[,"graphID"])
xn.m=xn.m[,mapping$probesetID]
   
index = intersect(mapping[,'probesetID'],colnames(xn.m))    
x <- x[,index]
colnames(xn.m) <- map2entrez[index]
ex.sum = t(avereps(t(xn.m), ID=map2entrez[index]))   

int= intersect(int, colnames(ex.sum))
ex.sum=ex.sum[,int]         ## GEP which matched to PPI network
Gsub=Gsub[int,int]            ## PPI network which matched to GEP


2.  Run FrSVM program
##################################################
# You need install for flowing packages for run FrSVM.R programs:
#    library(ROCR)
#    library(Matrix)
#    library(kernlab)
#
## If you want to running parallelly, you also need  to load:
#    library(multicore)
#  
## Here is an expale for 5 times 10-folds Cross-Validtaion
source("../FrSVM.R")
res <- frSVM.cv(x=ex.sum, y=y, folds=10,Gsub=Gsub, repeats=5, parallel = FALSE, cores = 2, DEBUG=TRUE,d=0.5,top.uper=0.95,top.lower=0.9)
## the AUC values for 5*10-folds CV
AUC= res$auc
