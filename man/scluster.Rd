% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scluster.R
\docType{package}
\name{scluster}
\alias{scluster}
\alias{scluster-package}
\title{Supervised Sparse k-means Clustering}
\description{
Supervised (and unsupervised) sparse k-means clustering, with 
single-cell RNA-seq data in mind.  The clustering can be performed using 
either a traditional (squared) euclidean distance, or via a dissimilarity 
matrix that takes into account the drop-out phenomenom in single-cell  
RNA-seq data. The dissimilarity measure is obtained through the 
single-cell analysis package CIDR. License: GPL (>=3)
}
\examples{
par(ask=FALSE)
## Generate simulated single-cell RNA-Seq tags.
N=3 ## 3 cell types
k=50 ## 50 cells per cell type
sData <- scSimulator(N=N, k=k)
## tags - the tag matrix
tags <- as.matrix(sData$tags)
cols <- c(rep("RED",k), rep("BLUE",k), rep("GREEN",k))
## Standard principal component analysis.
ltpm <- log2(t(t(tags)/colSums(tags))*1000000+1)
pca <- prcomp(t(ltpm))
plot(pca$x[,c(1,2)],col=cols,pch=1,xlab="PC1",ylab="PC2",main="prcomp")
## use scluster to cluster the data - unsupervised, euclidean distance
## The input for scluster should be a tag matrix.
nCluster <- 3
scSim <- sclDataConstructor(tags)
scSim <- sclust(scSim, nCluster=nCluster, nf=200, 
                max.nIter=10, markers=NA, supervise=FALSE, alpha=0)
## Use Adjusted Rand Index to measure the accuracy of the clustering output by cidr.
scluster::adjustedRandIndex(scSim@cluster,cols)
## 0.98

## Now redo the clustering using the CIDR "distance"
scSim <- sclDataConstructor(tags)
scSim <- sclust(scSim, nCluster=nCluster, nf=200, dist="cidr",
                max.nIter=10, markers=NA, supervise=FALSE, alpha=0)
## Use Adjusted Rand Index to measure the accuracy of the clustering output by cidr.
scluster::adjustedRandIndex(scSim@cluster,cols)
## 1.0

## Now perform clustering using the parameter "s" instead of "nf".
## NOTE that we use either but not both of the parameters "s" or "nf".
## s - A tuning parameter determines the sum of weights for all features.
## nf - Parameter that determines the number of features (with non-zero weight) 
## used for clustering.
scSim <- sclDataConstructor(tags)
scSim <- sclust(scSim, nCluster=nCluster, s=10, 
                max.nIter=10, markers=NA, supervise=FALSE, alpha=0)
## Use Adjusted Rand Index to measure the accuracy of the clustering 
##output by cidr.
scluster::adjustedRandIndex(scSim@cluster,cols)
## 1.0

## Input tags are cpm. If using cidr distance, it is preferable to use raw 
## tags, unless they are not available.
tags_cpm <- t(t(tags)/colSums(tags))*1000000
scSim <- sclDataConstructor(tags, tagType = "cpm")
scSim <- sclust(scSim, nCluster=nCluster,max.nIter=10, 
                markers=NA, supervise=FALSE, nf=300, alpha=0)
scluster::adjustedRandIndex(scSim@cluster,cols)
## 0.98
}
\author{
Tingting Gong, Michael Troup
}

