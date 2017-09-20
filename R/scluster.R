#' @title Supervised Sparse k-means Clustering
#'
#' @description Supervised (and unsupervised) sparse k-means clustering, with 
#' single-cell RNA-seq data in mind.  The clustering can be performed using 
#' either a traditional (squared) euclidean distance, or via a dissimilarity 
#' matrix that takes into account the drop-out phenomenom in single-cell  
#' RNA-seq data. The dissimilarity measure is obtained through the 
#' single-cell analysis package CIDR. License: GPL (>=3)
#'
#' @author Tingting Gong, Michael Troup
#'
#' @docType package
#' @name scluster
#' @useDynLib scluster
#'
#' @examples
#' par(ask=FALSE)
#' ## Generate simulated single-cell RNA-Seq tags.
#' N=3 ## 3 cell types
#' k=50 ## 50 cells per cell type
#' sData <- scSimulator(N=N, k=k)
#' ## tags - the tag matrix
#' tags <- as.matrix(sData$tags)
#' cols <- c(rep("RED",k), rep("BLUE",k), rep("GREEN",k))
#' ## Standard principal component analysis.
#' ltpm <- log2(t(t(tags)/colSums(tags))*1000000+1)
#' pca <- prcomp(t(ltpm))
#' plot(pca$x[,c(1,2)],col=cols,pch=1,xlab="PC1",ylab="PC2",main="prcomp")
#' ## use scluster to cluster the data - unsupervised, euclidean distance
#' ## The input for scluster should be a tag matrix.
#' nCluster <- 3
#' scSim <- sclDataConstructor(tags)
#' scSim <- sclust(scSim, nCluster=nCluster, nf=200, 
#'                 max.nIter=10, markers=NA, supervise=FALSE, alpha=0)
#' ## Use Adjusted Rand Index to measure the accuracy of the clustering output by cidr.
#' scluster::adjustedRandIndex(scSim@cluster,cols)
#' ## 0.98
#'
#' ## Now redo the clustering using the CIDR "distance"
#' scSim <- sclDataConstructor(tags)
#' scSim <- sclust(scSim, nCluster=nCluster, nf=200, dist="cidr",
#'                 max.nIter=10, markers=NA, supervise=FALSE, alpha=0)
#' ## Use Adjusted Rand Index to measure the accuracy of the clustering output by cidr.
#' scluster::adjustedRandIndex(scSim@cluster,cols)
#' ## 1.0
#'
#' ## Now perform clustering using the parameter "s" instead of "nf".
#' ## NOTE that we use either but not both of the parameters "s" or "nf".
#' ## s - A tuning parameter determines the sum of weights for all features.
#' ## nf - Parameter that determines the number of features (with non-zero weight) 
#' ## used for clustering.
#' scSim <- sclDataConstructor(tags)
#' scSim <- sclust(scSim, nCluster=nCluster, s=10, 
#'                 max.nIter=10, markers=NA, supervise=FALSE, alpha=0)
#' ## Use Adjusted Rand Index to measure the accuracy of the clustering 
#' ##output by cidr.
#' scluster::adjustedRandIndex(scSim@cluster,cols)
#' ## 1.0
#'
#' ## Input tags are cpm. If using cidr distance, it is preferable to use raw 
#' ## tags, unless they are not available.
#' tags_cpm <- t(t(tags)/colSums(tags))*1000000
#' scSim <- sclDataConstructor(tags, tagType = "cpm")
#' scSim <- sclust(scSim, nCluster=nCluster,max.nIter=10, 
#'                 markers=NA, supervise=FALSE, nf=300, alpha=0)
#' scluster::adjustedRandIndex(scSim@cluster,cols)
#' ## 0.98
NULL
#' @title Adjusted Rand Index
#' @author Chris Fraley, Adrian Raftery, Luca Scrucca.
#' @rdname adjustedRandIndex
#' @name adjustedRandIndex
#' @description Calculates the Adjusted Rand Index which meansures the accuracy of clustering when the ground truth is known.
#' @importFrom mclust adjustedRandIndex
#' @export adjustedRandIndex
#' @details Imported from the package \emph{mclust}; see \code{adjustedRandIndex::mclust} help page for more details.
#' @references
#' Chris Fraley, Adrian E. Raftery, T. Brendan Murphy, and Luca Scrucca (2012) mclust Version 4 for R: Normal Mixture Modeling for Model-Based Clustering, Classification, and Density Estimation Technical Report No. 597, Department of Statistics, University of Washington
#'
#' Chris Fraley and Adrian E. Raftery (2002) Model-based Clustering, Discriminant Analysis and Density Estimation Journal of the American Statistical Association 97:611-631
NULL
#' @title Single-cell RNA-seq data simulator
#' @rdname scSimulator
#' @name scSimulator
#' @export scSimulator
#' @importFrom cidr scSimulator
#' @details Single cell simulation function from the CIDR package.
NULL
## class sclData - single-cell RNA-Seq data object with attributes relevant to
## clustering through imputation and dimensionality reduction
setClass("sclData", representation(tags="matrix",
                                  tagType="character",
                                  tags.lcpm="matrix",
                                  sampleSize="numeric",
                                  nCluster="numeric",
                                  s="numeric",
                                  nf="numeric",
                                  max.nIter="numeric",
                                  markers="vector",
                                  supervise="logical",
                                  alpha="numeric",
                                  seed="numeric",
                                  dist="character",
                                  priorTPM="numeric",
                                  cluster="vector",
                                  weight="matrix",
                                  w.bcss="vector",
                                  bcss="matrix",
                                  totFeatDist="matrix",
                                  sumWeight="vector",
                                  nf.zw="vector"))
         
#' @title sclData Constructor
#'
#' @rdname sclDataConstructor
#'
#' @description
#' \code{sclDataConstructor} creates a new sclData class object from a tag table.
#'
#' @details
#' Creates an object in sclData (single-cell RNA-Seq dataset) class.
#' Attributes of the class include scalar, vector and matrix
#' data types necessary for the \emph{sclust} analysis - such as tag table, 
#' library  sizes, dropout candidates, imputation weighting threshold. The 
#' tags can be  raw tags (default) or tags per million (cpm).  Raw tags 
#' are preferrable as the individual library sizes, as determined by the raw 
#' tags, are used to determine dropout candidates.  
#'
#' @param tags A matrix of tags where the rows correspond to features (genes, transcripts, etc) and the columns correspond to cells.
#' @param tagType - \code{raw} for when tags are raw tags ; \code{cpm} when tags are tags per million ; default is \code{raw}.
#' @export
#' @return an sclData class object.
#' @examples
#' example(scluster)
sclDataConstructor <- function(tags, tagType="raw"){
    validTagTypes <- c("raw", "cpm")
    if (!(tagType %in% validTagTypes)) {
        stop("Invalid tagType parameter supplied: ", tagType, ".  Valid Tags: ", 
             paste(validTagTypes, collapse=", ")) 
    }        
    object <- new("sclData", tags=tags, tagType=tagType)
    object@sampleSize <- ncol(tags)
    return(object)
}


setGeneric("sclust", function(object, nCluster, s=NA_real_, nf=NA_integer_,
                 max.nIter=10, markers, supervise=F, alpha, seed=123, 
                 dist="euclidean", priorTPM=1){
    standardGeneric("sclust")
})

#' @title sclust
#' @rdname sclust
#' @name sclust
#' @importFrom Rcpp evalCpp
#' @import RcppParallel
#' @importFrom cidr scDataConstructor determineDropoutCandidates wThreshold
#' @importFrom cluster pam
#' @param nCluster Number of clusters.
#' @param s A tuning parameter determines the sum of weights for all features.  NOTE that only one of the parameters \code{nf} and \code{s} should be set; if neither are set, then a default value of \code{nf=300} is applied.
#' @param nf Determines the number of features (with non-zero weight) used for clustering.
#' @param max.nIter The number of iterations.
#' @param markers A list of names of marker genes, not used when supervise=FALSE.
#' @param supervise Boolean; If TRUE, implement supervised clustering.
#' @param alpha A tuning parameter determines the sum of weights of marker genes.
#' @param seed Used for random number generation reproducability.
#' @param dist Either "ecuclidean" or "cidr".  If "euclidean", then the squared euclidean distance is used, otherwise the CIDR dissimilarity measure is used.
#' @param priorTPM Count to add to deal with 0 tags for log function.
#' @export sclust
#' @return An updated sclData object with the following additional attributes:
#' 
#' \item{cluster}{A vector with the clustering result.}
#' \item{weight}{nIter x p matrix; the weight for each feature in each iteration.}
#' \item{w.bcss}{The weighted between cluster sum of squares in each iteration.}
#' \item{bcss}{nIter x p matrix; the between cluster sum of squares for each feature in each iteration.}
#' \item{totFeatDist}{1 x p matrix; the total distance calculated for each feature.}
#' \item{sumWeight}{nIter x 1 matrix; the sum of weight in each iteration.}
#' \item{nf.zw}{nIter x 1 matrix; the number of features with non-zero weight in each iteration.}
#'
#' @examples
#' example(scluster)
setMethod("sclust", "sclData", function(object, nCluster, s, nf, max.nIter, 
                                       markers, supervise, alpha, seed, 
                                       dist, priorTPM){
    # check for valid "dist" values
    validDistTypes <- c("euclidean", "cidr")
    if (!(dist %in% validDistTypes)) {
        stop("Invalid 'dist' parameter supplied: ", dist, ".  Valid dist values: ", 
             paste(validDistTypes, collapse=", ")) 
    }        
    object@nCluster <- nCluster
    object@s <- s
    object@max.nIter <- max.nIter
    object@markers <- markers
    object@supervise <- supervise
    object@alpha <- alpha
    object@seed <- seed
    object@dist <- dist
    object@priorTPM <- priorTPM
    ## small number for stopping criterion
    EPSILON <- 10^(-4)
    set.seed(seed)
    ## preprocess the input data depending on tag type (raw vs CPM), and 
    ## also do some CIDR preprocessing to allow for dropouts if the distance 
    ## measure is "cidr"
    if (dist=="cidr") {
        ## use the CIDR package to calculate an alternative distance to the default
        ## squared euclidean.  CIDR uses implicit imputation to allow for dropouts
        cData <- scDataConstructor(as.matrix(object@tags), tagType=object@tagType)
        cData <- determineDropoutCandidates(cData)
        cData <- wThreshold(cData)
        cidr_dout <- cData@dropoutCandidates    ## matrix of dropout candidate truth values
        cidr_threshold <- cData@wThreshold      ## cidr threshold
        object@tags.lcpm <- cData@nData         ## log cpm - only include where rowsums > 0
        rm(cData)   ## remove cidr object
    } else {
        tags_feat_gt_0 <- object@tags[rowSums(object@tags)>0,]
        if (object@tagType=="cpm") {
            object@tags.lcpm <- log2(tags_feat_gt_0 + priorTPM)
        } else {
            object@tags.lcpm <- log2(t(t(tags_feat_gt_0)/
                                       colSums(tags_feat_gt_0))*1000000+priorTPM)
        }
    }
    if(supervise==TRUE){
        index_markers <- match(markers, rownames(object@tags.lcpm))
        index_markers <- index_markers[!is.na(index_markers)]
        w.initial <- rep(0,nrow(object@tags.lcpm))
        w.initial[index_markers] <- 1/sqrt(length(index_markers))
    }else{
        w.initial <- rep((1/sqrt(nrow(object@tags.lcpm))), nrow(object@tags.lcpm))
        index_markers <- c(1:nrow(object@tags.lcpm))
    }
    p <- nrow(object@tags.lcpm)
    n <- ncol(object@tags.lcpm)
    ## Total distance for each feature
    totFeatDist <- matrix(rep(0), ncol=p)
    for (j in 1:n){
        if (dist=="cidr") {
            ## CIDR distance measure
            temp_d <- parallelMatrixCidrDist(object@tags.lcpm, cidr_dout, j, cidr_threshold)
        } else {
            ## squared euclidean distance
            temp_d <- parallelMatrixSquareDist(object@tags.lcpm, j)
        }
        totFeatDist <- totFeatDist + rowSums(temp_d)
    }
    ## weight for features in each iteration
    w.sparse <- matrix(NA, nrow=max.nIter, ncol=p)
    ## BCSS (Between Cluster Sum of Squares) for features in each iteration
    bcss <- matrix(NA, nrow=max.nIter, ncol=p)
    ## Weighted BCSS in each iteration
    w.bcss <- matrix(NA, nrow=max.nIter)
    ## sum of weight for each iteration
    sumW <- matrix(NA, nrow=max.nIter)
    ## the number of features in each iteration
    final_nf <- matrix(NA, nrow=max.nIter)
    ## initialize the weight for sparse clustering
    w.sparse[1,] <- w.initial
    ## calculate the sum of initial weight
    sumW[1] <- sum(w.initial)
    ## calculate the number of features initially used
    final_nf[1] <- sum(w.initial>0)
    weighted.dissim <- matrix(rep(0), nrow=n, ncol=n)
    for(nIter in 1: max.nIter){
        for (j in 1:n){
            ## calculating distance as needed
            if (dist=="cidr") {
                temp_d <- parallelMatrixCidrDist(object@tags.lcpm, cidr_dout, j, cidr_threshold)
            } else {
                temp_d <- parallelMatrixSquareDist(object@tags.lcpm, j)
            }
            weighted.dissim[j,] <- (vectorMatrixProduct(temp_d,w.sparse[nIter,])$result)[1,]
        }
        start <- proc.time()
        pam <- pam(weighted.dissim, k=nCluster, diss=TRUE)
        ## WCSS (Within Cluster Sum of Squares) for all clusters
        total.WCSS <- matrix(rep(0),nrow=nCluster, ncol=p)
        for (k in 1: nCluster){
            ## WCSS with respect to each feature in each cluster
            WCSS <- matrix(rep(0), ncol=p)
            for(j in which(pam$clustering==k)){
                ## calculating distance as needed
                if (dist=="cidr") {
                    temp_d <- parallelMatrixCidrDist(object@tags.lcpm, cidr_dout, j, cidr_threshold)
                } else {
                    temp_d <- parallelMatrixSquareDist(object@tags.lcpm, j)
                }
                WCSS <- WCSS + rowSums(as.matrix(temp_d[,which(pam$clustering==k)]))
            }
            total.WCSS[k,] <- (1/sum(pam$clustering==k))*WCSS
        }
        ## BCSS calculation
        bcss[nIter,] <- (1/n)*totFeatDist - colSums(total.WCSS)
        ## Weighted BCSS calculation
        w.bcss[nIter] <- sum(w.sparse[nIter,]*bcss[nIter,])
        ## If already up to maximum number of iteration, no need to update weight
        if(nIter==max.nIter){break}
        ## Calculate the new weights according to BCSS in last step
        w.sparse[nIter+1,] <- bcss[nIter,]/sqrt(sum(bcss[nIter,]^2)) 
        # supply a default nf value if neither nf or s are supplied
        if(is.na(nf) & is.na(s)){
            nf=300 ##Temporary default value of nf
            warning("Neither parameters 'nf' or 's' were supplied.  A default value of nf=300 has been applied")
        }
        if(is.na(nf) & !is.na(s)){
            ## Check the constraints on w for s
            if(sum(w.sparse[(nIter+1),]) < s){
                delta1 <- NULL
                diff <- bcss[nIter,]
            }else{
                delta1 <- 1
                while(sum(w.sparse[(nIter+1),]) > s){
                    diff <- bcss[nIter,] - delta1
                    diff[which(diff<0)] <- 0
                    w.sparse[nIter+1,] <- diff/sqrt(sum(diff^2))
                    delta1 <- delta1 + 1
                }
            }
        }else{
            if (!is.na(nf) & !is.na(s)) {
                stop("Set only one of the parameters 'na' and 's'")
            }
            ## Check the constraints on w for nf
            if(sum(w.sparse[(nIter+1),] > 0) < nf){
                delta1 <- NULL
                diff <- bcss[nIter,]
            }else{
                delta1 <- bcss[nIter,][(order(bcss[nIter,], decreasing=TRUE)[nf+1])]
                diff <- bcss[nIter,] - delta1
                diff[which(diff<0)] <- 0
                w.sparse[nIter+1,] <- diff/sqrt(sum(diff^2))
            }
        }
        diff3 <- matrix(NA, ncol=p)
        if(sum(w.sparse[(nIter+1),index_markers])/sum(w.sparse[nIter+1,]) >= alpha){
            delta <- NULL
        }else{
            delta <- 1
            while(sum(w.sparse[(nIter+1),index_markers])/sum(w.sparse[nIter+1,]) < alpha){
                diff1 <- diff[index_markers] + delta
                diff2 <- diff[-index_markers] - delta

                diff1[which(diff1<0)] <- 0
                diff2[which(diff2<0)] <- 0

                diff3[index_markers] <- diff1
                diff3[-index_markers] <- diff2
                w.sparse[nIter+1,] <- diff3/sqrt(sum(diff3^2))
                delta <- delta + 1
                if(sum(w.sparse[(nIter+1),-index_markers])==0){break}
            }
        }
        ## calculate the sum of weight and store in the vector
        sumW[nIter+1] <- sum(w.sparse[nIter+1,])
        ## calculate the number of features with non-zero weights
        final_nf[nIter+1] <- sum(w.sparse[nIter+1,]>0) 
        ## Check the stopping criterion
        if((sum(abs(w.sparse[nIter+1,]-w.sparse[nIter,]))/sum(w.sparse[nIter,]))< EPSILON){break}
    }
    ## Return the result
    object@cluster <- pam$clustering
    object@weight <- w.sparse[1:nIter,]
    object@w.bcss <- w.bcss[1:nIter]
    object@bcss <- bcss[1:nIter,]
    object@totFeatDist <- totFeatDist
    object@sumWeight <- sumW[1:nIter]
    object@nf.zw <- final_nf[1:nIter]
    return(object)
})


