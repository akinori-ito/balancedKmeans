#library(clue)

class_seq <- function(n,k) {
  if (n < k) {
    stop("Number of samples should be larger than the number of clusters")
  }
  mag <- ceiling(n/k)
  cseq <- rep(1:k,each=mag)
  ind <- sort(sample(1:length(cseq),n))
  cseq[ind]
}

#'The balanced k-means algorithm
#'
#'@param x a matrix or data frame to be clustered, each row represents an item
#'@param k number of clusters
#'@param iter.max maximum number of iteration
#'@return A list of following items:
#'cluster: cluster numbers of the each item
#'centers: the centroids
#'size: size of the clusters
#'iter: final number of iteration
#'@export
balanced_kmeans <- function(x,k=2,iter.max=10) {
  n <- nrow(x)
  if (n < k) {
    stop("Number of samples should be larger than the number of clusters")
  }
  centers <- x[sample(1:n,k),]
  iter <- 0
  dmat <- matrix(0,nrow=n,ncol=n)
  cseq <- class_seq(n,k)
  cluster <- rep(0,n)
  while (iter < iter.max) {
    for (i in 1:n) {
      for (j in 1:n) {
        dmat[i,j] <- sum((x[i,]-centers[cseq[j],])^2)
      }
    }
    a <- as.numeric(clue::solve_LSAP(dmat))
    newcluster <- cseq[a]
    if (all(newcluster == cluster)) {
      break
    }
    cluster <- newcluster
    for (i in 1:k) {
      centers[i,] <- colMeans(x[cluster==i,])
    }
    iter <- iter+1
  }
  list(cluster=cluster,
       centers=centers,
       size=as.numeric(table(cluster)),
       iter=iter)
}
