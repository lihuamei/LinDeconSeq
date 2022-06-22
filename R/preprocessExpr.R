#' @name prerpocessExpr
#' @description Expression data of each cell type are intergrated into a single data matrix, which rows are genes and columns are cell typs.
#' @param X Normalized gene expression profile of each cell type, which rows represent genes and columns represent pure samples.
#' @param phenotypes Phenotype classes. value 1 = indicate membership of the reference sample,
#' value 2 = indicates the class that the sample will be compared against.
#' @param method Merge replicates by specified method, default: mean.
#' @param cv.cutoff Filtered out CV greater than threshold genes, default: 1.5.
#' @return list of normalized expression matrix and weights [list]
#' @export

prerpocessExpr <- function(X, phes, method = c('mean', 'median'), cv.cutoff = 1.5) {
    ref.grouped  <- c()
    if (method == 'mean') {
        ref.grouped <- X %*% t(phes) %*% diag(1.0 / rowSums(phes))
    } else {
        for (cell in rownames(phes)) {
            ref.grouped <- X[, phes[cell, ] == 1] %>% (matrixStats::rowMedians) %>% cbind(ref.grouped, .)
        }
    }
    ref.grouped <- addNames(ref.grouped, col.names = rownames(phes), row.names = rownames(X))
    ref.grouped <- ref.grouped[matrixStats::rowSds(ref.grouped) != 0, ]
    X.sub <- X[rownames(ref.grouped), ]
	
    max.indexes <- apply(ref.grouped, 1, which.max)
    cvs <- sapply(1 : dim(X.sub)[1], function(ind) {
	    tmp.v <- X.sub[ind, which(phes[max.indexes[ind], ] == 1)]
        cv <- sd(tmp.v) / mean(tmp.v) 
    })
    ref.grouped <- ref.grouped[which(cvs < cv.cutoff), ]
    return(ref.grouped)
}
