#' @name readData
#' @description Read table file as a matrix or convert data.frame to matrix.
#' @param X Input data, which can be text table file or data.frame.
#' @param header A logical value indicating whether the file contains the names of the variables as its first line, default: TRUE.
#' @param sep The field separator character, default: '\t'.
#' @param row.names Row index of table file serves as row names, default: NULL.
#' @return df (matrix).
#' @export

readData <- function(X, header = TRUE, sep = '\t', row.names = NULL) {
    if (is.matrix(X) || is.data.frame(X)) {
        df <- X
    } else if (is.character(X)) {
        tryCatch(
            {
                df <- read.table(X, header = header, sep = sep)
            }, error = function(e) {
                println('[ERROR] Table file cannot be detected')
            }
        )
    } else {
        stop('[ERROR] X must be a matrix, data.frame or correct table file path')
    }
    if (!is.null(row.names)) {
        rownames(df) <- make.names(df[, row.names], unique = TRUE)
        df <- df[, -row.names]
    }
    return(df)
}

#' @name println
#' @description Print out message if verbose equals TRUE.
#' @param infos Message that need to be printed.
#' @param verbose Bool value, print out message or not, default: FALSE.
#' @return NULL.
#' @export

println <- function(X, verbose = FALSE, ...) {
    infos <- do.call(sprintf, c(list(X, ...)))
    if (verbose) cat(paste0(infos, '\n'))
}

#' @name addNames
#' @description Add row and column names of the object.
#' @param X Data.frame or matrix of the input data.
#' @param row.names A vector of row names, default: NULL.
#' @param col.names A vector of column names, default: NULL.
#' @return Data.frame that added names.
#' @export

addNames <- function(X, row.names = NULL, col.names = NULL) {
    if (!is.null(row.names)) rownames(X) <- row.names
    if (!is.null(col.names)) colnames(X) <- col.names
    return(X)
}

#' @name calcLogFC
#' @description Calculate log-fold change distance of vector
#' @param vals A vector of values.
#' @return Log-fold change distance result.
#' @export

calcLogFC <- function(vals) {
    max.idx <- which.max(vals)
    max.val <- vals[max.idx]
    avg.other.vals <- mean(vals[-max.idx])
    logFC <- log2(max.val / avg.other.vals)
    return(list(logFC = logFC, idx = max.idx))
}

#' @name selectSubsets
#' @description Select top N rows of each subsets from data list
#' @param data.lst Data list that includes various phenotype class subsets.
#' @param top.num Top N rows need to be selected.
#' @param fields Fields of each subsets.
#' @return Selected subsets from data list.
#' @export

selectSubsets <- function(data.lst, top.num, fields = NULL) {
    selected.sets <- c()
    for (name in names(data.lst)) {
        tmp.df <- data.lst[[name]]
        if(is.null(fields)) fields <- colnames(tmp.df)
        selected.sets <- rbind(selected.sets, tmp.df[1 : top.num, fields])
    }
    return(selected.sets)
}

#' @name resortMatrix
#' @description Resort expression matrix by cell names.
#' @param X Gene expression data matrix.
#' @param maxFC.idx Maximum gene expression index.
#' @param cells sorted cell names.
#' @return Sorted gene expression matrix.
#' @export

resortMatrix <- function(X, maxFC.idx, cells) {
    results <- c()
    for (cell in cells) {
        idx <- which(cells == cell)
        sub.names <- names(maxFC.idx[maxFC.idx == idx])
        results <- rbind(results, X[sub.names, ])
    }
    return(data.frame(results))
}
