#' @name deconSeq
#' @title Deconvolution using w-RLM method.
#' @description Cell type deoovolution for bulk samples based on weighted RLM method.
#' @param bulks bulk gene expression (row = Genes, column = Samples).
#' @param signature Signature matrix, i.e. cell type-specific gene expression matrice of cell types.
#' @param weight Weight gene or not, default: TRUE.
#' @param intercept Add intercept when using robust linear regression, default: TRUE.
#' @param scale Scale bulk gene expression data or not, default: FALSE.
#' @param QN Quantile normalize bulk expression profile or not, default: FALSE
#' @param verbose logical, to print the detailed information for each iteration.
#' @return A data.frame of fractions for bulk samples.
#' @export deconSeq

deconSeq <- function(bulks, signature, weight = TRUE, intercept = TRUE, scale = FALSE, QN = FALSE, verbose = TRUE) {
    if (missing(bulks) || missing(signature)) {
        stop("[ERROR] bulk and pures must be fed simultaneously")
    }
    if (class(signature) != 'data.frame' && class(signature) != 'matrix') {
        stop('[ERROR] Signature must be data.frame or matrix, exiting...')
    }
    if (class(bulks) != 'data.frame' && class(bulks) != 'matrix') {
        stop('[ERROR] bulks must be data.frame or matrix, exiting...')
    }
	if (max(signature) < 50) signature <- 2^signature
    if (max(bulks) < 50) bulks <- 2^bulks
    if (min(signature) == 0) signature <- as.matrix(signature + 1)
    if (min(bulks) == 0) bulks <- bulks + 1
    println('[INFO] There are %d bulk samples need to be deconvoluted', verbose, dim(bulks)[2])
	
	if (QN) {
		xx <- bulks
        bulks <- normalize.quantiles(bulks)
        bulks <- addNames(bulks, rownames(xx), colnames(xx))
    }
    if (scale) bulks <- scale(bulks)
	
	ovp.genes <- intersect(rownames(signature), rownames(bulks))
    signature <- signature[ovp.genes, ] %>% as.matrix
    bulks <- bulks %>% .[ovp.genes, ]

    est.res <- apply(bulks, 2, function(bulk) {
		intercept <- ifelse(intercept, 1, 0) 
		if (weight) {
			status <- try(w.coefs <- solveDampenedRLM(signature, bulk, intercept), silent = TRUE)
		} else {
			status <- 'Okay'
		}
		if ((class(status) == 'try-error') || (!weight)) {
			if (intercept) {
				w.coefs <- coef(rlm(bulk ~ signature, maxit = 10000))
			} else {
				w.coefs <- coef(rlm(bulk ~ signature - 1, maxit = 10000))			
			}
		}
		w.coefs[w.coefs < 0] <- 0
		coefs <- w.coefs[(intercept + 1) : length(w.coefs)] / sum(w.coefs[(intercept + 1) : length(w.coefs)])
    })

    fractions <- addNames(est.res, row.names = colnames(signature), col.names = colnames(bulks)) %>% t
    if (verbose) print(fractions)
    return(fractions %>% data.frame)
}
