#' @name .shannonEntropy
#' @description Identify high variable genes and filter out
#' @param X Expression matrix data of samples with the same cell type.
#' @return Estimated shannon entropy (data.frame)
#' @export 

.shannonEntropy <- function(X) {
    nrep.sns  <- dim(X)[2]
    X.norm <- X / rowSums(X)
    X.mu <- rowMeans(X.norm)
    weight.genes <- rowMaxs(X)
    weight.genes <- tanh(0.1 * weight.genes / median(weight.genes))
    shannon.dist <- 1 / nrep.sns * rowSums(X.norm / X.mu * log2(X.norm / X.mu)) * weight.genes
    return(shannon.dist)
}

#' @name .calcScores
#' @description Calculated co-linearity score based on seed genes.
#' @param cell Cell type name.
#' @param idx Index number of cell type.
#' @return Calculated score
#' @export 

.calcScores <- function(x, seed.genes, cell, idx) {
    tmp.w <- log2(x[cell] / mean(x[-idx]))
    cor.v <- cor(x, seed.genes[[cell]])
    sign  <- ifelse(cor.v > 0, 1, -1)
    score <- cor.v^2 * 1 / (1 + exp(-tmp.w)) * sign
    return(score)
}

#' @name .fitNorm
#' @description Fit normal distribution for shannon-entropy distance.
#' @param x A vector of shannon-entropy distance.
#' @return Estimated mu and sd of normal distribution (list)
#' @export

.fitNorm <- function(x) {
	dor <- density(x, kernel = 'gaussian')
    dist.max <-dor$x[which.max(dor$y)]
    dist.offset <- x- dist.max
    tmp.dist <- c(dist.offset[dist.offset <= 0], abs(dist.offset[dist.offset < 0])) + dist.max
    dist.norm.fit <- fitdistr(tmp.dist, 'normal')
    est.params <- as.vector(dist.norm.fit$estimate)
	return(est.params)
}

#' @name deriveDEGenes
#' @description Derive cell type-specific genes based Shannon Entropy distance.
#' @param X Expression data of each cell type are intergrated into a single data matrix.
#' @param q.cut Threashold q.value to filter out cell type-specific genes, default: 0.01.
#' @param perm Number of permutation, default: 100.
#' @return X.filtered
#' @export

deriveDEGenes <- function(X, q.cut = 0.01, perm = 100) {
	set.seed(2019)
    shannon.dist  <- .shannonEntropy(X)
    random.params <- sapply(1 : perm, function(pcnt) {
        random.mat <- matrix(runif(dim(X)[1] * dim(X)[2], min = min(X), max = max(X)), nrow = dim(X)[1], ncol = dim(X)[2])
        shannon.dist.random <- .shannonEntropy(random.mat)
        est.params <- .fitNorm(shannon.dist.random)
        return(est.params)
    }) %>% t %>% colMeans
    
	p.raw  <- pnorm(shannon.dist, mean = random.params[1], sd = random.params[2], lower.tail = FALSE)
    q.vals <- p.adjust(p.raw, method = 'BH')
    X.degenes <- X[names(q.vals[q.vals <= q.cut]), ]
    logFC.idx <- apply(X.degenes, 1, calcLogFC) %>% do.call(rbind.data.frame, .) %>% as.matrix
    X.degenes <- X.degenes[rownames(logFC.idx), ]
    
	return(
		list(
			DEG = X.degenes, 
			adj.p = q.vals[rownames(X.degenes)], 
			logFC = logFC.idx[, 1], 
			topIdx = logFC.idx[, 2], 
			shannonDist = shannon.dist
		)
	)
}

#' @name assignedCellMarkers
#' @description Assign cell type specific genes to each cell type.
#' @param X Cell type specific gene expression matrix.
#' @param bg.genes.expr Background gene expression for computing random scores.
#' @param verbose  Print out information or not, default: FALSE.
#' @return A list of cell markers.
#' @export

assignedCellMarkers <- function(X, bg.genes.expr, p.cut = 0.05, verbose = FALSE) {  
    set.seed(2019)	
	rand.scores <- seed.genes <- list()
    cell.names <- colnames(X$DEG)[sort(unique(X$topIdx))]
    
	println('[INFO] Select seed genes and random permutating...', verbose = verbose)	
	bg.genes.expr <- matrix(runif(10000 * ncol(bg.genes.expr), min = min(bg.genes.expr), max = max(bg.genes.expr)), nrow = 10000, ncol = ncol(bg.genes.expr))
	rownames(bg.genes.expr) <- paste0('Gene.', 1 : dim(bg.genes.expr)[1])
	colnames(bg.genes.expr) <- colnames(X$DEG) 

	for (idx in unique(X$topIdx)) {        
		tar.cell  <- colnames(X$DEG)[idx]
		tmp.genes <- X$logFC[X$topIdx == idx] %>% sort %>% rev %>% names(.)
		seed.genes[[tar.cell]] <- X$DEG[tmp.genes[1], ]
		random.exprs <- bg.genes.expr[sample(rownames(bg.genes.expr), 5000), ]
		random.exprs <- random.exprs[rowSds(as.matrix(random.exprs)) != 0, ]
		rand.scores[[tar.cell]] <- apply(random.exprs %>% data.frame, 1, function(x) {
			score <- .calcScores(x, seed.genes, tar.cell, idx)
		}) %>% sort
	}
	rm(bg.genes.expr)

    println('[INFO] Assigning cell type-specific genes to each cell subset', verbose)
    pb <- progress_bar$new(total = dim(X$DEG)[1], clear = TRUE)
    results <- sapply(rownames(X$DEG), function(gene) {
        pb$tick()
        score.pval <- sapply(cell.names, function(cell) {
            idx    <- which(cell.names == cell)
            score  <- .calcScores(X$DEG[gene, ], seed.genes, cell, idx)
            pval   <- 1 - (which.min(abs(rand.scores[[cell]] - score)) / length(rand.scores[[cell]]))
            return(c(score, pval))
        })
    }) %>% t %>% data.frame

    assigned.df <- pvals.df <- addNames(
        results[, seq(2, dim(results)[2], 2)],
        row.names = rownames(X$DEG),
        col.names = cell.names
    )
    scores.df <- addNames(results[, seq(1, dim(results)[2], 2)], row.names = rownames(X$DEG), col.names = cell.names)
    assigned.df[pvals.df < p.cut] <- 1; assigned.df[pvals.df >= p.cut] <- 0
    X.marker.filtered <- resortMatrix(X$DEG[rownames(assigned.df), ], X$topIdx[rownames(assigned.df)], cell.names)
    scores.df <- scores.df[rownames(X.marker.filtered), ]
    return(
		list(
			score       = scores.df,
			markers     = assigned.df[rownames(X.marker.filtered), ],
			marker.expr = X.marker.filtered,
			logFC       = X$logFC[rownames(scores.df)],
			raw.marker  = X$DEG
		)
    )
}

#' @name optimizeSignatures
#' @description Optimize the number of cell type-specific genes of each phnotype class to ensure the stability of the signature matrices.
#' @param X A list of cell markers that generated by the above step (assignedCellMarkers).
#' @param cell.names Target cell names of each subset.
#' @param min.group Minimum number of each phnotype class, default: 50.
#' @param max.group Maximum number of each phnotype class, default: 200.
#' @return Derived signature matrices.
#' @export

optimizeSignatures <- function(X, cell.names, min.group = 50, max.group = 200) {
    X$markers <- sapply(rownames(X$markers), function(gene) {
        tmp.lab <- X$markers[gene, ]
        if (sum(tmp.lab) > 1) {
            idxes <- which(tmp.lab == 1)
            idx2  <- which.max(X$marker.expr[gene, idxes])
            tmp.lab[tmp.lab == 1] <- 0
            tmp.lab[idxes[idx2]]  <- 1
        }
        return(tmp.lab)
    }) %>% as.matrix %>% t

    X <- sapply(cell.names, function(cell) {
        genes.df <- X$markers %>% as.matrix %>% .[, cell]
        genes <- names(genes.df[genes.df == 1])
        sorted.lfc <- X$logFC[genes] %>% sort %>% rev %>% names     
        return(list(X$marker.expr[sorted.lfc, ]))
    })

    iter.number <- sapply(names(X), function(cell) { dim(X[[cell]])[1] }) %>% as.vector %>% min(min(.), max.group)
    if (iter.number > min.group) {
        kapppa.vals <- sapply(min.group : iter.number, function(idx) {
            as.matrix(selectSubsets(X, idx, cell.names)) %>% kappa
        })
    } else {
        kapppa.vals <- as.matrix(selectSubsets(X, iter.number, cell.names)) %>% kappa
    }
    group.size <- min(min.group, iter.number) + which.min(kapppa.vals) - 1
    sigmat.res <- selectSubsets(X, group.size, cell.names)
    return(list(sig.mat = sigmat.res, kappa = kapppa.vals, group.size = group.size))
}

#' @name findMarkers
#' @title Derive signature matrices based on purified gene expression profile.
#' @description Derive signature marker genes by a series steps, including normalizing, filtered low confidence genes and DEG analysis.
#' @param refs Gene expression profile of each cell type, which rows represent genes and columns represent pure samples.
#' @param phes Phenotype classes. Cell type (Row) x Sample name (column); value 1 = indicate membership of the reference sample; Otherwise, value = 0. 
#' @param QN Quantile normalization for expression profile or not, default: TRUE.
#' @param q.cut Cutoff q-value for filtering candidate marker genes. default: 0.01.
#' @param p.cut Cutoff p-value for assigning marker genes to cell type, default: 0.1.
#' @opt.sigmat Optimizing singature matrix for deconvolution based on minimum kappa value, default: TRUE.
#' @param min.group Minimum number of each phnotype class when optimizing signature matrix, default: 50.
#' @param max.group Maximum number of each phnotype class when optimizing signature matrix, default: 200.
#' @param QN Quantile normalization for expression profile or not, default: TRUE.
#' @param verbose logical, to print the detailed information.
#' @return A list of derived results.
#' @export findMarkers

findMarkers <- function(refs, phes, QN = TRUE, q.cut = 0.01, p.cut = 0.1, opt.sigmat = TRUE, min.group = 50, max.group = 200, verbose = TRUE) {
	if (!is.matrix(refs) && !is.data.frame(refs)) stop(sprintf("'%s' needs to be given as a matrix or data.frame", deparse(substitute(refs))))
	if (!is.matrix(phes) && !is.data.frame(phes)) stop(sprintf("'%s' needs to be given as a matrix or data.frame", deparse(substitute(phes))))
	if (nrow(refs) < 5000) println('[WARN] Number of genes less than 5000, the statistical power may not be sufficien.', println)
	if (nrow(phes) < 2) stop(println('[ERROR] Number of cell types must be greater than 1.'))
	if (opt.sigmat) {
		if (min.group > max.group) println("[WARN] 'min.group' size should be less than 'max.group' szie.")
	    if (min.group <= 0) stop(sprintf("min.group = %s must be greater than 0", deparse(substitute(min.group))))
		if (max.group <= 0) stop(sprintf("max.group = %s must be greater than 0", deparse(substitute(max.group))))
	}
	
	if (max(refs) < 50) refs <- 2^refs
    println('[INFO] %d samples and %d genes in the reference profile', verbose, dim(refs)[2], dim(refs)[1])
	
	rownames(phes) <- gsub('-', '.', rownames(phes))
	ovp.sns <- intersect(colnames(phes), colnames(refs))
	if (length(ovp.sns) > 0) {
		phes <- phes[, ovp.sns]
		refs <- refs[, ovp.sns]
	} else {
		stop('No overlapping samples between expression profile and phenotype data.')
	}
	refs <- refs[rowSums(refs) > 0, ] %>% as.matrix
		
	# Quantile normalization or not. 
	if (QN) {
		refs.norm <- normalize.quantiles(refs)
		refs.norm <- addNames(refs.norm, rownames(refs), colnames(refs))
	} else refs.norm <- refs
	rm(refs)

    ref.grouped <- prerpocessExpr(refs.norm, phes %>% as.matrix, method = 'mean', cv.cutoff = 2.0)
    ref.grouped[ref.grouped == 0] <- mean(ref.grouped)
	ctsgs.infos <- deriveDEGenes(ref.grouped, q.cut = q.cut) 
	println('[INFO] %d candidate cell type-specific genes detected with q.cut %g', verbose, nrow(ctsgs.infos$DEG), q.cut)

	bg.genes.expr <- ref.grouped[!(rownames(ref.grouped) %in% rownames(ctsgs.infos$DEG)), ]
    cell.markers  <- assignedCellMarkers(ctsgs.infos, bg.genes.expr, p.cut = p.cut, verbose = verbose)
   
	markers <- cell.markers$markers[rowSums(cell.markers$markers) > 0, ]
	marker.expr <- cell.markers$marker.expr[rownames(markers), ]
	println('[INFO] %d low confidence marker genes have been filtered out.', verbose = verbose, dim(cell.markers$markers)[1] - dim(markers)[1])

	if (opt.sigmat) {
		println('[INFO] Optimizing cell type-specific genes to derive signature matrix...', verbose)
		sig.marker.infos <- optimizeSignatures(
			cell.markers         ,
			colnames(ref.grouped),
			min.group = min.group,
			max.group = max.group
		)
		println('[INFO] Group size -> %d, condition number -> %f', verbose, sig.marker.infos$group.size, min(sig.marker.infos$kappa))
	} else {
		sig.marker.infos <- NULL
	}
	return(list(
		cellMarkers = list(markers = markers, marker.expr = marker.expr), 
		sigMatrix   = sig.marker.infos
	))
}
