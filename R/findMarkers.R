################################################################################
#' @description Identify high variable genes and filter out
#' @param X Expression matrix data of samples with the same cell type.
#' @return Shannon entropy data.frame

################################################################################

.shannonEntropy <- function(X) {
    nrep.sns  <- dim(X)[2]
    X.norm <- X / rowSums(X)
    X.mu <- rowMeans(X.norm)
    weight.genes <- matrixStats::rowMaxs(X)
    weight.genes <- tanh(0.1 * weight.genes / median(weight.genes))
    shannon.dist <- 1 / nrep.sns * rowSums(X.norm / X.mu * log2(X.norm / X.mu)) * weight.genes
    return(shannon.dist)
}
.calcScores <- function(x, seed.genes, cell, idx) {
    tmp.w <- log2(x[cell] / mean(x[-idx]))
    cor.v <- cor(x, seed.genes[[cell]])
    sign  <- ifelse(cor.v > 0, 1, -1)
    score <- cor.v^2 * 1 / (1 + exp(-tmp.w)) * sign
    return(score)
}

.fitNorm <- function(x) {
	dor <- density(x, kernel = 'gaussian')
    dist.max <-dor$x[which.max(dor$y)]
    dist.offset <- x- dist.max
    tmp.dist <- c(dist.offset[dist.offset <= 0], abs(dist.offset[dist.offset < 0])) + dist.max
    dist.norm.fit <- fitdistr(tmp.dist, 'normal')
    est.params <- as.vector(dist.norm.fit$estimate)
	return(est.params)
}

################################################################################
#' @name deriveDEGenes
#' @description Derive cell type-specific genes based Shannon Entropy distance.
#' @param X Expression data of each cell type are intergrated into a single data matrix.
#' @param q.cut Threashold q.value to filter out cell type-specific genes, default: 0.01.
#' @param verbose Print out filter information or not, default: FALSE.
#' @return X.filtered

################################################################################

deriveDEGenes <- function(X, q.cut = 0.01, perm = 100, verbose = FALSE) {
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
    logFC.idx <- apply(X.degenes, 1, calcLogFC) %>% do.call(rbind.data.frame, .) %>% .[.$logFC < 10, ] %>%  as.matrix
    X.degenes <- X.degenes[rownames(logFC.idx), ]
    println('[INFO] %d cell type specific genes have been detected with q.cut %f', verbose, dim(X.degenes)[1], q.cut)
    return(list(DEG = X.degenes, adj.p = q.vals[rownames(X.degenes)], logFC = logFC.idx[, 1], topIdx = logFC.idx[, 2], ShannonDist = shannon.dist))
}

################################################################################
#' @description Assign cell type specific genes to each cell type.
#' @param X Cell type specific gene expression matrix.
#' @param bg.genes.expr Background gene expression for computing random scores.
#' @param verbose  Print out information or not, default: FALSE.
#' @return A list of cell markers.
#'
################################################################################

assignedCellMarkers <- function(X, bg.genes.expr, p.cut = 0.05, do.perm = 10, verbose = FALSE) {  
    set.seed(2019)	
	rand.scores <- list(); seed.genes <- list()
	bak.genes.expr <- bg.genes.expr
    cell.names <- colnames(X$DEG)[sort(unique(X$topIdx))]
    
	println('[INFO] Select seed genes and random permutating...', verbose = verbose)	
	for (i in 1 : do.perm) {
		tmp.genes.expr <- t(apply(bak.genes.expr, 1, sample)) 
		bg.genes.expr  <- rbind(bg.genes.expr, tmp.genes.expr)
	}
	rownames(bg.genes.expr) <-  paste0('Gene.', 1 : dim(bg.genes.expr)[1])
	bg.genes.expr <- bg.genes.expr[sample(rownames(bg.genes.expr)), ]
	for (idx in unique(X$topIdx)) {        
		tar.cell  <- colnames(X$DEG)[idx]
		tmp.genes <- X$logFC[X$topIdx == idx] %>% sort %>% rev %>% names(.)
		seed.genes[[tar.cell]] <- X$DEG[tmp.genes[1], ]
		random.exprs <- bg.genes.expr[sample(rownames(bg.genes.expr), 10000), ]
		random.exprs <- random.exprs[matrixStats::rowSds(as.matrix(random.exprs)) != 0, ]
		rand.scores[[tar.cell]] <- apply(random.exprs, 1, function(x) {
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
    assigned.df <- assigned.df[rowSums(assigned.df) > 0, ]
    X.marker.filtered <- resortMatrix(X$DEG[rownames(assigned.df), ], X$topIdx[rownames(assigned.df)], cell.names)
    scores.df <- scores.df[rownames(X.marker.filtered), ]
    println('[INFO] %d low confidence marker genes have been filtered out.', verbose = verbose, dim(X$DEG)[1] - dim(X.marker.filtered)[1])
    
    return(list(
        score = scores.df,
        markers = assigned.df[rownames(X.marker.filtered), ],
        marker.expr = X.marker.filtered,
        logFC       = X$logFC[rownames(scores.df)],
        raw.marker = X$DEG)
    )
}

################################################################################
#' Optimize the number of cell type-specific genes of each phnotype class to ensure the stability of the signature matrices.
#' @param X A list of cell markers that generated by the above step (assignedCellMarkers).
#' @param cell.names Target cell names of each subset.
#' @param min.group Minimum number of each phnotype class, default: 50.
#' @param max.group Maximum number of each phnotype class, default: 200.
#' @return Derived signature matrices.

################################################################################

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

################################################################################
#' @title Derive signature matrices.
#' @description Derive signature marker genes by a series steps, including normalizing, filtered low confidence genes and DEG analysis.
#' @param refs Gene expression profile of each cell type, which rows represent genes and columns represent pure samples.
#' @param phes Phenotype classes. value 1 = indicate membership of the reference sample,
#' value 2 = indicates the class that the sample will be compared against.
#' @param min.group Minimum number of each phnotype class, default: 50.
#' @param max.group Maximum number of each phnotype class, default: 200.
#' @param norm.method Select method to normalize reference profile [QN, TPM, CPM, DESeq2 or NONE], default: TPM (Transcript per million).
#' @param data.type Data type, RNASeq, MA or PEAK, default: RNASeq.
#' @param verbose logical, to print the detailed information.
#' @export findMarkers
#' @return list
#'
################################################################################

findMarkers <- function(refs, phes, min.group = 50, max.group = 200, norm.method = 'TPM', q.cut = 0.01, p.cut = 0.05, data.type = 'RNASeq', verbose = TRUE) {
    if (data.type == 'PEAK') {
        rownames(refs) <- paste(refs[, 1], refs[, 2], refs[, 3], sep = '_')
        refs <- refs[, -c(1, 2, 3)]
    } else if (data.type != 'RNASeq' && data.type != 'MA' && data.type != 'PEAK') {
        stop('[ERROR] Data.type only supports gene and peak, exiting...')
    }
    if (max(refs) < 50) refs <- 2^refs
    phes[phes == 2] <- 0
    println('[INFO] There are %d samples and %d genes in the reference profile', verbose, dim(refs)[2], dim(refs)[1])
    refs.norm <- refs[rowSums(refs) > 0, ] %>% normalizeData(., phes, norm.method, data.type = data.type)

    ref.grouped   <- prerpocessExpr(refs.norm, phes, method = 'median')	
    ctsgs.infos   <- deriveDEGenes(ref.grouped, q.cut = q.cut, verbose = verbose)
    bg.genes.expr <- ref.grouped[!(rownames(ref.grouped) %in% rownames(ctsgs.infos$DEG)), ]
    cell.markers  <- assignedCellMarkers(ctsgs.infos, bg.genes.expr, p.cut = p.cut, verbose = verbose)
    println('[INFO] Optimizing cell type-specific genes to derive signature matrix...', verbose)
    sig.marker.infos <- optimizeSignatures(
        cell.markers         ,
        colnames(ref.grouped),
        min.group = min.group,
        max.group = max.group
    )
    println('[INFO] Group size -> %d, condition number -> %f', verbose, sig.marker.infos$group.size, min(sig.marker.infos$kappa))
    return(list(
		cellMarkers = cell.markers[c('score', 'markers', 'marker.expr')], 
		sigMatrix   = sig.marker.infos
	))
}
