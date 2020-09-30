################################################################################
#' @description Calculate TPM values of genes geneated from high-throughput sequencing.
#' @param X Input data, which must be matrix or data.frame.
#' @return TMP matrix.
#' @export calcTpm
################################################################################

calcTpm <- function(X) {
    load(system.file('inst', 'hg19_transcripts_length.Rda', package = 'LinDeconSeq'))
    genes <- intersect(rownames(gene_lengths), rownames(X))
    X <- X[genes, ]; gene.lengths <- as.numeric(as.vector(gene_lengths[genes, 'LENGTH']))
    rate <- X / gene.lengths
    X.norm <- addNames(t(t(rate) / colSums(rate)) * 1e6, rownames(X), colnames(X))
    return(X.norm)
}

################################################################################
#' @description Normalize counts data of gene expression matrix by DESeq2 method.
#' @param X Input data, which must be matrix or data.frame.
#' @param phes Phenotype class of all samples.
#' @return Normalized matrix.
################################################################################

deseq2Norm <- function(X, phes) {
    coditions.df <- data.frame(name = colnames(phes), condition = 1, row.names = colnames(phes))
    for(cell in rownames(phes)) {
        coditions.df[which(phes[cell, ] == 1), 'condition'] <- cell
    }
    coditions.df[, 'condition'] <- factor(coditions.df[, 'condition'], levels = rownames(phes))
    dds <- DESeqDataSetFromMatrix(X, colData = coditions.df, design = ~condition)
    dds <- estimateSizeFactors(dds)
    X.norm <- counts(dds, normalized = TRUE)
    return(X.norm)
}

################################################################################
#' @description Normalize counts data of gene expression matrix by DESeq2 method.
#' @param X Input data, which must be matrix or data.frame.
#' @return Normalized matrix.
################################################################################

VST2Norm <- function(X) {
    X.norm <- varianceStabilizingTransformation(X)
    return(X.norm)
}

################################################################################
#' @description Normalize data with specified method.
#' @param X Input data, which must be matrix or data.frame.
#' @param phes Phenotype class of all samples, default: NULL.
#' @param method Specified normalize method, support QN, TPM and None. default: QN.
#' @param data.type Data type, RNASeq, MA or PEAK, default: RNASeq.
#' @return Normalized data (matrix).
#' @export normalizeData
################################################################################

normalizeData <- function(X, phes =NULL, method = c('QN', 'NONE', 'TPM', 'DESeq2', 'CPM', 'VST'), data.type = 'RNASeq') {
    if (missing(X)) stop("[ERROR] argument 'X' is misiing, with no default")
    X <- as.matrix(X)
    if ('QN' == method) {
        X.norm <- preprocessCore::normalize.quantiles(X)   
        X.norm <- addNames(X.norm, rownames(X), colnames(X))
    } else if ('TPM' == method) {
        X.norm <- calcTpm(X)
    } else if ('DESeq2' == method) {
        X.norm <- deseq2Norm(X, phes)
    } else if ('CPM' == method) {
        X.norm <- cpm(X)
    } else if ('VST' == method) {
        X.norm <- VST2Norm(X)
    } else if ('NONE' == method) {
        X.norm <- X
    } else {
        stop(sprintf('[ERROR] Unused argument in method (method = %s)', method))
    }  
    return(X.norm + 1)
}
