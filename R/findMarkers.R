##########################################################################################################
#' @description Search marker genes on the basis of seed genes for each cell type.
#' @param ref.expr Refernce gene expression profile.
#' @param query.genes A list of query genes.
#' @param scale.data Scale gene expression or not, default: TRUE.
#' @param p.cut Threshold for filtering cell type-specific markers. Default: 0.1.
#' @return Seurat object 
#' @export 

##########################################################################################################

searchMarkersByCorr <- function(ref.expr, query.genes, scale.data = TRUE, p.cut = 0.1) {
	if (!inherits(query.genes, 'list')) stop(sprintf('%s must be a list of query genes, exiting...', deparse(substitute(query.genes))))
	if (scale.data) {
		ref.expr <- sweep(ref.expr, 1, apply(ref.expr, 1, mean), "-")
		ref.expr <- sweep(ref.expr, 1, sqrt(apply(ref.expr, 1, var)), "/")
	}
	ref.expr <- ref.expr[rowSums(is.na(ref.expr)) == 0, ]
	p1 <- prcomp_irlba(ref.expr %>% t, n = min(length(query.genes) * 2, min(dim(ref.expr))))
	data.u <- p1$rotation %>% as.data.frame %>% t
	colnames(data.u) <- rownames(ref.expr)
	marker.lst <- lapply(query.genes, function(QGs) {
		QGs <- names(QGs) 
		cov.mat <- t(data.u[, QGs]) %*% data.u
		diags <- sqrt(colSums((data.u)^2))
		cor.mat <- cov.mat / diags[QGs]
		cor.mat <- sweep(cor.mat, 2, diags, '/')
		colnames(cor.mat) <- colnames(data.u)
		rownames(cor.mat) <- QGs
		cor.vec <- colMeans(atanh(cor.mat[seq_len(length(QGs)), -c(which(colnames(cor.mat) %in% QGs))]))
		cor.vec[order(cor.vec) %>% rev]		
		pvals <- pnorm(cor.vec, mean = mean(cor.vec), sd = sd(cor.vec), lower.tail = FALSE)
		c(QGs, pvals[pvals < p.cut] %>% names) %>% unique
	})
	return(marker.lst)
}


