##########################################################################################################
#' @title println
#' @description Print out message if verbose equals TRUE.
#' @param infos Message that need to be printed.
#' @param verbose Bool value, print out message or not, default: FALSE.
#' @return NULL.
#' @export

println <- function(X, verbose = FALSE, ...) {
    infos <- do.call(sprintf, c(list(X, ...)))
    if (verbose) cat(paste0(infos, '\n'))
}

##########################################################################################################
#' @title createSeuratObjAndQC
#' @param ref.expr Reference gene expression profile.
#' @param ref.anno Cell types corrsponding to the single cells in the reference.
#' @return Seurat object.

##########################################################################################################

createSeuratObjAndQC <- function(ref.expr, ref.anno) {
	ref.obj <- CreateSeuratObject(
		counts = ref.expr, 
		project = 'Reference', 
		meta.data = ref.anno %>% as.data.frame(.) %>% `colnames<-`('CellType') %>% `rownames<-`(colnames(ref.expr))
	)
	Idents(ref.obj) <- ref.anno
	ref.obj[["percent.mt"]] <- PercentageFeatureSet(ref.obj, pattern = "^MT-")
	ref.obj[["percent.rb"]] <- PercentageFeatureSet(ref.obj, pattern = "^RP[SL]")
	ref.obj <- subset(
		ref.obj, 
		subset = ( 
			nFeature_RNA > min(quantile(ref.obj@meta.data$nFeature_RNA, 0.05), 300) & 
			nFeature_RNA < max(quantile(ref.obj@meta.data$nFeature_RNA, 0.95), 5000) & 
			nCount_RNA   > min(quantile(ref.obj@meta.data$nCount_RNA, 0.05), 500) & 
			nCount_RNA   < max(quantile(ref.obj@meta.data$nCount_RNA, 0.95), 10000) & 
			percent.mt   < max(quantile(ref.obj@meta.data$percent.mt, 0.95), 15) & 
			percent.rb   < max(quantile(ref.obj@meta.data$percent.rb, 0.95), 35)
		)
	)
	return(ref.obj)
}

##########################################################################################################
#' @title identSeedGenes
#' @param avg.expr Reference gene expression profile.
#' @param factor.size Factor size to scale the weight. Default: 0.1.
#' @param count Number of candidate query genes. Default: 10.
#' @return A list of seed genes for cell types.
#' @export 
#' @example	

##########################################################################################################

identSeedGenes <- function(avg.expr, factor.size = 0.1, count = 10) {
	genes.weight <- apply(avg.expr, 1, median) / median(apply(avg.expr, 1, median))
	genes.weight <- tanh(factor.size * genes.weight)
	seed.genes <- lapply(1 : ncol(avg.expr), FUN = function(idx) {
		xx <- (avg.expr[, idx] / rowMeans(avg.expr[, -idx]) * genes.weight) %>% .[!is.na(.)]
		xx[order(xx) %>% rev][1 : count]
	}) %>% `names<-`(colnames(avg.expr %>% data.frame))
	return(seed.genes)
}

##########################################################################################################
#' @title buildSigMatrix
#' @param markers.lst A list of markers for cell types.
#' @param avg.expr Average gene expression of reference.
#' @param min.size Minimum size of each subset. Default: 10.
#' @param max.size Maximum size of each subset. Default: 200.
#' @return Signature matrix.
#' @export 
#' @example	

##########################################################################################################

buildSigMatrix <- function(markers.lst, avg.expr, min.size = 10, max.size = 200, subset.size = NULL) {
	if (!is.null(subset.size)) {
		sigmat.res <- lapply(markers.lst, function(genes) genes[1 : min(subset.size, length(genes))]) %>% unlist(.) %>% unique(.) %>% avg.expr[., ]
		return(sigmat.res)
	}
	max.size <- lapply(markers.lst, length) %>% unlist(., use.names = FALSE) %>% min(.) %>% min(., max.size)
	if (max.size < 10) {
		message('Warning: too few marker genes to be used accurately for deconvolution!')
		return(avg.expr[unlist(markers.lst) %>% unique, ])
	}	
	kappa.vals <- lapply(min.size : max.size, function(it.num) {		
		kp.val <- lapply(markers.lst,function(genes) genes[1 : it.num]) %>% 
			unlist(., use.names = FALSE) %>% 
			unique(.) %>% 
			avg.expr[., ] %>% 
			kappa(.) 
		kp.val 
	})
    group.size <- min.size + which.min(kappa.vals) - 1
	sigmat.res <- lapply(markers.lst, function(genes) genes[1 : group.size]) %>% unlist(.) %>% unique(.) %>% avg.expr[., ]
	return(sigmat.res)
}

##########################################################################################################
#' @title selectFeatures
#' @description Apply Lasso method to filter the features of logistic regression.
#' @param ref.obj Seurat object.
#' @param markers Cell type-specific genes.
#' @param sig.mat Signature matrix.
#' @param subset.size Marker size of each cell type. Default: NULL.
#' @param rand.seed Set random seed. Default: 1234.
#' @return signature gene list
#' @export

##########################################################################################################

selectFeatures<-function(ref.obj, markers, sig.mat, subset.size = NULL, rand.seed = 1234) {
	if (!is.null(subset.size)) {
		sig.genes <- lapply(markers, function(mm) {
			mm[1 : min(subset.size, length(mm))]
		}) %>% unlist %>% unique
		return(sig.genes)
	}
    set.seed(rand.seed)
	markers.uniq <- markers %>% unlist(.) %>% unique
	if(length(levels(ref.obj)) > 15) {
		return(rownames(sig.mat))
	} else {
		fitcv <- glmnet::cv.glmnet(
			ref.obj[markers.uniq, ] %>% GetAssayData(.) %>% as.matrix %>% t, 
			Idents(ref.obj) %>% as.vector,
			family = "multinomial", 
			alpha = 1,
			nfolds = 5
		)
	}
	coefs <- coef(fitcv, s = fitcv$lambda.min) %>% as.data.frame(.)
    sig.genes <- names(which(rowSums(coefs) != 0))
    return(sig.genes)
}

##########################################################################################################
#' @title downSamplSeurat
#' @param obj Seurat object
#' @param cnt Sample size for each ident, default: 2000.
#' @param seed Randon seed, default: 123.
#' @return Subset of seurat object.
#' @export

##########################################################################################################

downSamplSeurat <- function(obj, cluster.col = NULL, cnt = 200, seed = 123, percent = NULL) {
    set.seed(seed)
	if (!is.null(cluster.col)) Idents(obj) <- obj@meta.data[, cluster.col]
    cells <- Idents(obj) %>% table
    sub.cells <- sapply(names(cells), function(xx) {
        sub.cells <- Idents(obj)[Idents(obj) == xx] %>% names
        cnt <- ifelse(is.null(percent), cnt, length(sub.cells) * percent)
        if (length(sub.cells) > cnt) sub.cells <- sample(sub.cells, cnt, replace = FALSE)
        return(sub.cells)
    }) %>% unlist(use.names = F)
    subset(obj, cells = sub.cells)
}

