##########################################################################################################
#' @title downSamplSeurat
#' @param obj Seurat object
#' @param seed set random seed, default: 123
#' @param percent Sampling ratio
#' @return Seurat sub object
#' @export 
#' @example	

##########################################################################################################

downSamplSeurat <- function(obj, seed = 123, percent = 0.5){
  set.seed(seed)
  cells <- Idents(obj) %>% table
  sub.cells <- sapply(names(cells), function(xx) {
    sub.cells <- Idents(obj)[Idents(obj) == xx] %>% names
    cnt <- ifelse(is.null(percent), cnt, length(sub.cells) * percent)
    if (length(sub.cells) > cnt) 
      sub.cells <- sample(sub.cells, cnt, replace = FALSE)
    return(sub.cells) }) %>% unlist(use.names = F)
  subset(obj, cells = sub.cells)
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

buildSigMatrix <- function(markers.lst, avg.expr, min.size = 10, max.size = 200) {
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


#' @title signature_gene
#' @description Apply Lasso method to filter the features of logistic regression
#' @param ref.obj Seurat object
#' @param markers Cell type-specific genes
#' @return signature gene list
#' @export
#' @example

signature_gene<-function(ref.obj,markers){
    set.seed(100)
	if(length(names(table(Idents(ref.obj))))>15){
	fitcv = glmnet(ref.obj[markers %>% unlist(.) %>% unique,] %>% GetAssayData(.) %>% as.matrix %>% t, Idents(ref.obj) %>% as.vector,family = "multinomial")
	}
	else{
	fitcv=cv.glmnet(ref.obj[markers %>% unlist(.) %>% unique,] %>% GetAssayData(.) %>% as.matrix %>% t, Idents(ref.obj) %>% as.vector,family="multinomial", alpha=1,nfolds=5)
	}
	coefs <- coef(fitcv,s = fitcv$lambda.min) %>% as.data.frame()
    signature.gene<-names(which(rowSums(coefs)!=0))
    return(signature.gene)
}


##########################################################################################################
#' @title searchMarkersByROC
#' @description Search marker genes of each cell type using ROC method.
#' @param obj Seurat object.
#' @return A list of markers for cell types.
#' @export

##########################################################################################################

searchMarkersByROC <- function(obj, ...) {
	options(future.globals.maxSize = 40000 * 1024^2)
	plan("multiprocess", workers = detectCores())
	roc.markers <- FindAllMarkers(obj, only.pos = TRUE, test.use = 'roc', ...)	
	return(split(roc.markers$gene, roc.markers$cluster))
}