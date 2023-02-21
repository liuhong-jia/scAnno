##########################################################################################################
#' @title rlmDewcon 
#' @param ref Cell type-specific gene reference expression profiles,the reference profile is the mean value of gene expression in the cell type
#' @param bulk Gene expression profiling,bulk is the gene expression profile after clustering.
#' @return A data.frame of coefficients for cell types.
#' @return 

##########################################################################################################
options(warn = -1)

rlmDewcon <- function(ref, bulk){
  over.genes <- intersect(rownames(ref),rownames(bulk))
  bulk <- bulk[over.genes, ,drop = F] %>% `colnames<-`(paste0('C', colnames(bulk)))
  ref <- ref[over.genes,]
  colnames(ref) <- gsub(' |&|\\+|\\.|-|/|\\(|\\)', '_',colnames(ref))
  coefs <- c()
  for(i in 1 : ncol(bulk)){
    data.test <- cbind.data.frame(ref, mix = bulk[, i])
    lmfit <- MASS::rlm(mix ~ . - 1, data = data.test, maxit = 110)
    coef <- lmfit$coefficients
    coef[coef < 0] <- 0
    coefs <- rbind(coefs, coef / sum(coef))
  } 
   return (t(coefs) %>% data.frame %>% `colnames<-`(colnames(bulk)))
}
  
##########################################################################################################
#' @title lrPredModel 
#' @param obj.ref Seurat object of reference gene expression profile.
#' @param obj.seu Seurat object of target gene expression profile.
#' @return A list of predicted results.
#' @export

##########################################################################################################
  
lrPredModel <- function(obj.ref, obj.seu, down.sample = 100) {
	python.path <- py_config()$python
	use_python(python.path)
	if (down.sample != -1) {
		obj.ref <- downSamplSeurat(obj.ref, cnt = down.sample, seed = 123)
	}
	com.genes <- intersect(rownames(obj.ref), rownames(obj.seu))
	obj.ref.sub <- obj.ref[com.genes, ]
	obj.seu.sub <- obj.seu[com.genes, ]
	
	sklearn <- import('sklearn.linear_model')
	lr <- sklearn$LogisticRegression(max_iter = 100000, penalty = 'l2', C = 0.2)
	lr$fit(obj.ref.sub %>% GetAssayData(.) %>% as.matrix %>% t, Idents(obj.ref.sub) %>% as.vector)
	
    pred.prob <- lr$predict_proba(obj.seu.sub %>% GetAssayData(.) %>% as.matrix %>% t)
	pred.class <- lr$predict(obj.seu.sub %>% GetAssayData(.) %>% as.matrix %>% t)
	anno.cnt <- table(obj.seu.sub$seurat_clusters, pred.class)
	anno.frac <- sweep(anno.cnt, 1, rowSums(anno.cnt), '/') %>% t
	return(list(model = lr, prob = pred.prob, cellType = anno.frac, features = com.genes))
}

##########################################################################################################
#' @title mergeModelScore 
#' @return A list of predicted results.
#' @export

##########################################################################################################
 
mergeModelScore <- function(coef.vals, lr.pred.score, ref.obj, obj.seu, cluster.col) {
	rownames(coef.vals) <- gsub('_',' ', rownames(coef.vals))
	names <- cbind(rownames(coef.vals), names(table(Idents(ref.obj))))
	rownames(lr.pred.score) <- gsub("_|&|\\+|\\.|-|/|\\(|\\)", " ", lr.pred.score %>% as.matrix %>% rownames(.))
	com.names <- intersect(rownames(coef.vals), rownames(lr.pred.score))
	score <- -log10(1- coef.vals[com.names, ] * lr.pred.score[com.names, colnames(coef.vals)] + 10^(-10))
	rownames(score) <- names[match(rownames(score), names), 2]
	score[score < 0] <- 0
	prop <- sweep(score, 2, colSums(score), '/')
	prop[is.na(prop)] <- 0
	label <- apply(prop, 2, function(t) rownames(prop)[which.max(t)]) %>% as.vector
	names(label) <- paste("C", levels(obj.seu@meta.data[, cluster.col]), sep = "")
	return(list(score = prop, label = label))
}
