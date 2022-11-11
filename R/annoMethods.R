##########################################################################################################
#' @title rlmDewcon 
#' @param ref Cell type-specific gene reference expression profiles,the reference profile is the mean value of gene expression in the cell type
#' @param bulk Gene expression profiling,bulk is the gene expression profile after clustering.
#' @param weight Assign the weight of each genes, default: NULL.
#' @return A data.frame of coefficients for cell types.
#' @return 

##########################################################################################################

rlmDewcon <- function(ref, bulk, weight = NULL){
  over.genes <- intersect(rownames(ref),rownames(bulk))
  bulk <- bulk[over.genes, ]
  ref <-ref[over.genes,]
  colnames(ref) <- gsub(' |&|\\+|\\.|-|/|\\(|\\)', '_',colnames(ref))
  coefs <- c()
  for(i in 1 : ncol(bulk)){
    data.test <- cbind.data.frame(ref, mix = bulk[, i])
    lmfit <- rlm(mix ~ . - 1, data = data.test)
    coef <- lmfit$coefficients
    coef[coef < 0] <- 0
    coefs <- rbind(coefs, coef / sum(coef))
  } 
   return (t(coefs) %>% data.frame)
}
  
##########################################################################################################
#' @title lrPredModel 
#' @param obj.ref Seurat object of reference gene expression profile.
#' @param obj.seu Seurat object of target gene expression profile.
#' @return Predicted results.
#' @export

##########################################################################################################
  
lrPredModel <- function(obj.ref, obj.seu) {
	python.path <- py_config()$python
	use_python(python.path)
	com.genes <- intersect(rownames(obj.ref), rownames(obj.seu))
	obj.ref.sub <- obj.ref[com.genes, ]
	obj.seu.sub <- obj.seu[com.genes, ]
	sklearn <- import('sklearn.linear_model')
	lr <- sklearn$LogisticRegression(max_iter = 100000, penalty = 'l2', C = 0.2)
	lr$fit(obj.ref.sub %>% GetAssayData(.) %>% as.matrix %>% t, Idents(obj.ref.sub) %>% as.vector)
    pred.prob <- lr$predict_proba(obj.seu.sub %>% GetAssayData(.) %>% as.matrix %>% t)
	pred.class <- lr$predict(obj.seu.sub %>% GetAssayData(.) %>% as.matrix %>% t)
	anno.cnt <- table(obj.seu.sub$seurat_clusters,pred.class)
	#anno.cnt <- table(Idents(obj.seu.sub), pred.class)
	anno.frac <- sweep(anno.cnt, 1, rowSums(anno.cnt), '/') %>% t
	return(list(model = lr, prob = pred.prob, cellType = anno.frac))
}  
