##########################################################################################################
#' @title getRandomScore
#' @Description Obtain cell type background scores for random permutation tests
#' @param obj.seu Seurat object, which need to be annotated.
#' @param ref.obj Reference Seurat object.
#' @return Cell Type Score from Permutation Test.
#' @export 
#' @example	

##########################################################################################################

getRandomScore <- function(obj.seu, lr.model, sig.mat, sig.features, cluster.col, sample.times = 100){
	bulk.data <- GetAssayData(obj.seu, slot = 'data')
	pb <- progress::progress_bar$new(total = sample.times, clear = TRUE)
	rand.score <- lapply(1 : sample.times, function(idx) {
		pb$tick()
		set.seed(idx)
		rand.cells <- sample(colnames(bulk.data), min(floor(dim(obj.seu)[2] * 1 / 3), 1000))
		rand.expr <- bulk.data[, rand.cells] %>% as.matrix	
		coef.vals <- rlmDewcon(sig.mat, rand.expr %>% rowMeans %>% as.data.frame()) %>% `rownames<-`(colnames(sig.mat))
		clus.data <- obj.seu@meta.data[rand.cells, cluster.col]
		pred.class <- lr.model$predict(rand.expr[sig.features, ] %>% t)
		
		lr.prop <- rep(0, ncol(sig.mat)) %>% `names<-`(rownames(coef.vals))
		tmp.prop <- table(pred.class) / sum(table(pred.class))
		lr.prop[names(tmp.prop)] <- tmp.prop
		
		#Combine the scores of the two models
		score <- 0.436 * coef.vals[, 1] + 0.564 * lr.prop
		prop <- score / sum(score)
		prop[is.na(prop)] <- 0
		return(prop)
	}) %>% do.call(rbind, .) %>% as.matrix	
	return(rand.score)
}

##########################################################################################################
#' @title getPvalues
#' @Description Estimate p-values based on the random permutation test.
#' @param obj.seu Seurat object, which need to be annotated.
#' @param ref.obj Reference Seurat object.
#' @return Cell Type Score from Permutation Test.
#' @export 
#' @example	

##########################################################################################################

getPvalues <- function(obs.score, rand.score) {
	anno.label <- apply(obs.score, 2, which.max) %>% rownames(obs.score)[.]
	pvals <- lapply(1 : ncol(obs.score), function(idx) {
		cell.name <- anno.label[idx]
		p.val  <- pnorm(
			max(obs.score[, idx]), 
			mean = mean(rand.score[, cell.name]), 
			sd = sd(rand.score[, cell.name]) + 1e-5, 
			lower.tail = FALSE
		)
	}) %>% unlist %>% `names<-`(colnames(obs.score))
	return(pvals)
}
