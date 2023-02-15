options(warn = -1)

##########################################################################################################
#' @title scAnno
#' @param obj.seu Seurat object, which need to be annotated.
#' @param ref.expr Reference gene expression profile. 
#' @param ref.anno Cell type information of reference profile, corresponding to the above `ref.expr`.
#' @param save.markers Specified the filename of makers need to be saved. 
#' @param cluster.col Column name of clusters to be annotated in meta.data slot of query Seurat object. Default: seurat_clusters.
#' @param factor.size Factor size for scaling the weight of gene expression. Default: 0.1.
#' @param method Specified method for selecting candidate marker genes, co.exp [co-expression with seed genes] or diff.exp [differential expressed genes by Seurat method]. Default: co.exp.
#' @param seed.num Number of seed genes of each cell type for recognizing candidate markers, only used when method = 'co.exp'. Default: 10.
#' @param subset.size Marker size of each cell type from reference profile. Default: NULL.  
#' @param redo.markers Re-search candidate markers or not. Default: FALSE.
#' @param gene.anno Gene annotation data.frame. Default: gene.anno.
#' @param permut.num  Number of permutations for estimating p-values of annotations. Default: 100.
#' @param show.plot Show annotated results or not. Default: TRUE.
#' @param verbose Show running messages or not. Default: TRUE.
#' @param tcga.data.u Bulk data of pan-cancer in TCGA. 
#' @return Annotated Seurat object.Default: tcga.data.u
#' @export 
#' @example	

##########################################################################################################

scAnno <- function(
	query        = obj.seu, 
	ref.expr     = ref.expr, 
	ref.anno     = ref.anno,
	save.markers = NULL,
	cluster.col  = 'seurat_clusters',
	factor.size  = 0.1,
	pvalue.cut   = 0.01,
	method       = c('co.exp', 'diff.exp'),
	seed.num     = 10,
	subset.size  = NULL,
	redo.markers = FALSE,
	gene.anno    = gene.anno,
	permut.num   = 100, 
	show.plot    = TRUE,
	verbose      = TRUE,
	tcga.data.u  = tcga.data.u
) {	
	println('[INFO] Checking the legality of parameters', verbose = verbose)
	checkPython() && checkInputParams(obj.seu, ref.expr, ref.anno, save.markers, cluster.col)
	ref.obj <- createSeuratObjAndQC(ref.expr, ref.anno)
	ref.obj <- ref.obj[which(rownames(ref.obj) %in% gene.anno[, 'gene_name']), ]	
	println(
		sprintf(
			'[INFO] %g cell types in reference, %g clusters in query objects', 
			length(levels(ref.obj)), 
			obj.seu@meta.data[, cluster.col] %>% unique %>% length
		), 
		verbose = verbose
	)
	select.method <- match.arg(method)
	avg.expr.ref <- AverageExpression(ref.obj, slot = "data")$RNA
	marker.file <- file.path(tempdir(), paste0(save.markers, '.rds'))
	if (redo.markers || (!file.exists(marker.file))) {
		println('[INFO] Searching candidate marker genes...', verbose = verbose)
		ref.markers <- switch(
			select.method,
			co.exp = {
				seed.genes <- identSeedGenes(avg.expr.ref, factor.size, count = seed.num)
				ref.markers_single <- searchMarkersByCorr(
					GetAssayData(ref.obj) %>% as.matrix, 
					seed.genes, 
					scale.data = TRUE,
					p.cut = pvalue.cut,
					seed.num = seed.num,
					tcga.data.u = NULL
				)
				ref.markers_bulk <- searchMarkersByCorr(
					GetAssayData(ref.obj) %>% as.matrix, 
					seed.genes, 
					scale.data = TRUE,
					p.cut = pvalue.cut,
					seed.num = seed.num,
					tcga.data.u = tcga.data.u
				)
                ref.markers <- list()
				for (i in 1:length(ref.markers_single)){
				ref.markers[[i]] <- intersect(ref.markers_single[[i]],ref.markers_bulk[[i]])
				ifelse (length(ref.markers[[i]]) < 20,ref.markers[[i]] <- ref.markers_single[[i]],ref.markers[[i]] <- ref.markers[[i]])
				}
				names(ref.markers) <- names(ref.markers_single)
				ref.markers = ref.markers
			},
			diff.exp = {
				ref.markers <- searchMarkersByASeurat(ref.obj)
			}
		)
		saveRDS(ref.markers, file = marker.file)
	} else {
		ref.markers <- readRDS(marker.file)
	}
	# Deconvolution for predicting.
	println('[INFO] Deconvolution by using RLM method', verbose = verbose)
	sig.mat <- buildSigMatrix(
		ref.markers, 
		avg.expr.ref, 
		min.size = 10, 
		max.size = 200,
		subset.size = subset.size
	)
	coef.vals <- rlmDewcon(
		sig.mat, 
		AverageExpression(obj.seu, group.by = cluster.col, slot = "data")$RNA
	)
	# Logistic regression for predicting.
	println('[INFO] Logistic regression for cell-type predictions, waiting...', verbose = verbose)
	sig.features <- selectFeatures(
		ref.obj, ref.markers, 
		sig.mat, 
		subset.size = subset.size
	)
	lr.pred.res <- lrPredModel(
		ref.obj[sig.features, ], 
		obj.seu, 
		down.sample = 200
	)
	lr.pred.score <- lr.pred.res$cellType
	colnames(lr.pred.score) <- paste("C", colnames(lr.pred.score) %>% as.numeric, sep = "")
	
	#Merge the scores of the two models, and assign annotations to clusters.
	println('[INFO] Merging the scores of both models, and assign annotations to clusters', verbose = verbose)
	anno.res <- mergeModelScore(
		coef.vals, 
		lr.pred.score, 
		ref.obj, 
		obj.seu, 
		cluster.col
	)
	anno.label <- apply(anno.res$score, 2, which.max) %>% rownames(anno.res$score)[.]
	anno.score <- apply(anno.res$score, 2, max) %>% as.vector  
	obj.seu <- AddMetaData(
		obj.seu, 
		metadata = obj.seu@meta.data[, cluster.col] %>% as.vector, 
		col.name = 'scAnno'
	)
	for(i in 0 : length(levels(obj.seu@meta.data[, cluster.col])) - 1) {
		obj.seu$scAnno[obj.seu@meta.data[, cluster.col] == i] <- anno.label[i + 1]
	}
	
	println('[INFO] Estimating p-values for annotations...', verbose = verbose)
	rand.scores <- getRandomScore(
		obj.seu, 
		lr.pred.res$model, 
		sig.mat, 
		lr.pred.res$features, 
		cluster.col, 
		sample.times = permut.num
	)
	pvals <- getPvalues(anno.res$score, rand.scores)
	
	for (i in 1:length(names(pvals))){
	   ifelse(pvals[i] < 0.01 ,anno.res$label[i] <- anno.label[i],anno.res$label[i] <- "UNKNOWN")
	}
	
	if (show.plot) {
		(DimPlot(obj.seu, group.by = cluster.col, label = TRUE) | DimPlot(obj.seu, group.by = 'scAnno', label = TRUE)) %>% print
	}
	println('[INFO] Finish!', verbose = verbose)
	return(list(query = obj.seu, pred.label = anno.res$label, pred.score = anno.res$score, pvals = pvals))
}

