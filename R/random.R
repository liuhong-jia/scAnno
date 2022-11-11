
##########################################################################################################
#' @title get_random_score
#' @Description Obtain cell type background scores for random permutation tests
#' @param obj.seu Seurat object, which need to be annotated.
#' @param ref.obj Reference Seurat object.
#' @return Cell Type Score from Permutation Test.
#' @export 
#' @example	

##########################################################################################################




get_random_score<-function(obj.seu , ref.obj , ref.markers , sig.mat , max.names , second.names){
    
	    lr <- lrPredModel(ref.obj[ref.markers %>% unlist(.) %>% unique, ],obj.seu)$model
	    score <- lapply(1:100,function(idx){
			prop <- c()
			set.seed(idx)
			cluster <- obj.seu$seurat_clusters %>%as.vector
			obj.seu$seurat_clusters <- sample(cluster)
			coef.vals <- rlmDewcon(sig.mat, AverageExpression(obj.seu[rownames(sig.mat),],group.by="seurat_clusters")$RNA, weight = NULL)
			obj.seu.sub <- obj.seu[ref.markers %>% unlist(.) %>% unique, ]
			pred.class <- lr$predict(obj.seu.sub %>% GetAssayData(.) %>% as.matrix %>% t)
			anno.cnt <- table(obj.seu.sub$seurat_clusters,pred.class)
			lr.pred.score <- sweep(anno.cnt, 1, rowSums(anno.cnt), '/') %>% t
			colnames(lr.pred.score) <- colnames(lr.pred.score) %>% as.numeric +1
			colnames(lr.pred.score) <- paste("X",colnames(lr.pred.score),sep="")
			#Combine the scores of the two models
			rownames(coef.vals) <- gsub('_',' ' ,rownames(coef.vals))
			rownames(lr.pred.score) <- gsub("_|&|\\+|\\.|-|/|\\(|\\)"," ",lr.pred.score %>%as.matrix %>% rownames(.))
			com.names <- intersect(rownames(coef.vals),rownames(lr.pred.score))
			score = -log10(1- coef.vals[com.names,]%>%as.matrix * lr.pred.score[com.names,colnames(coef.vals)]+10^(-10))
			score[score < 0] <- 0
			prop <- rbind(prop, score/rowSums(score)) 
			prop[is.na(prop)] <- 0
		    
			diff <- vector() 
			for(i in 1:ncol(prop)){
				diff<- c(diff,prop[,i][max.names[i]]-prop[,i][second.names[i]]) %>% as.vector
			}
			diff
			})
		
			random_score <- do.call(rbind,score) %>% as.matrix
			#random_score[random_score <0 ]<-0
			
			cutoff <- random_score[nrow(random_score)*0.05,]
			return(list(random_score=random_score,cutoff=cutoff))
	}



 
 
 