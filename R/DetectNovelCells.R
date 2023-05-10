#' Detect novel cells.
#'
#' This function allows you to detect any novel cell type emerged over time.
#' @param CASiList The CASi list.
#' @param corr.cutoff The cutoff correlation value to determine novel cells. Default to 0.4. Cells with the correlation larger than 0.4 will not be considered novel cells.
#' @param p.cutoff The cutoff p value of t-tests to determine novel cells. Default to 0.01. Cells with the p value smaller than 0.01 will be considered novel cells.
#' @param DrawUMAP Draw UMAP plots per cell type-specific correlation and cell labels.
#' @export

DetectNovelCells <- function(CASiList, corr.cutoff=0.4, p.cutoff=0.01, DrawUMAP = TRUE){
  N <- dim(CASiList[[1]])[2]

  set.seed(42)

  CASiList_6 <- CASiList[[6]]

  df <- data.frame(table(CASiList[[2]]))
  types <- c(df$Var1)
  # temp <- CASiList_6
  for (i in 1:length(CASiList_6)) {
    idx = as.character(CASiList_6[i])
    CASiList_6[i] = match(idx,types)
  }


  # CASiList_6 <- gsub("C", "", CASiList_6)
  CASiList_6 <- as.numeric(CASiList_6)

  train_temp <-
    cbind(as.matrix(CASiList[[1]]),
          as.numeric(as.factor(CASiList[[2]])))

  train_temp <- train_temp[!duplicated(row.names(train_temp)),]

  train_temp <- train_temp[sample(nrow(train_temp)),]

  colnames(train_temp)[colnames(train_temp) == ""] <- "label"

  nClass <- nlevels(as.factor(train_temp[,N+1]))

  nCells <- ncol(train_temp)-1 # t0 cell number


  # combine t0 and t1
  t0_time <- rep(0,dim(train_temp)[1])
  t0_all_time <- cbind(train_temp, t0_time)

  t1_time <- rep(1,dim(CASiList[[3]])[1])
  t1_all_pred_time <- cbind(CASiList[[3]],CASiList_6, t1_time)

  t0t1_all_pred <- rbind(t0_all_time,t1_all_pred_time)

  colnames(t0t1_all_pred)[colnames(t0t1_all_pred) == "t0_time"] <- "time"

  gene_names <- rownames(t0t1_all_pred)

  t0t1_umap_results <- umap::umap(t0t1_all_pred[,1:N], n_components=20)

  n_t0 <- dim(t0_all_time)[1] # t0 cell number

  # t0 umap

  t0_umap <- data.frame(gene_names[1:n_t0],t0t1_umap_results$layout[1:n_t0,],train_temp[,'label'])

  # for (i in 1:nClass){
  #   nam <- paste("t0_C", i, sep = "")
  #   val <- t0_umap[t0_umap[,21] == i,]
  #   assign(nam, val)
  # }

  lst_t0_mean <- vector(mode='list', length=nClass)

  for (m in 1:nClass){
    # nam <- paste("t0_mean_C", i, sep = "")
    lst_t0_mean[[m]] <- colMeans((t0_umap[t0_umap[,22] == m,])[,2:21])
  }



  ## t1 umap

  nClass_2 <- nlevels(as.factor(CASiList_6))

  n1_t1 <- n_t0+1
  n2_t1 <- dim(t0t1_umap_results$layout)[1]

  t1_umap <- data.frame(gene_names[n1_t1:n2_t1],t0t1_umap_results$layout[n1_t1:n2_t1,],CASiList_6)

  ## t1 clusters
  lst_t1 <- vector(mode='list', length=nClass)


  for (j in unique(CASiList_6)){
    # nam <- paste("t1_C", i, sep = "")
    lst_t1[[j]] <- t1_umap[t1_umap[,22] == j,]
  }


  # find correlation
  for (k in unique(CASiList_6)){
    corr <- cor(t(lst_t1[[k]][,2:21]),lst_t0_mean[[k]],use="everything")
    lst_t1[[k]] <- cbind(lst_t1[[k]], corr)
  }

  lst_t1 <- lst_t1[lengths(lst_t1) != 0]

  ##############################################################
  lst_t1_2 <- do.call(rbind, lst_t1)
  lst_t1_ordered <- lst_t1_2[match(rownames(CASiList[[3]]), lst_t1_2[,1]), ]

  # combine labels & correlation

  t1_label_corr <- data.frame(CASiList_6,lst_t1_ordered[,23])
  row.names(t1_label_corr) <- row.names(lst_t1_ordered)
  colnames(t1_label_corr) <- c("labels","corr")

  lst_t1_label_novel <- vector(mode='list', length=nClass_2)
  t1_cluster_genes <- vector(mode='list', length=nClass_2)
  # clusters <- vector()

  for (i in 1:nClass_2){

    t1_cluster_genes[[i]] <- CASiList[[3]][row.names(CASiList[[3]]) %in% lst_t1[[i]][,1], ]

    if (is.null(dim(t1_cluster_genes[[i]]))){
      t1_cluster_genes[[i]] <- t(matrix(t1_cluster_genes[[i]]))
    }

    # table(row.names(t1_cluster_genes[[i]])==lst_t1[[i]][,1])
    t1_cluster_genes[[i]] <- data.frame(t1_cluster_genes[[i]],lst_t1[[i]][,23])

    # if (dim(t1_cluster_genes[[i]])[1]>20) {
    # clusters <- c(clusters, i)
    # } else
    #   clusters <- clusters
  }

  ###################### t-test

  temp <- subset(t1_label_corr, select = "labels")
  temp <- data.frame(temp)
  temp$names <- row.names(temp)
  temp$names <- gsub("-", ".", temp$names)
  temp$labels <- as.numeric(temp$labels)



  for (i in 1:nClass_2){

    # t1_cluster_genes[[i]] <- CASiList[[3]][row.names(CASiList[[3]]) %in% lst_t1[[i]][,1], ]
    #
    # # table(row.names(t1_cluster_genes)==lst_t1[[i]][,1])
    # t1_cluster_genes[[i]] <- data.frame(t1_cluster_genes[[i]],lst_t1[[i]][,23])

    if (dim(t1_cluster_genes[[i]])[1]<=20) {
      break }
    else if (dim(t1_cluster_genes[[i]])[1]>20){
      sim <-
        Seurat::CreateSeuratObject(
          counts = t(t1_cluster_genes[[i]][, 1:N]),
          project = "sim",
          min.cells = 0,
          min.features = 0
        )
      sim <- Seurat::NormalizeData(sim)
      sim <- Seurat::ScaleData(sim, features = colnames(t1_cluster_genes[[i]]))
      sim <- Seurat::RunPCA(sim, features = colnames(t1_cluster_genes[[i]]),npcs = 20)
      sim <- Seurat::FindNeighbors(sim, dims = 1:20)
      reso_vec <- seq(0.001, 1, by = 0.01)
      K = 1
      j = 1
      while (K != 2) {
        sim <- Seurat::FindClusters(sim, resolution = reso_vec[j], algorithm = 1)
        j = j + 1
        K = length(unique(Seurat::Idents(sim)))
      }

      recluster_res <- data.frame(sim$seurat_clusters)
      names(recluster_res) <- "cluster"
      recluster_res_0 <- subset(recluster_res, cluster==0)
      recluster_res_1 <- subset(recluster_res, cluster==1)

      t1_recluster_0 <- lst_t1[[i]][lst_t1[[i]][,1] %in% row.names(recluster_res_0),]
      t1_recluster_1 <- lst_t1[[i]][lst_t1[[i]][,1] %in% row.names(recluster_res_1),]

      x <- t1_recluster_0[,23]
      y <- t1_recluster_1[,23]
      t1_tt_res <- t.test(x, y, paired = FALSE, alternative = "two.sided", var.equal = FALSE)


      if (t1_tt_res$p.value <= p.cutoff & mean(abs(x))<corr.cutoff | mean(abs(y))<corr.cutoff) {

        if (abs(mean(x)) < abs(mean(y))) {
          recluster_res <- data.frame(recluster_res_0)
        } else {
          recluster_res <- data.frame(recluster_res_1)
        }

        recluster_res$names <- row.names(recluster_res)
        recluster_res$names <- gsub("-", ".", recluster_res$names)

        for (id in 1:nrow(recluster_res)) {
          temp$labels[temp$names %in% recluster_res$names[id]] <- 0
        }

        novel_labels <- as.vector(temp$labels)

        lst_t1_label_novel <- data.frame(t1_label_corr, novel_labels)
        # colnames(lst_t1_label_novel) <- c("labels","corr", "labels_novel")

      } else if (t1_tt_res$p.value > p.cutoff | mean(abs(x))>=corr.cutoff & mean(abs(y))>=corr.cutoff) {
        novel_labels <- as.vector(temp$labels)
        lst_t1_label_novel <- data.frame(t1_label_corr, novel_labels)
        # colnames(lst_t1_label_novel) <- c("labels","corr", "labels_novel")

      }

    } # table(lst_t1_label_novel[,1]==lst_t1_label_novel[,3])

  }

  lst_t1_label_novel_ordered <- lst_t1_label_novel[match(rownames(CASiList[[3]]), rownames(lst_t1_label_novel)), ]
  colnames(lst_t1_label_novel_ordered) <- c("labels","corr", "labels_novel")

  df <- data.frame(table(CASiList[[2]]))
  types <- c(df$Var1)
  for (i in 1:length(lst_t1_label_novel_ordered$labels)) {
    idx = as.numeric(lst_t1_label_novel_ordered$labels[i])
    lst_t1_label_novel_ordered$labels[i] = as.character(types[idx])
  }


  temp2 <- as.numeric(lst_t1_label_novel_ordered$labels_novel)

  for (i in 1:length(temp2)) {
    idx = temp2[i]
    if (idx != 0) {
      lst_t1_label_novel_ordered$labels_novel[i] = as.character(types[idx])
    }

    else if (idx ==0){
      lst_t1_label_novel_ordered$labels_novel[i] == "0"
    }
  }

  lst_t1_label_novel_ordered$labels_novel <- replace(lst_t1_label_novel_ordered$labels_novel, lst_t1_label_novel_ordered$labels_novel=="0", "Novel")


  CASiList[[7]] <- lst_t1_label_novel_ordered

  names(CASiList)[7] <- "Data: t1 cell types (potential novel cells)"

  if (DrawUMAP==TRUE){
    set.seed(42)
    t1_umap_results <- umap::umap(CASiList[[3]])
    ## Corr
    t1_umap_plot <- data.frame(t1_umap_results$layout,corr = abs(CASiList[[7]][,2]))

    p<-ggplot2::ggplot(t1_umap_plot,
                       ggplot2::aes(
                x = X1,
                y = X2 # label points with different colors for each `subgroup`
              )) +
      ggplot2::geom_point(size=0.6,ggplot2::aes(colour = corr)) +ggplot2::labs(color = 'Correlation') + ggplot2::theme_bw()

    p1 <- p+ ggplot2::scale_colour_gradient(low = "red", high = "grey", na.value = "red",limits=c(0, 1), breaks=seq(0,1,by=0.25)) + ggplot2::ggtitle("Correlation for checking any novel cell type") +ggplot2::theme(legend.text=ggplot2::element_text(size=11), plot.title = ggplot2::element_text(hjust = 0.5,size=16), legend.background = ggplot2::element_rect(fill = "white", color = "white"))+ggplot2::xlab("UMAP_1")+ggplot2::ylab("UMAP_2")

    print(p1)


    ## novel vs existed

    if (sum(stringr::str_detect(CASiList[[7]][,3], '^Novel$')) > 0) {

      # table(CASiList[[7]][,3])

      type <- CASiList[[7]][,3]
      type <- replace(type, type!="Novel", "Existed")

      t1_umap_plot <- data.frame(t1_umap_results$layout,cluster = as.factor(type))

      # lock in factor level order
      t1_umap_plot$cluster <- factor(t1_umap_plot$cluster, levels = unique(t1_umap_plot$cluster))

      p<-ggplot2::ggplot(t1_umap_plot,
                         ggplot2::aes(
                  x = X1,
                  y = X2 # label points with different colors for each `subgroup`
                )) +
        ggplot2::geom_point(size=0.6,ggplot2::aes(color = cluster)) +ggplot2::labs(color = 'Cell type') + ggplot2::theme_bw()

      p2<- p + ggplot2::ggtitle("Potential novel cells in t1 data") +ggplot2::theme(legend.text=ggplot2::element_text(size=11), plot.title = ggplot2::element_text(hjust = 0.5,size=16), legend.background = ggplot2::element_rect(fill = "white", color = "white"))+ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=1.5,shape=16))) +ggplot2::xlab("UMAP_1")+ggplot2::ylab("UMAP_2")

      print(p2)

    }

    else {
      print("No novel cells have been detected under the current setting of p.cutoff and corr.cutoff.")
    }


    ## Label
    t1_umap_plot <- data.frame(t1_umap_results$layout,cluster = as.factor(CASiList[[7]][,3]))

    # lock in factor level order
    t1_umap_plot$cluster <- factor(t1_umap_plot$cluster, levels = unique(t1_umap_plot$cluster))

    p<-ggplot2::ggplot(t1_umap_plot,
                       ggplot2::aes(
                x = X1,
                y = X2 # label points with different colors for each `subgroup`
              )) +
      ggplot2::geom_point(size=0.6,ggplot2::aes(color = cluster)) +ggplot2::labs(color = 'Cell type') + ggplot2::theme_bw()

    p3 <- p + ggplot2::ggtitle("Predicted cell types of t1 cells") +ggplot2::theme(legend.text=ggplot2::element_text(size=11), plot.title = ggplot2::element_text(hjust = 0.5,size=16), legend.background = ggplot2::element_rect(fill = "white", color = "white"))+ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=1.5,shape=16))) +ggplot2::xlab("UMAP_1")+ggplot2::ylab("UMAP_2")

    print(p3)

  }

  return(CASiList)
}
