#' Identify temporal differentially expressed genes (tEGs).
#'
#' This function allows you to find tDEGs.
#' @param CASiList The CASi list.
#' @param metadata The dataframe containing the Group variable and the Time variable.
#' @param pct.cutoff The percentage of cells that the gene is expressed. Default to 0.3.
#' @param p.cutoff The cutoff p value of t-tests to determine tDEGs. Default to 0.05. Genes with the p value smaller than 0.05 will be considered tDEGs.
#' @export

FindTemporalDEGs <- function(CASiList, metadata, pct.cutoff=0.3, p.cutoff=0.05){

  if("Time" %in% colnames(metadata)==FALSE)
  {
    print(cat("Please name the numerical time variable as 'Time'."))
  }
  else if ("Time" %in% colnames(metadata)==TRUE){

    if ("Group" %in% colnames(metadata)==FALSE){
      print(cat("Please name the binary group variable as 'Group'."))
    }
    else if ("Group" %in% colnames(metadata)==TRUE){
      N <- dim(CASiList[[1]])[2]

      counts_all <- rbind(as.matrix(CASiList[[1]]),as.matrix(CASiList[[3]]))

      counts_all <- log(counts_all+1)

      CASiList[[6]]<- as.character(CASiList[[6]])
      CASiList[[2]]<- as.character(CASiList[[2]])

      labels_all <- as.numeric(as.factor(c(CASiList[[2]],CASiList[[6]])))
      all_data <- data.frame(counts_all,label=labels_all,status=metadata$Group, time=metadata$Time)

      ## subset each cell type
      data_by_label <- split(all_data, all_data$label)
      nClass <- nlevels(as.factor(CASiList[[2]]))
      tDEGs <- vector(mode='list', length=nClass)

      for (i in 1:length(data_by_label)){
        idx1 <- dim(data_by_label[[i]])[1]

        genes.percent.expression <- colMeans(data_by_label[[i]][1:idx1,1:N]>0)

        genes.filter <- names(genes.percent.expression[genes.percent.expression>pct.cutoff])  #select genes expressed in at least 10% of cells

        if (length(genes.filter) == 0) {
          break }

        else if (length(genes.filter) != 0) {

          data_by_label_sub <- data_by_label[[i]][,genes.filter]


          meta <- data.frame(status=as.factor(data_by_label[[i]][,N+3]),time=as.factor(data_by_label[[i]][,N+2]))

          idx2 <- dim(data_by_label_sub)[2]

          expression <- Seurat::CreateSeuratObject(counts=t(data_by_label_sub))
          p_val <- matrix(NA,idx2)

          suppressWarnings(for (j in 1:idx2) {
            gene_name = Seurat::GetAssayData(object = expression,
                                     assay = "RNA",
                                     slot = "counts")[j,]

            meta$y <- gene_name

            meta <- na.omit(meta)
            tryCatch({
              glm_nb <-
                lm(
                  y ~ status + as.numeric(time) + as.numeric(time) * status,
                  data = meta
                )

              p_val[j] <- as.numeric(coef(summary(glm_nb))[, 4][4])}, error=function(e){})


          })

          names(p_val) <- colnames(data_by_label_sub[,1:idx2])
          p_val <- p_val[which(p_val < 0.05)]
          genes.filter.2 <- names(p_val)
          n2 <- length(names(p_val))

          if (n2!=0) {

            expression_filtered <- subset(x=expression, features = genes.filter.2)

            p_val_2 <- matrix(NA,n2)

            suppressWarnings(for (j in 1:n2) {
              gene_name = Seurat::GetAssayData(object = expression_filtered,
                                       assay = "RNA",
                                       slot = "counts")[j,]

              meta$y_2 <- gene_name

              tryCatch({
                glm_nb <-
                  glm.nb(
                    y_2 ~ status + as.numeric(time) + as.numeric(time) * status,
                    data = meta
                  )

                p_val_2[j] <- as.numeric(coef(summary(glm_nb))[, 4][4])}, error=function(e){})
            })

            names(p_val_2) <- genes.filter.2
            p_val_2 <- p_val_2[which(p_val_2 < p.cutoff)]

            tDEGs[[i]] <- p_val_2


          }
        }
      }

      CASiList[[8]] <- tDEGs

      names(CASiList)[8] <- "Cell-type-specific p-value of temporal DEGs"

      return(CASiList)

    }
  }
}
