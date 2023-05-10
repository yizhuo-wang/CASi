#' Build LASiC object.
#'
#' This function allows you to build the LASiC list and all other functions will be applied to this list.
#' @param t0_count Gene expression matrix of the initial data point/t0.
#' @param t1_count Gene expression matrix of the following data points/t1.
#' @param t0_label Cell labels of the initial data point.
#' @param nfeature The number of top highly variable genes. Default to 2000.
#' @export


CreateCASiList <- function(t0_count, t1_count, t0_label=NULL, nfeature=2000) {


  train_count <- t0_count
  test_count <- t1_count

  # ### log normalization
  #
  # train_count <- log(train_count+1)

  ### sample normalization
  cmax <- apply(train_count, 2, max)
  cmin <- apply(train_count, 2, min)
  tmp1 <- sweep(train_count, 2, cmin, "-")
  tmp2 <- sweep(tmp1, 2, cmax-cmin, "/")

  ### let's try some feature selection here
  oosd <- apply(tmp2, 1, sd)
  tmp3 <- tmp2[order(oosd, decreasing = TRUE)[1:nfeature], ]
  train_count_out <- t(tmp3)
  rm(tmp1, tmp2, tmp3)

  ### sample normalization
  cmax <- apply(test_count, 2, max)
  cmin <- apply(test_count, 2, min)
  tmp1 <- sweep(test_count, 2, cmin, "-")
  tmp2 <- sweep(tmp1, 2, cmax-cmin, "/")

  ### let's try some feature selection here
  mmidx <- match(colnames(train_count_out), rownames(tmp2))
  tmp3 <- tmp2[mmidx, ]
  test_count_out <- t(tmp3)
  rm(tmp1, tmp2, tmp3)

  t0_count <- t(train_count_out)
  t1_count <- t(test_count_out)

  set.seed(42)


  if (!is.null(t0_label)) {

    temp <- rbind(as.matrix(t0_count),t0_label)

    temp <- temp[!duplicated(row.names(temp)), ]
    temp <- temp[,!duplicated(colnames(temp))]

    temp_idx <- dim(t0_count)[1]
    t0_label <- temp[temp_idx+1,]
    t0_count <- temp[1:temp_idx,]

    t1_count <- t1_count[!duplicated(row.names(t1_count)), ]
    t1_count <- t1_count[,!duplicated(colnames(t1_count))]

    CASiList <- list(t0_count = t(t0_count),
                      t0_label = as.factor(t0_label),
                      t1_count = t(t1_count))

    class(CASiList[[1]]) <- "numeric"

    names(CASiList)[1] <- "Data: t0 cells gene expression"
    names(CASiList)[2] <- "Data: t0 cell types"
    names(CASiList)[3] <- "Data: t1 cell gene expression"

  }

  else if (is.null(t0_label)) {

    t0_count <- t0_count[!duplicated(row.names(t0_count)), ]
    t0_count <- t0_count[,!duplicated(colnames(t0_count))]

    t1_count <- t1_count[!duplicated(row.names(t1_count)), ]
    t1_count <- t1_count[,!duplicated(colnames(t1_count))]

    t0_count <-
      Seurat::CreateSeuratObject(counts = t0_count,
                         min.cells = 3,
                         min.features = 200)

    # normalize data
    t0_count <-
      Seurat::NormalizeData(t0_count,
                    normalization.method = "LogNormalize",
                    scale.factor = 10000)

    t0_count <-
      Seurat::FindVariableFeatures(t0_count,
                           selection.method = "vst")
    # scale data
    all.genes <- rownames(t0_count)
    t0_count <- Seurat::ScaleData(t0_count, features = all.genes)
    # pca
    suppressMessages(t0_count <-
                       Seurat::RunPCA(t0_count, features = Seurat::VariableFeatures(object = t0_count)))

    t0_count <- Seurat::FindNeighbors(t0_count, dims = 1:10)
    t0_count <- Seurat::FindClusters(t0_count, resolution = 0.5)

    train_meta <- as.data.frame(t0_count@meta.data)
    train_label <- train_meta[, "seurat_clusters"]
    train_label <- paste("C", as.character(as.numeric(train_label)), sep="")

    train_count <- Seurat::GetAssayData(object = t0_count, slot = "counts")

    CASiList <- list(t0_count = t(train_count),
                      t0_label = as.factor(train_label),
                      t1_count = t(t1_count))

    names(CASiList)[1] <- "Data: t0 cells gene expression"
    names(CASiList)[2] <- "Data: t0 cell types"
    names(CASiList)[3] <- "Data: t1 cell gene expression"


  }
  return(CASiList)
}
