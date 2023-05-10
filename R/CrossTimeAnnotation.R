#' Cross-timepoint cell annotation.
#'
#' This function allows you to annotate t1 cells.
#' @param CASiList The CASi list.
#' @param DrawUMAP Draw UMAP plots per time points and cell labels. Default to TRUE.
#' @importFrom dplyr %>%
#' @export

# BigNN <- function(train_count, train_label, nClass,
#                   lossname = 'categorical_crossentropy',
#                   last_act = "softmax",
#                   verbose = TRUE) {
#
#   ### Build a big neural network:
#   model_big <- keras::keras_model_sequential() %>%
#     keras::layer_dense(units = 256, activation = "relu", input_shape = ncol(train_count)) %>%
#     keras::layer_dropout(0.6) %>%
#     keras::layer_dense(units = 128, activation = "relu") %>%
#     keras::layer_dropout(0.6) %>%
#     keras::layer_dense(units = 64, activation = "relu") %>%
#     keras::layer_dropout(0.6) %>%
#     keras::layer_dense(units = nClass, activation = last_act) %>%
#     keras::compile(
#       loss = lossname,
#       optimizer = "adam",
#       metrics = c('accuracy')
#     )
#   history <- model_big %>%
#     keras::fit(
#       x = train_count,
#       y = train_label,
#       epochs = 100,
#       batch_size = 256,
#       validation_split =0.2,
#       callbacks = list(
#         keras::callback_early_stopping(monitor="accuracy",patience = 10,mode="max"),
#         keras::callback_reduce_lr_on_plateau(monitor="accuracy",factor = 0.05,patience=5,mode="max")
#       ),
#       verbose = verbose,
#       return_best_model=TRUE
#     )
#
#   return(list(model_big,history))
# }




CrossTimeAnnotation <- function(CASiList, DrawUMAP = TRUE){

  N <- dim(CASiList[[1]])[2]

  train_labels <- as.numeric(as.factor(CASiList[[2]]))

  train_temp <-
    cbind(as.matrix(CASiList[[1]]),
          train_labels)

  # train_temp <- train_temp[!duplicated(row.names(train_temp)),]
  #
  # train_temp <- train_temp[sample(nrow(train_temp)),]

  nClass <- nlevels(as.factor(train_temp[,N+1]))

  nCells <- ncol(train_temp)-1

  t0_count <- train_temp[, 1:nCells]
  t0_label <- keras::to_categorical(train_temp[, ncol(train_temp)] - 1)

  mod_nn <- suppressWarnings(BigNN(train_count = t0_count, train_label =t0_label, nClass = nClass))

  CASiList[[4]] <- mod_nn[[1]]
  CASiList[[5]] <- mod_nn[[2]]

  # predict t1 labels

  pred_t1 <- mod_nn[[1]] %>% predict(CASiList[[3]]) %>% keras::k_argmax()
  pred_t1 <- as.array(pred_t1)+1

  df <- data.frame(table(CASiList[[2]]))
  types <- c(df$Var1)
  for (i in 1:length(pred_t1)) {
    idx = as.numeric(pred_t1[i])
    pred_t1[i] = as.character(types[idx])
  }

  # pred_t1 <- paste("C", as.character(as.numeric(pred_t1)), sep="")

  CASiList[[6]] <- pred_t1


  if (DrawUMAP==TRUE){
    set.seed(42)
    train_temp <-
      data.frame(as.matrix(CASiList[[1]]),
                 label=as.factor(CASiList[[2]]))

    # combine t0 and t1
    t0_time <- rep(0,dim(train_temp)[1])
    t0_all_time <- data.frame(train_temp, time=t0_time)

    t1_time <- rep(1,dim(CASiList[[3]])[1])
    t1_all_pred_time <- data.frame(CASiList[[3]],label=as.factor(CASiList[[6]]), time=t1_time)

    t0t1_all_pred <- rbind(t0_all_time,t1_all_pred_time)

    t0t1_umap_results <- umap::umap(t0t1_all_pred[,1:N])


    ## UMAP time
    time_umap_plot <- data.frame(t0t1_umap_results$layout,time = as.factor(t0t1_all_pred$time))

    # lock in factor level order
    time_umap_plot$cluster <- factor(time_umap_plot$time, levels = unique(time_umap_plot$time))

    p<-ggplot2::ggplot(time_umap_plot,
                       ggplot2::aes(
                x = X1,
                y = X2 # label points with different colors for each `subgroup`
              )) +
      ggplot2::geom_point(size=0.6,ggplot2::aes(color = time),alpha = 5/10) +ggplot2::labs(color = 'Time point') + ggplot2::theme_bw()

    p1 <- p + ggplot2::ggtitle("UMAP: time point") +ggplot2::theme(legend.text=ggplot2::element_text(size=11), plot.title = ggplot2::element_text(hjust = 0.5,size=16), legend.background = ggplot2::element_rect(fill = "white", color = "white"))+ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=1.5,shape=16))) +ggplot2::xlab("UMAP_1")+ggplot2::ylab("UMAP_2")

    print(p1)
    # ggsave("MCL_all_time.png",width = 170, height = 130,units='mm')

    ## label UMAP

    type_umap_plot <- data.frame(t0t1_umap_results$layout,type = as.factor(t0t1_all_pred$label))

    # lock in factor level order
    type_umap_plot$type <- factor(type_umap_plot$type, levels = unique(type_umap_plot$type))

    p <-ggplot2::ggplot(type_umap_plot,
                        ggplot2::aes(
                 x = X1,
                 y = X2 # label points with different colors for each `subgroup`
               )) +
      ggplot2::geom_point(size=0.6,ggplot2::aes(color = type),alpha = 5/10) +ggplot2::labs(color = 'Cell type') + ggplot2::theme_bw()

    p2 <- p + ggplot2::ggtitle("UMAP: cell types") +ggplot2::theme(legend.text=ggplot2::element_text(size=11), plot.title = ggplot2::element_text(hjust = 0.5,size=16), legend.background = ggplot2::element_rect(fill = "white", color = "white"))+ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=1.5,shape=16))) +ggplot2::xlab("UMAP_1")+ggplot2::ylab("UMAP_2")

    print(p2)
    # ggsave("MCL_all_label.png",width = 170, height = 130,units='mm')

  }


  names(CASiList)[1] <- "Data: t0 cells gene expression"
  names(CASiList)[2] <- "Data: t0 cell types"
  names(CASiList)[3] <- "Data: t1 cell gene expression"

  names(CASiList)[4] <- "Neural network: setting"
  names(CASiList)[5] <- "Neural network: training history"
  names(CASiList)[6] <- "Data: t1 cell types"

  return(CASiList)

}
