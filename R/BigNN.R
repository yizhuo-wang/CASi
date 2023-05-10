BigNN <- function(train_count, train_label, nClass,
                  lossname = 'categorical_crossentropy',
                  last_act = "softmax",
                  verbose = TRUE) {

  ### Build a big neural network:
  model_big <- keras::keras_model_sequential() %>%
    keras::layer_dense(units = 256, activation = "relu", input_shape = ncol(train_count)) %>%
    keras::layer_dropout(0.6) %>%
    keras::layer_dense(units = 128, activation = "relu") %>%
    keras::layer_dropout(0.6) %>%
    keras::layer_dense(units = 64, activation = "relu") %>%
    keras::layer_dropout(0.6) %>%
    keras::layer_dense(units = nClass, activation = last_act) %>%
    keras::compile(
      loss = lossname,
      optimizer = "adam",
      metrics = c('accuracy')
    )
  history <- model_big %>%
    keras::fit(
      x = train_count,
      y = train_label,
      epochs = 100,
      batch_size = 256,
      validation_split =0.2,
      callbacks = list(
        keras::callback_early_stopping(monitor="accuracy",patience = 10,mode="max"),
        keras::callback_reduce_lr_on_plateau(monitor="accuracy",factor = 0.05,patience=5,mode="max")
      ),
      verbose = verbose,
      return_best_model=TRUE
    )

  return(list(model_big,history))
}
