penffr1 <- function(Y.mat, t.mat = "all.ok", X_fd.list, X.scal = data.frame(), nbasis = "all.ok", n.order = 4, pen = T){

  '
  Y.mat       : matrix response with of n rows and m = max m_i columns with na values after the end of the observation for each sample
  t.mat       : a matrix of time with the corresponded values of Y.mat
  X_fd.list   : list of functional predictors
  X.scal      : data frame of the non functional predictors
  nbasis      : integer for the number of basis for functional parameter
  Pen         : boolean value for the penalization
  '

  # Check if the t.mat is not provided ----
  n <- nrow(Y.mat)
  m <- ncol(Y.mat)
  d <- length(X_fd.list)
  if (all(t.mat == "all.ok")) {
    obs.grid <- ((1:m)-1)/(m-1)
    t.mat <- matrix(rep(obs.grid, n), nrow = n, byrow = T)
  }

  # Check if nbasis is provided ----
  if (all(nbasis == "all.ok")){
    nbasis = min(20, floor((n*m)/(d+1)))
  }

  # Build the dataframe and penalty matrix for the model ----
  source("R/func-to-mat1.R")
  data <- get.data1(X_fd.list = X_fd.list,
                    Y.mat = Y.mat, nbasis = nbasis,
                    t.mat = t.mat, n.order = n.order)

  if (ncol(X.scal) > 0){
    # Formatting the non-functional data ----
    X.scal.new <- data.frame()
    for (i in 1:n) {
      ni <- length(data$time[data$id == i])
      tmp <- matrix(rep(as.vector(unlist(X.scal[i,])), each = ni),
                    nrow = ni, byrow = F)
      X.scal.new <- rbind(X.scal.new, tmp)
    }
    X.scal.new <- data.frame(X.scal.new)
    colnames(X.scal.new) <- paste("Z", 1:ncol(X.scal.new), sep = ".")
    data <- cbind.data.frame(data, X.scal.new)
  }
  #saveRDS(data, file = "data.rds")

  if (pen == F){

    ## Run the model ----
    my_formula <- as.formula(paste("output ~ 0", paste(colnames(data)[-(1:3)],
                                                       collapse = " + "), sep = " + "))
    model <- lm(my_formula, data = cbind.data.frame(data))

    ## get the functional parameters ----
    params <- model$coefficients
    params[is.na(params)] <- 0
    beta.fd <- lapply(1:(d+1), function(l){
      if (l != 1){
        fd(coef = params[(nbasis*(l-1)+1):(nbasis*l)],
           basisobj = create.bspline.basis(rangeval = range(t.mat, na.rm = T),
                                           nbasis = nbasis, norder = 4),
           fdnames = list(main = paste("beta_",(l-1), sep = ""),
                          xlab = "Time (millisec.)",
                          ylab = paste("beta_",(l-1),sep = "")))
      } else{
        fd(coef = params[(nbasis*(l-1)+1):(nbasis*l)],
           basisobj = create.bspline.basis(rangeval = range(t.mat, na.rm = T),
                                           nbasis = nbasis, norder = 4),
           fdnames = list(main = paste("beta_",(l-1),sep = ""),
                          xlab = "Time (millisec.)",
                          ylab = paste("beta_",(l-1),sep = "")))
      }
    })

    if (ncol(X.scal) > 0) {
      ## get the non functional parameters ----
      beta.scal <- model$coefficients[-c(1:(nbasis*(d+1)))]
      ## results ----
      return(list(model = model, beta.fd = beta.fd, beta.scal = beta.scal,
                  n.grid = m))
    } else {
      ## results ----
      return(list(model = model, beta.fd = beta.fd, n.grid = m))
    }

  } else {
    source("R/Pensim1.R")

    ## Build the penalty matrix ----
    vals <- my_penmat1(obs.grid = t.mat, d = length(X_fd.list)+1,
                       nbasis = nbasis, n.order = n.order, deg = 2)

    # Run the functional model ----
    n.lam <- min(10, floor(exp(log(100)/(d+1))))
    fofreg <- Pensim1(data = data, lams = seq(0.5, 5.0,length.out = n.lam),
                      vals = vals, d = d, d.z = ncol(X.scal),
                      t.mat = t.mat)

    model <- fofreg$model
    lams <- fofreg$lambda

    # ## Run the model ----
    # library(glmnet)
    # tmp.model <- cv.glmnet(x = as.matrix(data[,-c(1:3)]),
    #                        y = data$output,
    #                        alpha = 0,
    #                        intercept = FALSE)
    # #print(tmp.model$lambda.min)
    #
    # model <- glmnet(x = as.matrix(data[,-c(1:3)]),
    #                 y = data$output, intercept = FALSE,
    #                 alpha = 0, lambda = tmp.model$lambda.min)

    ## get the functional parameters ----
    params <- model$coefficients #as.vector(coef(model))[-1]
    params[is.na(params)] <- 0
    beta.fd <- lapply(1:(d+1), function(l){
      if (l != 1){
        fd(coef = params[(nbasis*(l-1)+1):(nbasis*l)],
           basisobj = create.bspline.basis(rangeval = range(t.mat, na.rm = T),
                                           nbasis = nbasis, norder = n.order),
           fdnames = list(main = paste("beta_",(l-1), sep = ""),
                          xlab = "Time (millisec.)",
                          ylab = paste("beta_",(l-1),sep = "")))
      } else{
        fd(coef = params[(nbasis*(l-1)+1):(nbasis*l)],
           basisobj = create.bspline.basis(rangeval = range(t.mat, na.rm = T),
                                           nbasis = nbasis, norder = n.order),
           fdnames = list(main = paste("beta_",(l-1),sep = ""),
                          xlab = "Time (millisec.)",
                          ylab = paste("beta_",(l-1),sep = "")))
      }
    })

    if (ncol(X.scal) > 0) {
      ## get the non functional parameters ----
      beta.scal <- params[-c(1:(nbasis*(d+1)))]
      ## results ----
      return(list(model = model, beta.fd = beta.fd,
                  beta.scal = beta.scal, lambda = lams, n.grid = m))
    } else {
      ## results ----
      return(list(model = model, beta.fd = beta.fd, lambda = lams,
                  n.grid = m))
    }
  }
}


pred.penffr1 <- function(model, newX_fd.list = list(), newX.scal = data.frame(), t.mat = "all.ok"){

  '
  t.mat must be the time matrix of predictions with of n rows and m columns
  newX_fd.list is the list of functional predictors
  newX.scal is the data frame of the non functional predictors
  model is the penffr model
  '

  # Check the providing of new functional predictors ----
  if (length(newX_fd.list) == 0) {
    if (nrow(newX.scal) == 0) {
      pred <- predict(model$model)
    } else {
      return(NA)
    }

  } else {

    # observation grid of predictions ----
    n <- ncol(newX_fd.list[[1]][["coefs"]])
    if (all(t.mat == "all.ok")) {
      obs.grid <- seq(newX_fd.list[[1]][["basis"]][["rangeval"]][1],
                      newX_fd.list[[1]][["basis"]][["rangeval"]][2],
                      length.out = model$n.grid)
      t.mat <- matrix(rep(obs.grid, n), nrow = n, byrow = T)
    }

    d <- length(newX_fd.list)

    # Build the dataframe and penalty matrix for the model ----
    source("R/func-to-mat1.R")
    nbasis <- length(as.vector(model$beta.fd[[1]][["coefs"]]))
    n.order <- nbasis - length(model$beta.fd[[1]][["basis"]][["params"]])
    data <- get.data1(X_fd.list = newX_fd.list,
                      Y.mat = t.mat, nbasis = nbasis,
                      t.mat = t.mat, n.order = n.order)

    if (ncol(newX.scal) > 0){
      # Formatting the non-functional data ----
      X.scal.new <- data.frame()
      for (i in 1:n) {
        ni <- length(data$time[data$id == i])
        tmp <- matrix(rep(as.vector(unlist(newX.scal[i,])), each = ni),
                      nrow = ni, byrow = F)
        X.scal.new <- rbind(X.scal.new, tmp)
      }
      X.scal.new <- data.frame(X.scal.new)
      colnames(X.scal.new) <- paste("Z", 1:ncol(X.scal.new), sep = ".")
      data <- cbind.data.frame(data, X.scal.new)
    }

    # prediction ----
    pred <- predict(model$model, newdata = data, type = "response")
  }

  # result ----
  return(pred)

}




