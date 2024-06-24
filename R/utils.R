# Mean Square Error (MSE) ----
my.MSE <- function(betaact.vect, obs.grid, betahat.fd){
  b <- eval.fd(evalarg = obs.grid, fdobj = betahat.fd)
  tmp <- mean((betaact.vect - b)**2)
  return(tmp**(0.5))
}

# Functionnal R² ----
fR2 <- function(Y, Yhat){
  Ybar <- matrix(rep(colMeans(Y),n), nrow = n, byrow = T)
  tmp1 <- colSums((Y - Yhat)**2)
  tmp2 <- colSums((Y - Ybar)**2)

  return(1 - mean(tmp1/tmp2))
}

# Mean Square Prediction Error (MSPE) ----
my.MSPE <- function(Y.act, Y.hat, t.mat){
  Y.act <- as.vector(t(Y.act))
  Y.hat <- as.vector(t(Y.hat))
  len <- sapply(1:nrow(t.mat), function(u){
    length(as.vector(na.omit(as.vector(unlist(t.mat[u,])))))
  })
  tmp.1 <- sum((Y.act[1:len[1]] - Y.hat[1:len[1]])**2)/sum(Y.act[1:len[1]]**2)
  if (nrow(t.mat)>1) {
    tmp.all <- sapply(1:(nrow(t.mat)-1), function(u){
      tmp1 <- Y.act[(sum(len[1:u])+1):(sum(len[1:(u+1)]))]
      tmp2 <- Y.hat[(sum(len[1:u])+1):(sum(len[1:(u+1)]))]

      sum((tmp1 - tmp2)**2)/sum(tmp1**2)
    })

    return(c(tmp.1, tmp.all))
  } else {
    return(tmp.1)
  }
}

# Integrated Squared Error ----
my_ISE <- function(Y.act, Y.hat){

  # len <- sapply(1:nrow(t.mat), function(u){
  #   length(as.vector(na.omit(as.vector(unlist(t.mat[u,])))))
  # })
  #
  # if (nrow(t.mat)>1) {
  #   tmp.1 <- sum((Y.act[1:len[1]] - Y.hat[1:len[1]])**2)
  #   tmp.all <- sapply(1:(nrow(t.mat)-1), function(u){
  #     tmp1 <- Y.act[(sum(len[1:u])+1):(sum(len[1:(u+1)]))]
  #     tmp2 <- Y.hat[(sum(len[1:u])+1):(sum(len[1:(u+1)]))]
  #
  #     sum((tmp1 - tmp2)**2)
  #   })
  #
  #   return(c(tmp.1, tmp.all))
  #
  # } else {
  #   tmp.1 <- sum((Y.act[1:len] - Y.hat[1:len])**2)
  #   return(tmp.1)
  # }
  return(rowSums((Y.act - Y.hat)**2))
}

# Predict in mixture of experts ---
predict.moe <- function(params, params.gated, data, data.gated, m){

  if (ncol(params.gated) != 1) {
    K <- ncol(params)

    ## Prediction in all classes
    pred <- sapply(1:K, function(k){
      as.matrix(data[,-c(1:3)]) %*% params[,k]
    })

    ## multinomial logit model
    tmp <- sapply(1:K, function(k){
      exp(cbind(1, as.matrix(data.gated)) %*% as.vector(params.gated[,k]))
    })
    pik <- tmp/rowSums(tmp)

    # ## posterior probabilities
    # pik <- sapply(1:K, function(k){
    #   sapply(1:nrow(pred), function(i){
    #     pi_ik[i,k] * dnorm(x = pred[i,k], mean = pred[i,k], sd =  sqrt(sig[k]))
    #   })
    # })
    # pik <- pik/rowSums(pik)

    ## Deduce the belonging class
    z <- apply(pik, 1, function(x) which.max(x))

    clust <- data.frame(x = z)
    tab <- table(clust$x)
    if (length(tab) == 1) {
      #print("one class")
      Y_hat <- matrix(pred[,clust$x[1]], ncol = m, byrow = T)
    } else {
      clust$x <- as.factor(clust$x)
      pik <- data.frame(model.matrix(~x - 1, data = clust))
      #Y_hat <- matrix(rowSums(pik * pred), ncol = m, byrow = T)
      Y_hat <- rowSums(pik * pred)
    }

    return(list(pred = Y_hat, class = z))
  } else {
    ## Prediction in all classes
    Y_hat <- as.matrix(data[,-c(1:3)]) %*% c(params)
    ## Prediction
    #Y_hat <- matrix(pred, ncol = m, byrow = T)

    return(list(pred = Y_hat, class = rep(1, length(pred))))
  }
}

# Random effect matrix
Rdeff <- function(n,m){
  v <- matrix(1, nrow = m)
  tmp.mat <- v
  if (n>1) {
    for (i in 1:(n-1)) {
      tmp.mat <- rbind(cbind(tmp.mat, matrix(0, nrow = nrow(tmp.mat))),
                       cbind(matrix(0, nrow = m, ncol = ncol(tmp.mat)),v))
    }
  }
  return(tmp.mat)
}
