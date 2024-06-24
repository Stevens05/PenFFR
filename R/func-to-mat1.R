# Categorial variable for taking the repeated observation of the response ----
my_funmat1 <- function(t.mat, Y.mat){

  time <- data <- ind <- c()
  for (u in 1:nrow(Y.mat)) {
    tmp.time <- as.vector(na.omit(as.vector(unlist(t.mat[u,]))))
    tmp.data <- as.vector(na.omit(as.vector(unlist(Y.mat[u,]))))
    if (length(tmp.data) != length(tmp.time)){
      k <- min(length(tmp.data), length(tmp.time))
      tmp.data <- tmp.data[1:k]
      tmp.time <- tmp.time[1:k]
    }
    ind <- c(ind, rep(u, length(tmp.data)))
    time <- c(time, tmp.time)
    data <- c(data, tmp.data)
  }

  return(cbind.data.frame(output = data,
                          id = as.factor(ind),
                          time = time))
}

# Functionnal representation as matrix : Functionnal covariates ----
my_funmat2 <- function(X_fd.list, basis.beta, t.mat, Y.mat){

  d <- length(X_fd.list)
  m <- ncol(t.mat)
  n <- ncol(X_fd.list[[1]][["coefs"]])

  # Build the response ----
  data <- my_funmat1(t.mat = t.mat, Y.mat = Y.mat)

  tmp.basis <- function(x){
    d <- length(X_fd.list)
    tmp <- lapply(1:d, function(l){
      tmp_phi <- eval.basis(evalarg = x, basisobj = X_fd.list[[l]][["basis"]])
      tmp_psi <- eval.basis(evalarg = x, basisobj = basis.beta)
      list(phi = tmp_phi, psi = tmp_psi)
    })
    return(tmp)
  }

  ## Basis functions of parameters
  basis_psi <- function(t){
    tmp_basis <- tmp.basis(t)

    b <- matrix(tmp_basis[[1]][["psi"]], nrow = 1)
    if (d>1){
      for (l in 2:d) {
        tmp <- tmp_basis[[l]][["psi"]]
        cpt <- length(tmp)
        b <- rbind(cbind(b, matrix(0, nrow = nrow(b), ncol = cpt)),
                   c(rep(0,ncol(b)), tmp))
      }
    }

    return(b)
  }

  ## Basis functions of covariates
  basis_phi <- function(t){
    tmp_basis <- tmp.basis(t)

    b <- matrix(tmp_basis[[1]][["phi"]], nrow = 1)
    if (d>1){
      for (l in 2:d) {
        tmp <- tmp_basis[[l]][["phi"]]
        cpt <- length(tmp)
        b <- rbind(cbind(b, matrix(0, nrow = nrow(b), ncol = cpt)),
                   c(rep(0,ncol(b)), tmp))
      }
    }
    return(b)
  }


  ## Coefficients of basis expansion of covariates
  coefs <- c()
  for (l in 1:d) {
    coefs <- rbind(coefs, X_fd.list[[l]][["coefs"]])
  }
  coefs <- rbind(1, coefs)

  ## Design matrix
  R.mat <- c()
  for (i in 1:n) {
    tmp.time <- data$time[data$id == i]
    tmp <- c()
    for (tj in tmp.time) {
      tmp.phi <- basis_phi(tj)
      tmp.psi <- basis_psi(tj)
      tmp.b <- eval.basis(evalarg = tj, basisobj = basis.beta)

      tmp.phi <- rbind(cbind(matrix(1), matrix(0, ncol = ncol(tmp.phi))),
                       cbind(matrix(0, ncol = 1, nrow = nrow(tmp.phi)), tmp.phi))
      tmp.psi <- rbind(cbind(tmp.b, matrix(0, ncol = ncol(tmp.psi))),
                       cbind(matrix(0, ncol = ncol(tmp.b), nrow = nrow(tmp.psi)),
                             tmp.psi))

      tmp <- rbind(tmp, t(t(tmp.psi) %*% tmp.phi %*% coefs[,i]))

    }
    R.mat <- rbind(R.mat, tmp)
  }

  data <- cbind.data.frame(data, data.frame(R.mat))
  colnames(data) <- c("output", "id", "time",
                      paste("X", 1:(ncol(data)-3), sep = "."))
  data$id <- as.factor(data$id)

  return(data)
}


get.data1 <- function(X_fd.list, Y.mat, t.mat, nbasis, n.order){

  basis.beta <- create.bspline.basis(rangeval = range(t.mat, na.rm = T),
                                     nbasis = nbasis, norder = n.order)
  data <- my_funmat2(X_fd.list = X_fd.list,
                     basis.beta = basis.beta,
                     Y.mat = Y.mat,
                     t.mat = t.mat)
  return(data)
}


my_penmat1 <- function(obs.grid, nbasis, n.order, d, deg=2){

  basis.beta <- create.bspline.basis(rangeval = range(obs.grid, na.rm = T),
                                     nbasis = nbasis, norder = n.order)
  t0 <- min(obs.grid, na.rm = T)
  tf <- max(obs.grid, na.rm = T)
  L.beta <- nbasis

  vals <- lapply(1:1, function(l){

    tmp.b1 <- diag(rep(1,L.beta))

    tmp.fd <- fd(coef = tmp.b1, basisobj = basis.beta)
    tmp.fd2 <- deriv.fd(tmp.fd, deg)

    sapply(1:L.beta, function(j){
      t(sapply(1:L.beta, function(i){
        f <- function(x){
          tmp.eval <- as.vector(t(eval.fd(evalarg = x, fdobj = tmp.fd2)))
          return(abs(tmp.eval[i]*tmp.eval[j]))
        }
        integrate(Vectorize(f), lower = t0, upper = tf,
                  subdivisions=2000)[["value"]]
      }))
    })
  })
  for (l in 2:d) {
    vals[[l]] <- vals[[1]]
  }

  return(vals)
}


