# B-Splines basis representation of multidimensionnal observation
my_functrep <- function(data, argvals = data.frame(), n.order = 4, nbasis = 0, name.var = "", name.obs = ""){

  '
  argvals is a matrix of n rows and m = max m_i columns with na values after the end of the observation for each sample
  data is a matrix with the corresponded values of argvals
  n.order is the order of the splines we will use for the representation (cubic splines by default)
  nbasis is the number of basis use
  name.var is the name of the variable
  name.obs is the name of any individuals
  '

  # Default parameters calibration ----
  if (nbasis == 0) {
    nbasis <- max(n.order + 2, ceiling(ncol(data)/10))
  }
  if (name.var == "") {
    name.var <- deparse(substitute(data))
  }
  if (all(name.obs == "")) {
    name.obs <- rownames(data)
  }

  # Functional expansion
  if (nrow(argvals) == 0) {
    # argvals <- data.frame(matrix(NA, nrow = nrow(data), ncol = ncol(data)))
    # for (u in 1:nrow(data)) {
    #   k <- length(na.omit(as.vector(unlist(data[u,]))))
    #   argvals[u,1:k] <- seq(0,1,length.out=k)
    # }
    basis <- create.bspline.basis(rangeval = c(0,1),
                                  norder = n.order,
                                  nbasis = nbasis)
    tmp.sm.cv <- Data2fd(argvals = seq(0,1, length.out=ncol(data)),
                         y = t(as.matrix(data)),
                         basisobj = basis)
  } else {

    # T = [a,b]
    grd <- range(argvals, na.rm = T)

    #
    data.Par <- lapply(1:nrow(data), function(u){
      #print(u)
      k <- length(na.omit(as.vector(unlist(argvals[u,]))))
      n.step <- max(2, ceiling(k/65))

      lim_a <- length(na.omit(as.vector(unlist(argvals[u,]))))
      lim_d <- length(na.omit(as.vector(unlist(data[u,]))))

      ## Set the basis
      breaks <- na.omit(as.vector(unlist(argvals[u,])))[(1:k)%%n.step == 0]
      if (grd[1] < breaks[1]){
        breaks <- c(grd[1], breaks)
      }
      if (breaks[length(breaks)] < grd[2]){
        breaks <- c(breaks, grd[2])
      }
      breaks <- sort(breaks)
      if (lim_d > lim_a){
        data[u,(1+lim_a):lim_d] <- rep(NA,lim_d-lim_a)
      }
      else if (lim_d < lim_a){
        argvals[u,(lim_d+1):lim_a] <- rep(NA,lim_a-lim_d)
      }
      basis <- create.bspline.basis(rangeval = range(argvals, na.rm = T),
                                    breaks = breaks,
                                    norder = n.order)

      ## Set the labels
      var.labels <- list(Time = "Time (sec.)",
                         obs = name.obs[u],
                         name.var = name.var)

      # search the optimum lambda by cross-validation
      lam <- 1e-10
      grid.lam <- sapply(1:40, function(v){
        lam <- lam*(1.301**v) # Increase the value of lambda
        s.fdPar <- fdPar(fdobj = basis, lambda = lam)
        s.Par <- smooth.basis(argvals = sort(as.vector(na.omit(as.vector(unlist(argvals[u,]))))),
                              y = as.vector(na.omit(as.vector(unlist(data[u,])))),
                              fdParobj = s.fdPar,
                              fdnames = var.labels,
                              method = "chol")
        gcv_save <- sum(s.Par$gcv)
        df_save <- s.Par$df

        list(params = c(gcv_save, df_save), smooth = s.Par)
      })
      lam_min <- which.min(sapply((1:40)[1:40%%2 == 1],
                                  function(v){grid.lam[[v]]})[1,])

      grid.lam[[(2*lam_min)]]
    })

    # Eval all the fd object at a unique grid
    tmp.eval <- lapply(1:length(data.Par), function(u){
      eval.fd(evalarg = seq(range(argvals, na.rm = T)[1], range(argvals, na.rm = T)[2],
                            length.out = 1000),
              fdobj = data.Par[[u]]$fd,
              returnMatrix = T)
    })

    # Build a n*m matrix for the fdata
    tmp.mdata <- matrix(0, ncol = 1000, nrow = nrow(data))
    for (i in 1:nrow(data)) {
      tmp.mdata[i,] <- tmp.eval[[i]]
    }

    # Set the Basis
    basis <- create.bspline.basis(rangeval = range(argvals, na.rm = T),
                                  nbasis = nbasis,
                                  norder = n.order)

    # Now we build the final fd object
    tmp.sm.cv <- smooth.basis(argvals = seq(range(argvals, na.rm = T)[1],
                                            range(argvals, na.rm = T)[2],
                                            length.out = 1000),
                              y = t(tmp.mdata),
                              fdParobj = basis,
                              fdnames = list(Time = "Day of Year",
                                             stations = name.obs,
                                             name.var = name.var),
                              method = "chol")$fd

  }

  return(tmp.sm.cv)
}
