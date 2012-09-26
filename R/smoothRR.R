Ln <- function(beta, other) {
  ##  other <- list(y, x, delta, clsize, sigma)
  Y <- other[[1]]
  X <- other[[2]]
  delta <- other[[3]]
  clsize <- other[[4]]
  sigma <- other[[5]]
  weight <- other[[6]]
  Z <- other[[7]]
  p  <- ncol(X)
  N <- nrow(X)
  n <- length(clsize)
  ln <- double(1)
  .C("lfun", as.double(beta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize), as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(weight), as.double(Z), out=as.double(ln), PACKAGE = "aftgee") $out
}

abarlogfun <- function(beta, Y, X, delta, clsize, sigma, weight) {
    p <- ncol(X)
    N <- nrow(X)
    n <- length(clsize)
    a <- vector("double", p * p)
    matrix(.C("abarlogfun", as.double(beta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize), as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(weight), out = as.double(a), PACKAGE = "aftgee")$out, nrow = p)
}

abargehanfun <- function(beta, Y, X, delta, clsize, sigma, weight) {
    p <- ncol(X)
    N <- nrow(X)
    n <- length(clsize)
    a <- vector("double", p * p)
    matrix(.C("abargehanfun", as.double(beta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize), as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(weight), out = as.double(a), PACKAGE = "aftgee")$out, nrow = p)
}

omegaFun <- function(beta, Y, X, delta, clsize, weight) {
  p <- ncol(X)
  N <- nrow(X)
  n <- length(clsize)
  omega <- vector("double", p * p)
  matrix(.C("omegafun", as.double(beta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize), as.integer(n), as.integer(p), as.integer(N), as.double(weight), out = as.double(omega), PACKAGE = "aftgee")$out, nrow = p)
}

uiFun <- function(beta, Y, X, delta, clsize, sigma, n, Z, weight, smooth = TRUE, constant = 0) {
  N <- nrow(X)
  p <- ncol(X)
  ans <- numeric(p)
  sn <- vector("double", p)
  ## ans <- .C("unsfun", as.double(beta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize), as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(Z), as.double(weight), out = as.double(sn), PACKAGE = "aftgee")$out
  ans <- .C("ufun", as.double(beta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize), as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(Z), as.double(weight), out = as.double(sn), PACKAGE = "aftgee")$out
  ans <- ans - constant
}

uilogFun <- function(beta, Y, X, delta, clsize, sigma, n, Z, weight, smooth = TRUE, constant = 0, rankweight = "logrank") {
  N <- sum(clsize)
  p <- ncol(X)
  sn <- vector("double", p)
  ans <- numeric(p)
  pw <- rep(1, N)
  if (rankweight == "PW" | rankweight == "Prentice-Wilcoxon" | rankweight == "GP") {
    en <- Y - X %*% beta
    ik <- rep(1:N, each=N)
    jl <- rep(1:N, N)
    Shat <- NULL
    dummy <- 1:N
    ord <- order(en)
    ei <- en[ord]
    weighti <- weight[ord]
    deltai <- delta[ord]
    repeats <- table(ei)
    Shati <- survfit(Surv(ei, deltai) ~ 1, weights = weighti)$surv
    Shati <- rep(Shati, repeats)
    Shat[dummy[ord]] <- Shati
    pw <- Shat
  }
  if (rankweight == "GP") {
      pw <- pw ^ p
  }
  if (smooth == TRUE) {
      ans <- .C("ulogfun", as.double(beta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize), as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(Z), as.double(weight), as.double(pw), out = as.double(sn), PACKAGE = "aftgee")$out
  }
  if (smooth != TRUE) {
      ans <- .C("ulognsfun", as.double(beta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize), as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(Z), as.double(weight), as.double(pw), out = as.double(sn), PACKAGE = "aftgee")$out
  }
  ans - constant
}

viEmp <- function(beta, Y, delta, X, id, weight = rep(1, nrow(X)), nres = 500, mb = TRUE, zbeta = FALSE, smooth = TRUE, rankweight = "gehan"){
  p <- ncol(X)
  clsize <- unlist(lapply(split(id, id), length))
  n <- length(clsize)
  N <- sum(clsize)
  sigma <- diag(p)
  UnV <- zmat <- matrix(0, ncol = nres, nrow = p)
  for (i in 1:nres) {
    if ( mb == TRUE) {
      Z <- rep(rexp(length(clsize)), clsize)
    }
    if ( mb != TRUE) {
      Z <- rep(1, N)
    }
    if (zbeta == TRUE) {
      zb <- rnorm(p)
      newbeta <- beta + n ^ (-0.5) * zb
      zmat[,i] <- zb
    }
    if (zbeta != TRUE) {
      newbeta <- beta
    }
    sn <- vector("double", p)
    if (smooth == TRUE) {
        if (rankweight == "gehan") {
            UnV[,i] <- as.vector(.C("ufun", as.double(newbeta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize), as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(Z), as.double(weight), out = as.double(sn), PACKAGE = "aftgee")$out) / n
        }
        if (rankweight != "gehan") {
            UnV[,i] <- uilogFun(newbeta, Y, X, delta, clsize, sigma, n, Z, weight, smooth, constant = 0, rankweight) / n
            ## as.vector(.C("ulogfun", as.double(newbeta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize), as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(Z), as.double(weight), out = as.double(sn), PACKAGE = "aftgee")$out) / n
        }
    }
    if (smooth == FALSE) {
        if (rankweight == "gehan") {
            UnV[,i] <- as.vector(.C("unsfun", as.double(newbeta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize), as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(Z), as.double(weight), out = as.double(sn), PACKAGE = "aftgee")$out) / n
        }
        if (rankweight != "gehan") {
            UnV[,i] <- uilogFun(newbeta, Y, X, delta, clsize, sigma, n, Z, weight, smooth, constant = 0, rankweight) / n
            ## as.vector(.C("ulognsfun", as.double(newbeta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize), as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(Z), as.double(weight), out = as.double(sn), PACKAGE = "aftgee")$out) / n
        }
    }
}
  vi <- var(t(UnV))
  list(vi = vi, zmat = zmat, UnV = UnV)
}

get_si <- function(beta, Y, delta, X, id, weight = rep(1, nrow(X)), nres = 500, rankweight = "gehan", given = TRUE) {
    p <- ncol(X)
    clsize <- unlist(lapply(split(id, id), length))
    N <- sum(clsize)
    n <- length(clsize)
    en <- Y - X %*% beta
    ik <- rep(1:N, each=N)
    jl <- rep(1:N, N)
    Shat <- NULL
    dummy <- 1:N
    ord <- order(en)
    ei <- en[ord]
    weighti <- weight[ord]
    deltai <- delta[ord]
    repeats <- table(ei)
    Shati <- survfit(Surv(ei, deltai) ~ 1, weights = weighti)$surv
    Shati <- rep(Shati, repeats)
    Shat[dummy[ord]] <- Shati
    xdif <- X[ik,] - X[jl,]
    edif <- en[ik] - en[jl]
    ind <- ifelse(edif <= 0, 1, 0)
    minEn <- ifelse(edif <= 0, ik, jl)
    si <- s <- NULL
    if (rankweight == "gehan") {
        s <- weight[ik] * weight[jl] * delta[ik] * ind * xdif + weight[ik] * weight[jl] * log(Shat[minEn]) * xdif ## need second weight[ik]?
        s <- ifelse(s == Inf, 0, s)
        s <- ifelse(s == -Inf, 0, s)
        s <- ifelse(is.na(s) == TRUE, 0, s)
        ##  s <- s / N
        s <- rowsum(s, ik)  ## N X p
    }
    if (rankweight == "logrank") {
        haz <- -1 * log(Shat)
        haz <- ifelse(haz == Inf, 0, haz)
        haz <- ifelse(haz == -Inf, 0, haz)
        haz <- ifelse(is.na(haz) == TRUE, 0, haz)
        gamma1 <- rowsum(ind * X[jl,] * weight[ik] * weight[jl], ik)    ## weight[ik]
        gamma0 <- as.numeric(rowsum(ind * weight[ik] * weight[jl], ik)) ## weight[ik]
        si1i <-  si1 <- s2 <- matrix(0, nrow = N, ncol = p)
        si1 <- (X[1:N,] - gamma1 / gamma0)
        si1i <- si1[ord,]
        si1dif <- rbind(rep(0,p), diff(si1i))
        s2s <- apply(si1dif * haz[ord], 2, cumsum)
        s2[dummy[ord],] <- s2s
        s <- delta[1:N] * si1 - s2
    }
    s <- subset(s, given)
}

get_vi <- function(s, id, delta, weight, n, N, pn) {
    s1 <- rowsum(s, group = id)
    si1 <- lapply(split(s1, 1:n), function(x) x %o% x)
    s11 <- mapply("*", si1, weight, SIMPLIFY = FALSE)
    vi <- Reduce("+", s11) / N
    p <- ncol(s)
    v2 <- matrix(0, ncol = p, nrow = p)
    if (pn <1) {
        s2 <- rowsum(s * (1 - delta), group = id)
        si2 <- lapply(split(s2, 1:n), function(x) x %o% x)
        s22 <- mapply("*", si2, weight, SIMPLIFY = FALSE)
        v21n <- Reduce("+", s22) / N
        s3 <- s * weight * (1 - delta)
        v22n <- apply(s3, 2, sum) / N
        v22n <- v22n %o% v22n
        v2 <- ( (1 - pn) / pn) * (v21n - v22n)
        ## vi <- vi+ ( (1 - pn) / pn) * (v21n - v22n)
    }
    vi <- vi / n
    v2 <- v2 / n
    list(vi = vi, v2 = v2)
}

viClo <- function(beta, Y, delta, X, id, weight = rep(1, nrow(X)), nres = 500, p0 = 1, p1 = 1, rankweight = "gehan") {
    p <- ncol(X)
    clsize <- unlist(lapply(split(id, id), length))
    n <- length(clsize)
    N <- sum(clsize)
    s <- get_si(beta, Y, delta, X, id, weight, nres, rankweight)

    xi <- get_si(beta, Y, delta, X, id, weight, nres, rankweight, rep(delta[unique(id)] == 0, clsize))
    xiCase <- get_si(beta, Y, delta, X, id, weight, nres, rankweight, rep(delta[unique(id)] == 1, clsize))

    control <- which(delta[unique(id)] == 0)
    case <- which(delta[unique(id)] == 1)

    v1 <- get_vi(s, id, delta, weight, n, N, 1)$vi / N
    v21 <- get_vi(xi, id[control], delta[control], weight[control], length(control), N, p0)$v2 / length(control)
    v22 <- get_vi(xiCase, id[case], delta[case] - 1, weight[case], length(case), N, p1)$v2 / length(case)

    ## v21 <- get_vi(xi, id[control], delta[control], weight[control], length(control), length(control) + length(case), p0)$v2
    ## v22 <- get_vi(xiCase, id[case], delta[case], weight[case], length(case), length(control) + length(case), p1)$v2
    prDelta <- table(delta[unique(id)])
    ## vi <- v1 + v21 * p0 + v22 * p1
    ## vi <- v1 + v21 * (length(control) / N) + v22 * (length(case) / N)
    vi <- v1 + v21 * (prDelta[1] / sum(prDelta)) + v22 * (prDelta[2] / sum(prDelta))
    list(vi = as.matrix(vi), s = s)
}

isFun <- function(beta, Y, delta, X, id, weight = rep(1, nrow(X)), sigma, nres = 500, vClose = FALSE, p0 = 1, p1 = 1, rankweight = "gehan", omega = FALSE) {
  p <- ncol(X)
  clsize <- unlist(lapply(split(id, id), length))
  n <- length(clsize)
  N <- sum(clsize)
  UnMat <- zmat <- ahat <- NULL
  An.inv <- 1
  if (omega == TRUE && rankweight == "gehan") {
      vi <- omegaFun(beta, Y, X, delta, clsize, weight) * n## / n ^ 2
  }
  if (omega != TRUE && vClose == TRUE && rankweight == "gehan") {
      vi <- viClo(beta, Y, delta, X, id, weight, nres, p0, p1)$vi
  }
  if (omega != TRUE && vClose == TRUE && rankweight == "logrank") {
      vi <- viClo(beta, Y, delta, X, id, weight, nres, p0, p1, rankweight)$vi
  }
  if (omega != TRUE && vClose != TRUE) {
      vi <- viEmp(beta, Y, delta, X, id, weight, nres, mb = TRUE, zbeta = FALSE, rankweight = rankweight)$vi
  }
  if (rankweight == "gehan") {
      An <- abargehanfun(beta, Y, X, delta, clsize, sigma, weight) / n
  }
  if (rankweight == "logrank") {
      An <- abarlogfun(beta, Y, X, delta, clsize, sigma, weight) / n
  }
  if (qr(An)$rank != p) {
      covmat <- ginv(An) %*% vi %*% ginv(An)
      An.msg <- "An is singular"
      An.inv <- 0
  }
  if (qr(An)$rank == p) {
      covmat <- solve(An) %*% vi %*% solve(An)
      An.msg <- "An is nonsingular"
      An.inv <- 1
  }
  covmat <- matrix(as.numeric(covmat), p)
  list(covmat = covmat, vi = vi, An.msg = An.msg, An.inv = An.inv)
}

zlFun <- function(beta, Y, delta, X, id, weight = rep(1, nrow(X)), nres = 500, vClose = FALSE, p0 = 1, p1 = 1, rankweight = "gehan") {
  p <- ncol(X)
  clsize <- unlist(lapply(split(id, id), length))
  n <- length(clsize)
  N <- sum(clsize)
  UnMat <- zmat <- ahat <- unTime <- NULL
  An.inv <- 1
  sigma <- diag(p)
  UnV <- viEmp(beta, Y, delta, X, id, weight, nres, mb = FALSE, zbeta = TRUE, smooth = FALSE, rankweight = rankweight)
  zmat <- UnV$zmat
  UnV <- UnV$UnV
  if (vClose == TRUE) {
    vi <- viClo(beta, Y, delta, X, id, weight, nres, p0, p1, rankweight)$vi
  }
  if (vClose != TRUE) {
    vi <- viEmp(beta, Y, delta, X, id, weight, nres, mb = TRUE, zbeta = FALSE, smooth = FALSE, rankweight = rankweight)$vi
  }
  An <- matrix(0, nrow = p, ncol = p)
  for (i in 1:p) {
    An[i,] <- lm( UnV[i,] ~ matrix(t(zmat), ncol = p) - 1)$coef
  }
  if (qr(An)$rank != p) {
      covmat <- ginv(An) %*% vi %*% t(ginv(An))
      An.msg <- "An is singular"
      An.inv <- 0
  }
  if (qr(An)$rank == p) {
      covmat <- solve(An) %*% vi %*% solve(An)
      An.msg <- "An is nonsingular"
      An.inv <- 1
  }
  covmat <- covmat / n
  covmat <- matrix(as.numeric(covmat), p)
  list(covmat = covmat, vi = vi, An.msg = An.msg, An.inv = An.inv)
}

huangFun <- function(beta, Y, delta, X, id, weight = rep(1, nrow(X)), nres = 500, vClose = TRUE, p0 = 1, p1 = 1, rankweight = "gehan") {
  p <- ncol(X)
  clsize <- unlist(lapply(split(id, id), length))
  n <- length(clsize)
  N <- sum(clsize)
  betaTemp <- NULL
  sigma <- diag(p)
  UnMatV <- NULL
  if (vClose == TRUE) {
      vi <- viClo(beta, Y, delta, X, id, weight, nres, p0, p1, rankweight)$vi
  }
  if (vClose != TRUE) {
      vi <- viEmp(beta, Y, delta, X, id, weight, nres, mb = TRUE, zbeta = FALSE, smooth = FALSE, rankweight = rankweight)$vi
  }
  qq <- chol(vi)
  ## options(warn = -1)
  ## qq <- chol(vi, pivot = TRUE)
  ## options(warn = 0)
  ## if (max(abs(crossprod(qq) - vi)) > 1) {
  ##     stop("vi is not positive definite; suggest removing intercept.")
  ## }

  qq <- t(qq)
  newBeta <- NULL
  Z <- rep(1, N)
  newBeta <- matrix(0, ncol = p, nrow = p)
  for ( i in 1:p) {
      if (rankweight == "gehan") {
          bb <- BBsolve(beta, uiFun, Y = Y, X = X, delta = delta, clsize = clsize, sigma = sigma, n = n, Z = Z, weight = weight, smooth = TRUE, constant = qq[,i] * n ^ 0.5, quiet = TRUE)
          newBeta[, i] <- bb$par
      }
      if (rankweight == "logrank") {
          bb <- BBsolve(beta, uilogFun, Y = Y, X = X, delta = delta, clsize = clsize, sigma = sigma, n = n, Z = Z, weight = weight, smooth = TRUE, constant = qq[,i] * n ^ 0.5, rankweight = rankweight, quiet = TRUE)
          newBeta[, i] <- bb$par
      }
      if (rankweight == "PW" | rankweight == "Prentice-Wilcoxon") {
          bb <- BBsolve(beta, uilogFun, Y = Y, X = X, delta = delta, clsize = clsize, sigma = sigma, n = n, Z = Z, weight = weight, smooth = TRUE, constant = qq[,i] * n ^ 0.5, rankweight = rankweight, quiet = TRUE)
          newBeta[, i] <- bb$par
      }
  }
  dd <- newBeta - beta
  covmat <- n * t(dd) %*% dd
##  covmat <- covmat / n
  list(covmat = covmat)
}

smoothrr <- function(formula, data, subset, id,
                     Sigma = NULL,
                     weight = NULL,
                     binit = NULL,
                     nres = 100, p0 = 1, p1 = 1,
                     contrasts = NULL,
                     variance = "ISMB",
                     rankweight = "gehan",
                     control = aftgee.control()) {
  if (sum(!(variance %in% c("MB", "ZLCF", "ZLMB", "sHCF", "sHMB", "ISCF", "ISMB", "js"))) == 8) stop("Invaid variance estimates.")
  if (!(rankweight %in% c("gehan", "logrank", "PW", "Prentice-Wilcoxon", "GP"))) stop("Invaid rankweight weight.")
  scall <- match.call()
  mnames <- c("", "formula", "data", "weight", "subset", "na.nation", "id")
  cnames <- names(scall)
  cnames <- cnames[match(mnames, cnames, 0)]
  mcall <- scall[cnames]
  mcall[[1]] <- as.name("model.frame")
  m <- eval(mcall, parent.frame())
  id <- model.extract(m, id)
  if (is.null(id)) stop("id variable not found.")
  y <- model.extract(m, "response")
  N <- NROW(y)
  mterms <- attr(m, "terms")
  x <- model.matrix(mterms, m, contrasts) ## this x has interception
  weight <- model.extract(m, weight)
  if (is.null(weight)) weight <- rep(1, N)
  xnames <- colnames(x) ## [-1]
  if(is.null(Sigma)) {Sigma = diag(ncol(x))} ## x[,-1]
  out <- smoothrr.fit(Y = log(y[,1]), delta = y[,2], X = as.matrix(x), id = id, weight = weight, binit = binit, Sigma = Sigma, variance = variance, nres = nres, p0 = p0, p1 = p1, rankweight = rankweight, control = control) ## x[,-1]
  out$call <- scall
  out$vari.name <- xnames
  class(out) <- "smoothrr"
  return(out)
}


smoothrr.fit <- function(Y, delta, X, binit = NULL,
                         Sigma, id, weight = rep(1, nrow(X)),
                         variance = c("MB", "ZLCF", "ZLMB", "sHCF", "sHMB", "ISCF", "ISMB", "js"),
                         nres = 100, p0 = 1, p1 = 1, rankweight = "gehan",
                         control = aftgee.control()) {
  p <- ncol(X)
  go <- 1
  step <- 0
  if (is.null(binit)) {
    beta <- lm(Y ~ X - 1)$coef
  }
  beta.temp1 <- beta
  Sigma.temp1 <- Sigma
  beta.temp2 <- NULL
  Sigma.temp2 <- NULL
  clsize <- unlist(lapply(split(id, id), length))
  order <- unlist(lapply(clsize, function(x) 1:x))
  N <- sum(clsize)
  n <- length(clsize)
  Z <- rep(1, N)
  var.MB <- var.ZLCF <- var.ZLMB <- var.sHCF <- var.sHMB <- var.js <- var.ISCF <- var.ISMB <- NaN
  for (i in 1:control$maxiter) {
      ## Point Estimation
      if (rankweight %in% c("gehan", "PW", "Prentice-Wilcoxon", "GP") & i == 1) {
          tbeta <- system.time(beta.temp2 <- nlm(Ln, p = beta.temp1, other = list(Y, X, delta, clsize, Sigma.temp1, weight, Z), fscale = 0.01)$estimate)
          beta.conv <- 0
      }
      if (rankweight == "logrank") {
          tbeta <- system.time(temp <- BBsolve(beta.temp1, uilogFun, Y = Y, X = X, delta = delta, clsize = clsize, sigma = Sigma.temp1, n = n, Z = Z, weight = weight, smooth = TRUE, constant = 0, rankweight = rankweight, quiet = TRUE))
          beta.temp2 <- temp$par
          beta.conv <- temp$convergence
      }
      if (rankweight == "PW" | rankweight == "Prentice-Wilcoxon") {
          beta.temp1 <- beta.temp2
          tbeta <- system.time(temp <- BBsolve(beta.temp1, uilogFun, Y = Y, X = X, delta = delta, clsize = clsize, sigma = Sigma.temp1, n = n, Z = Z, weight = weight, smooth = TRUE, constant = 0, rankweight = rankweight, quiet = TRUE))
          beta.temp2 <- temp$par
          print(beta.temp1)

          ## if (abs(beta.temp1 - beta.temp2) < control$abstol) {break}
          if (max(abs(beta.temp1 - beta.temp2)) >= control$abstol) {
              beta.temp1 <- beta.temp2
              next
          }
          beta.conv <- ifelse(i == control$maxiter, 1, 0)
      }
      if (rankweight == "GP") {
          beta.temp1 <- beta.temp2
          tbeta <- system.time(temp <- BBsolve(beta.temp1, uilogFun, Y = Y, X = X, delta = delta, clsize = clsize, sigma = Sigma.temp1, n = n, Z = Z, weight = weight, smooth = TRUE, constant = 0, rankweight = rankweight, quiet = TRUE))
          beta.temp2 <- temp$par
          print(beta.temp1)

          ## if (abs(beta.temp1 - beta.temp2) < control$abstol) {break}
          if (max(abs(beta.temp1 - beta.temp2)) >= control$abstol) {
              beta.temp1 <- beta.temp2
              next
          }
          beta.conv <- ifelse(i == control$maxiter, 1, 0)
      }
      ## Ends Point Estimation

      pass.MB <- pass.ZLCF <- pass.ZLMB <- pass.sHCF <- pass.sHMB <- pass.ISCF <- pass.ISMB <- FALSE
      ZLMB.An.inv <- ZLCF.An.inv <- ISMB.An.inv <- ISCF.An.inv <- js.An.inv <- 1
      if (nres == 0) {break}
      ## begin multipler resampling method
      if (sum(variance %in% "MB") > 0 && pass.MB == FALSE) {
          temp <- matrix(0, nrow = nres, ncol = p)
          if(nres == 0) {
              var.MB <- NULL
              pass.MB <- TRUE
              if (sum(c("ZLCF", "js", "ZLMB", "sHCF", "sHMB", "ISCF", "ISMB") %in% variance) == 0) {break}
          }
          for (i in 1:nres) {
              Z <- rep(rexp(length(clsize)), clsize)

              if (rankweight == "gehan") {
                  temp[i,] <- nlm(Ln, p= beta.temp2, other = list(Y, X, delta, clsize, Sigma.temp1, weight, Z), fscale = 0.01)$estimate
              }
              if (rankweight == "logrank") {
                  temp[i,] <- BBsolve(beta.temp2, uilogFun, Y = Y, X = X, delta = delta, clsize = clsize, sigma = Sigma.temp1, n = n, Z = Z, weight = weight, smooth = TRUE, constant = 0, rankweight = rankweight, quiet = TRUE)$par
              }
              if (rankweight == "PW" | rankweight == "Prentice-Wilcoxon") {
                  temp[i,] <- BBsolve(beta.temp2, uilogFun, Y = Y, X = X, delta = delta, clsize = clsize, sigma = Sigma.temp1, n = n, Z = Z, weight = weight, smooth = TRUE, constant = 0, rankweight = rankweight, quiet = TRUE)$par
              }
              if (rankweight == "GP") {
                  temp[i,] <- BBsolve(beta.temp2, uilogFun, Y = Y, X = X, delta = delta, clsize = clsize, sigma = Sigma.temp1, n = n, Z = Z, weight = weight, smooth = TRUE, constant = 0, rankweight = rankweight, quiet = TRUE)$par
              }
          }
          var.MB <- var(temp)
          pass.MB <- TRUE
          if (sum(c("ZLCF", "js", "ZLMB", "sHCF", "sHMB", "ISCF", "ISMB") %in% variance) == 0) {break}
      } ## end MBipler resampling method

      ## begin Zeng and Lin's approach, close Si
      if (sum(variance %in% "ZLCF") > 0 && pass.ZLCF == FALSE) {
          var.ZLCF <- zlFun(beta.temp2, Y, delta, X, id, weight, nres, vClose = TRUE, p0, p1, rankweight)
          ## if (var.ZLCF$An.inv == 0) {print (var.ZLCF$An.msg)}
          ZLCF.An.inv <- var.ZLCF$An.inv
          var.ZLCF <- var.ZLCF$covmat
          pass.ZLCF <- TRUE
          if (sum(c("js", "ZLMB", "sHCF", "sHMB", "ISCF", "ISMB") %in% variance) == 0) {break}
      } ## end Zeng and Lin's approach, close Si

      ## begin Zeng and Lin's approach, emprical Si
      if (sum(variance %in% "ZLMB") > 0 && pass.ZLMB == FALSE) {
          var.ZLMB <- zlFun(beta.temp2, Y, delta, X, id, weight, nres, vClose = FALSE, p0, p1, rankweight)
          ## if(var.ZLMB$An.inv) {print(var.ZLMB$An.msg)}
          ZLMB.An.inv <- var.ZLMB$An.inv
          var.ZLMB <- var.ZLMB$covmat
          pass.ZLMB <- TRUE
          if (sum(c("js", "sHCF", "sHMB", "ISCF", "ISMB") %in% variance) == 0) {break}
      } ## end Zeng and Lin's approach, emprical Si

    ## begin Huang's approach, close Si
    if (sum(variance %in% "sHCF") > 0 && pass.sHCF == FALSE){
      var.sHCF <- huangFun(beta.temp2, Y, delta, X, id, weight, nres, vClose = TRUE, p0, p1, rankweight)$covmat
      pass.sHCF <- TRUE
      if (sum(c("js", "sHMB", "ISCF", "ISMB") %in% variance) == 0) {break}
    } ## end Huang's approach, close Si

    ## begin Huang's approach, emprical Si
    if (sum(variance %in% "sHMB") > 0 && pass.sHMB == FALSE){
      var.sHMB <- huangFun(beta.temp2, Y, delta, X, id, weight, nres, vClose = FALSE, p0, p1, rankweight)$covmat
      pass.sHMB <- TRUE
      if (sum(c("js", "ISCF", "ISMB") %in% variance) == 0) {break}
    } ## end Huang's approach, emprical Si

    ## begin IS's appraoch, close Si
    if (sum(variance %in% "ISCF") > 0 && pass.ISCF == FALSE) {
        var.ISCF <- isFun(beta.temp2, Y, delta, X, id, weight, Sigma.temp1, nres, vClose = TRUE, p0, p1, rankweight)
        ## if(var.ISCF$An.inv == 0){print(var.ISCF$An.msg)}
        ISCF.An.inv <- var.ISCF$An.inv
        var.ISCF <- var.ISCF$covmat
        if (sum(c("js", "ISMB") %in% variance) == 0) {break}
        pass.ISCF <- TRUE
    }
    ## end IS's approach, close Si

    ## begin IS's appraoch, emprical Si
    if (sum(variance %in% "ISMB") > 0 && pass.ISMB == FALSE) {
        var.ISMB <- isFun(beta.temp2, Y, delta, X, id, weight, Sigma.temp1, nres, vClose = FALSE, p0, p1, rankweight)
        ## if(var.ISMB$An.inv == 0) {print(var.ISMB$An.msg)}
        ISMB.An.inv <- var.ISMB$An.inv
        var.ISMB <- var.ISMB$covmat
        if (!("js" %in% variance)) {break}
        pass.ISMB <- TRUE
    }
    ## end IS's approach, emprical Si

    ## begin JS's iterative approach
    if (sum(variance %in% "js") > 0) {
      var.js <- isFun(beta.temp2, Y, delta, X, id, weight, Sigma.temp1, omega = TRUE, rankweight = rankweight)
      ## if(var.js$An.inv == 0) {print(var.js$An.msg)}
      js.An.inv <- var.js$An.inv
      Sigma.temp2 <- var.js$covmat
      e.sigma <- abs(Sigma.temp1-Sigma.temp2)
      e.beta <- abs(beta.temp1 - beta.temp2)
      e.rel <- max(max(e.sigma/abs(Sigma.temp2)), max(e.beta/abs(beta.temp2)))
      e.abs <- max(e.sigma, e.beta)
      if (abs(e.rel) < control$reltol) break
      if (abs(e.abs) < control$abstol) break
      ## mat add an absolute tolerance check
      beta.temp1 <- beta.temp2
      Sigma.temp1 <- Sigma.temp2
      ## var.js <- Sigma.temp2 / n
    } ## end js
  }
  covmat <- list(MB = var.MB, ZLCF = var.ZLCF, ZLMB = var.ZLMB, sHCF = var.sHCF, sHMB = var.sHMB, ISCF = var.ISCF, ISMB = var.ISMB, js = Sigma.temp2 / n)
  convergence <- if ( i == control$maxiter) 1 else 0
  out <- list(beta = beta.temp2, covmat = covmat, convergence = convergence, tbeta = tbeta, beta.conv = beta.conv, beta.conv.at = i, var.meth = variance,
              ZLMB.An.inv = ZLMB.An.inv, ZLCF.An.inv = ZLCF.An.inv, ISMB.An.inv = ISMB.An.inv, ISCF.An.inv = ISCF.An.inv, js.An.inv = js.An.inv)
  ## if (ZLMB.An.inv == 0 | ZLCF.An.inv == 0 | ISMB.An.inv == 0 | ISCF.An.inv == 0 | js.An.inv == 0) {
  ##      print("An is singular")
  ## }
  class(out) <- "smoothrr.fit"
  return(out)
}



