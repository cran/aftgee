Ln <- function(beta, other) {
  ##  other <- list(y, x, delta, clsize, sigma)
  Y <- other[[1]]
  X <- other[[2]]
  delta <- other[[3]]
  clsize <- other[[4]]
  sigma <- other[[5]]
  weights <- other[[6]]
  Z <- other[[7]]
  p  <- ncol(X)
  N <- nrow(X)
  n <- length(clsize)
  ln <- double(1)
  .C("lfun", as.double(beta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize), as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(weights), as.double(Z), out=as.double(ln), PACKAGE = "aftgee") $out
}

abarlogfun <- function(beta, Y, X, delta, clsize, sigma, weights) {
    p <- ncol(X)
    N <- nrow(X)
    n <- length(clsize)
    a <- vector("double", p * p)
    matrix(.C("abarlogfun", as.double(beta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize), as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(weights), out = as.double(a), PACKAGE = "aftgee")$out, nrow = p)
}

abargehanfun <- function(beta, Y, X, delta, clsize, sigma, weights) {
    p <- ncol(X)
    N <- nrow(X)
    n <- length(clsize)
    a <- vector("double", p * p)
    matrix(.C("abargehanfun", as.double(beta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize), as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(weights), out = as.double(a), PACKAGE = "aftgee")$out, nrow = p)
}

omegaFun <- function(beta, Y, X, delta, clsize, weights) {
  p <- ncol(X)
  N <- nrow(X)
  n <- length(clsize)
  omega <- vector("double", p * p)
  matrix(.C("omegafun", as.double(beta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize), as.integer(n), as.integer(p), as.integer(N), as.double(weights), out = as.double(omega), PACKAGE = "aftgee")$out, nrow = p)
}

uiFun <- function(beta, Y, X, delta, clsize, sigma, n, Z, weights, smooth = TRUE, constant = 0) {
  N <- nrow(X)
  p <- ncol(X)
  ans <- numeric(p)
  sn <- vector("double", p)
  ## ans <- .C("unsfun", as.double(beta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize), as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(Z), as.double(weights), out = as.double(sn), PACKAGE = "aftgee")$out
  ans <- .C("ufun", as.double(beta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize), as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(Z), as.double(weights), out = as.double(sn), PACKAGE = "aftgee")$out
  ans <- ans - constant
}

uilogFun <- function(beta, Y, X, delta, clsize, sigma, n, Z, weights, smooth = TRUE, constant = 0, rankweights = "logrank") {
  N <- sum(clsize)
  p <- ncol(X)
  sn <- vector("double", p)
  ans <- numeric(p)
  pw <- rep(1, N)
  if (rankweights == "PW" | rankweights == "Prentice-Wilcoxon" | rankweights == "GP") {
    en <- Y - X %*% beta
    ik <- rep(1:N, each=N)
    jl <- rep(1:N, N)
    Shat <- NULL
    dummy <- 1:N
    ord <- order(en)
    ei <- en[ord]
    weightsi <- weights[ord]
    deltai <- delta[ord]
    repeats <- table(ei)
    Shati <- survfit(Surv(ei, deltai) ~ 1, weights = weightsi)$surv
    Shati <- rep(Shati, repeats)
    Shat[dummy[ord]] <- Shati
    pw <- Shat
  }
  if (rankweights == "GP") {
      pw <- pw ^ p
  }
  if (smooth == TRUE) {
      ans <- .C("ulogfun", as.double(beta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize), as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(Z), as.double(weights), as.double(pw), out = as.double(sn), PACKAGE = "aftgee")$out
  }
  if (smooth != TRUE) {
      ans <- .C("ulognsfun", as.double(beta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize), as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(Z), as.double(weights), as.double(pw), out = as.double(sn), PACKAGE = "aftgee")$out
  }
  ans - constant
}

viEmp <- function(beta, Y, delta, X, id, weights = rep(1, nrow(X)), B = 500, mb = TRUE, zbeta = FALSE, smooth = TRUE, rankweights = "gehan", sigma = diag(ncol(X))){
  p <- ncol(X)
  clsize <- unlist(lapply(split(id, id), length))
  ## n <- length(clsize)
  ## N <- sum(clsize)
  n <- sum(unlist(lapply(split(weights, id), unique)))
  N <- sum(weights)
  UnV <- zmat <- matrix(0, ncol = B, nrow = p)
  for (i in 1:B) {
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
        if (rankweights == "gehan") {
            UnV[,i] <- as.vector(.C("ufun", as.double(newbeta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize), as.double(sigma), as.integer(length(clsize)), as.integer(p), as.integer(sum(clsize)), as.double(Z), as.double(weights), out = as.double(sn), PACKAGE = "aftgee")$out) # / n
        }
        if (rankweights != "gehan") {
            UnV[,i] <- uilogFun(newbeta, Y, X, delta, clsize, sigma, length(clsize), Z, weights, smooth, constant = 0, rankweights) # / n
            ## as.vector(.C("ulogfun", as.double(newbeta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize), as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(Z), as.double(weights), out = as.double(sn), PACKAGE = "aftgee")$out) / n
        }
    }
    if (smooth == FALSE) {
        if (rankweights == "gehan") {
            UnV[,i] <- as.vector(.C("unsfun", as.double(newbeta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize), as.double(sigma), as.integer(length(clsize)), as.integer(p), as.integer(sum(clsize)), as.double(Z), as.double(weights), out = as.double(sn), PACKAGE = "aftgee")$out) # / n ## n and N?
        }
        if (rankweights != "gehan") {
            UnV[,i] <- uilogFun(newbeta, Y, X, delta, clsize, sigma, length(clsize), Z, weights, smooth, constant = 0, rankweights) # / n
            ## as.vector(.C("ulognsfun", as.double(newbeta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize), as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(Z), as.double(weights), out = as.double(sn), PACKAGE = "aftgee")$out) / n
        }
    }
}
  vi <- var(t(UnV))
  list(vi = vi, zmat = zmat, UnV = UnV)
}

get_si <- function(beta, Y, delta, X, id, weights = rep(1, nrow(X)), B = 500, rankweights = "gehan") {
    p <- ncol(X)
    clsize <- unlist(lapply(split(id, id), length))
    N <- sum(clsize)
    en <- Y - X %*% beta
    ik <- rep(1:N, each=N)
    jl <- rep(1:N, N)
    Shat <- NULL
    dummy <- 1:N
    ord <- order(en)
    ei <- en[ord]
    weightsi <- weights[ord]
    deltai <- delta[ord]
    repeats <- table(ei)
    Shati <- survfit(Surv(ei, deltai) ~ 1, weights = weightsi)$surv
    Shati <- rep(Shati, repeats)
    Shat[dummy[ord]] <- Shati
    ## Shati <- as.matrix(survfit(Surv(ei, deltai) ~ 1, weights = weightsi, etype = withinIdi)$surv)
    ## Shati <- Shati[rep(1:nrow(Shati), repeats), ]
    ## Shat <- Shati[order(dummy[ord]),]
    xdif <- X[ik,] - X[jl,]
    edif <- en[ik] - en[jl]
    ind <- ifelse(edif <= 0, 1, 0)
    minEn <- ifelse(edif <= 0, ik, jl)
    si <- s <- NULL
    if (rankweights == "gehan") {
        ## s <- weights[ik] * weights[jl] * delta[ik] * ind * xdif + weights[ik] * weights[jl] * log(Shat[minEn]) * xdif
        s <- weights[jl] * delta[ik] * ind * xdif + weights[jl] * log(Shat[minEn]) * xdif
        if (length(which(s == Inf)) > 0) {
            s <- ifelse(s == Inf, 0, s)
        }
        if (length(which(s == -Inf)) > 0) {
            s <- ifelse(s == -Inf, 0, s)
        }
        if (sum(is.na(s) > 0)) {
            s <- ifelse(is.na(s) == TRUE, 0, s)
        }
        s <- rowsum(s, ik)  ## N X p
        }
    if (rankweights == "logrank") {
        haz <- -1 * log(Shat)
        haz <- ifelse(haz == Inf, 0, haz)
        haz <- ifelse(haz == -Inf, 0, haz)
        haz <- ifelse(is.na(haz) == TRUE, 0, haz)
        ## gamma1 <- rowsum(ind * X[jl,] * weights[ik] * weights[jl], ik)
        ## gamma0 <- as.numeric(rowsum(ind * weights[ik] * weights[jl], ik))
        gamma1 <- rowsum(ind * X[jl,] * weights[jl], ik)
        gamma0 <- as.numeric(rowsum(ind * weights[jl], ik))
        si1i <-  si1 <- s2 <- matrix(0, nrow = N, ncol = p)
        si1 <- (X[1:N,] - gamma1 / gamma0)
        si1i <- si1[ord,]
        si1dif <- rbind(rep(0,p), diff(si1i))
        s2s <- apply(si1dif * haz[ord], 2, cumsum)
        s2[dummy[ord],] <- s2s
        s <- delta[1:N] * si1 - s2
    }
    ## s <- subset(s, given)
    s
}

get_vi <- function(s, id, delta, weights, n) {
    ## clsize <- unlist(lapply(split(id, id), length))
    ## nstra <- length(unique(weights))
    clweights <- as.numeric(unlist(lapply(split(weights, id), unique)))
    s1 <- rowsum(s, group = id) ## / 2 ##
    si1 <- lapply(split(s1, 1:nrow(s1)), function(x) x %o% x)
    s11 <- mapply("*", si1, clweights, SIMPLIFY = FALSE)
    v1 <- Reduce("+", s11) / n
    v2 <- apply(s1 * clweights , 2, sum) / n
    v2 <- v2 %o% v2
    if (length(unique(weights)) == 1) {
        p <-  unique(weights) - 1
        vi <- p * (v1 - v2)
    }
    if (length(unique(weights)) > 1) {
        cweights <- unique(weights)
        cweights <- cweights[cweights != 1 & cweights != 0]
        vi <- (cweights - 1)* (v1 - v2)
    }
    ## v1 <- v1 / n
    ## v2 <- v2 / n
    ## vi <- vi / n
    list(vi = vi, v1 = v1, v2 = v2)
}



viClo <- function(beta, Y, delta, X, id, weights = rep(1, nrow(X)), B = 500, rankweights = "gehan", stratify = TRUE) {
    s <- v1 <- vi <- v2i <- NULL
    ## n <- length(unlist(lapply(split(weights, id), unique)))
    n <- sum(unlist(lapply(split(weights, id), unique)))
    p <- ncol(X)
    clsize <- unlist(lapply(split(id, id), length))
    clid <- unlist(lapply(clsize, function(x) 1:x))
    stra <- match(weights, unique(weights))
    dim <- unique(clsize)
    s <- get_si(beta, Y, delta, X, id, weights, B, rankweights)
    if (stratify == TRUE) {
        v1 <- get_vi(s, id, delta, weights, n)$v1
        vi <- v1 ## / n
        v2 <- matrix(0, ncol = p, nrow = p)
        if (length(unique(stra)) > 1) {
            for (i in 1:length(unique(stra))) {
                ns <- sum(unlist(lapply(split(weights[stra == i], id[stra == i]), unique)))
                ## weights at cluster level
                v2i <- get_vi(s = s[stra == i, ], id = id[stra == i], delta = delta[stra == i],
                              weights = weights[stra == i],
                              n = ns
                              )$vi
                ## strPr <- sum(unlist(lapply(split(weights[stra == i], id[stra == i]), unique))) / N
                strPr <- ns / n
                v2 <- v2 + v2i * strPr
            }
            ## vi <- (v1 + v2) / n
            vi <- v1 + v2
        }
    }
    if (stratify != TRUE) {
        v1 <- get_vi(s, id, delta, weights, n)$v1 # / sum(clsize)
        vi <- v1
        if (length(unique(stra)) > 1) {
            v2 <- get_vi(s, id, delta, weights * (1 - delta), n)$vi # / sum(clsize)
            vi <- v1 + v2
        }
    }
    list(vi = as.matrix(vi), s = s)
}


isFun <- function(beta, Y, delta, X, id, weights = rep(1, nrow(X)), sigma, B = 500, vClose = FALSE, rankweights = "gehan", omega = FALSE, stratify = TRUE) {
  p <- ncol(X)
  clsize <- unlist(lapply(split(id, id), length))
  n <- sum(unlist(lapply(split(weights, id), unique)))
  ## n <- length(clsize)
  UnMat <- zmat <- ahat <- NULL
  An.inv <- 1
  if (omega == TRUE && rankweights == "gehan") {
      vi <- omegaFun(beta, Y, X, delta, clsize, weights) * n## / n ^ 2
  }
  if (omega != TRUE && vClose == TRUE && rankweights == "gehan") {
      vi <- viClo(beta, Y, delta, X, id, weights, B, stratify = stratify)$vi * n
  }
  if (omega != TRUE && vClose == TRUE && rankweights == "logrank") {
      vi <- viClo(beta, Y, delta, X, id, weights, B, "logrank", stratify = stratify)$vi * n
  }
  if (omega != TRUE && vClose != TRUE) {
      vi <- viEmp(beta, Y, delta, X, id, weights, B, mb = TRUE, zbeta = FALSE, rankweights = rankweights)$vi
  }
  if (rankweights == "gehan") {
      An <- abargehanfun(beta, Y, X, delta, clsize, sigma, weights)
  }
  if (rankweights == "logrank") {
      An <- abarlogfun(beta, Y, X, delta, clsize, sigma, weights)
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

zlFun <- function(beta, Y, delta, X, id, weights = rep(1, nrow(X)), B = 500, vClose = FALSE, rankweights = "gehan", stratify = TRUE, sigma = diag(ncol(X))) {
  p <- ncol(X)
  clsize <- unlist(lapply(split(id, id), length))
  n <- sum(unlist(lapply(split(weights, id), unique)))
  # n <- length(clsize)
  # N <- sum(clsize)
  UnMat <- zmat <- ahat <- unTime <- NULL
  An.inv <- 1
  UnV <- viEmp(beta, Y, delta, X, id, weights, B, mb = FALSE, zbeta = TRUE, smooth = FALSE, rankweights = rankweights, sigma = sigma)
  zmat <- UnV$zmat
  UnV <- UnV$UnV
  if (vClose == TRUE) {
    vi <- viClo(beta, Y, delta, X, id, weights, B, rankweights, stratify = stratify)$vi * n
  }
  if (vClose != TRUE) {
    vi <- viEmp(beta, Y, delta, X, id, weights, B, mb = TRUE, zbeta = FALSE, smooth = FALSE, rankweights = rankweights)$vi
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
  covmat <- covmat  / n
  covmat <- matrix(as.numeric(covmat), p)
  list(covmat = covmat, vi = vi, An.msg = An.msg, An.inv = An.inv)
}

huangFun <- function(beta, Y, delta, X, id, weights = rep(1, nrow(X)), B = 500, vClose = TRUE, rankweights = "gehan", sigma = diag(ncol(X)), stratify = TRUE) {
  p <- ncol(X)
  clsize <- unlist(lapply(split(id, id), length))
  n <- sum(unlist(lapply(split(weights, id), unique)))
  N <- sum(weights)
  ## n <- length(clsize)
  ## N <- sum(clsize)
  betaTemp <- NULL
  UnMatV <- NULL
  if (vClose == TRUE) {
      vi <- viClo(beta, Y, delta, X, id, weights, B, rankweights, stratify = stratify)$vi * n
  }
  if (vClose != TRUE) {
      vi <- viEmp(beta, Y, delta, X, id, weights, B, mb = TRUE, zbeta = FALSE, smooth = FALSE, rankweights = rankweights)$vi
  }
  qq <- chol(vi)
  qq <- t(qq)
  newBeta <- NULL
  Z <- rep(1, sum(clsize)) #
  newBeta <- matrix(0, ncol = p, nrow = p)
  for ( i in 1:p) {
      if (rankweights == "gehan") {
          bb <- BBsolve(beta, uiFun, Y = Y, X = X, delta = delta, clsize = clsize, sigma = sigma, n = length(clsize), Z = Z, weights = weights, smooth = TRUE, constant = qq[,i] * n ^ -0.5, quiet = TRUE)
          newBeta[, i] <- bb$par
      }
      if (rankweights == "logrank") {
          bb <- BBsolve(beta, uilogFun, Y = Y, X = X, delta = delta, clsize = clsize, sigma = sigma, n = length(clsize), Z = Z, weights = weights, smooth = TRUE, constant = qq[,i] * n ^ -0.5, rankweights = rankweights, quiet = TRUE)
          newBeta[, i] <- bb$par
      }
      if (rankweights == "PW" | rankweights == "Prentice-Wilcoxon") {
          bb <- BBsolve(beta, uilogFun, Y = Y, X = X, delta = delta, clsize = clsize, sigma = sigma, n = length(clsize), Z = Z, weights = weights, smooth = TRUE, constant = qq[,i] * n ^ -0.5, rankweights = rankweights, quiet = TRUE)
          newBeta[, i] <- bb$par
      }
  }
  dd <- newBeta - beta
  covmat <- n * t(dd) %*% dd
##  covmat <- covmat / n
  list(covmat = covmat)
}

smoothrr <- function(formula, data, subset, contrasts = NULL, id,
                     weights = NULL, rankweights = "gehan",
                     binit = "lm", sigmainit = NULL,
                     variance = "ISMB",
                     B = 100,
                     strataid = NULL,
                     control = aftgee.control()) {
  if (sum(!(variance %in% c("MB", "ZLCF", "ZLMB", "sHCF", "sHMB", "ISCF", "ISMB", "js"))) == 8) stop("Invaid variance estimates.")
  if (!(rankweights %in% c("gehan", "logrank", "PW", "Prentice-Wilcoxon", "GP"))) stop("Invaid rankweights weights.")
  scall <- match.call()
  mnames <- c("", "formula", "data", "weights", "subset", "na.nation", "id", "strataid")
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
  weights <- model.extract(m, weights)
  strataid <- model.extract(m, strataid)
  if (is.null(weights)) weights <- rep(1, N)
  stratify <- TRUE
  if (is.null(strataid)) stratify <- FALSE
  xnames <- colnames(x) ## [-1]
  if(is.null(sigmainit)) {sigmainit = diag(ncol(x))} ## x[,-1]
  out <- smoothrr.fit(Y = log(y[,1]), delta = y[,2], X = as.matrix(x), id = id, weights = weights, binit = binit, sigmainit = sigmainit, variance = variance, B = B, rankweights = rankweights, stratify = stratify, control = control) ## x[,-1]
  out$call <- scall
  out$vari.name <- xnames
  class(out) <- "smoothrr"
  return(out)
}


smoothrr.fit <- function(Y, delta, X, binit = "lm",
                         sigmainit, id, weights = rep(1, nrow(X)),
                         variance = c("MB", "ZLCF", "ZLMB", "sHCF", "sHMB", "ISCF", "ISMB", "js"),
                         B = 100, rankweights = "gehan", stratify = TRUE,
                         control = aftgee.control()) {
  p <- ncol(X)
  go <- 1
  step <- 0
  if (is.numeric(binit)) {
      if (length(binit) != p) {
          stop("Binit value length does not match with numbers of covariates.")
      }
      beta <- binit
  }
  if (!(is.numeric(binit))) {
      if (binit == "lm") {
          beta <- lm(Y ~ X - 1)$coef
      }
  }
  beta.temp1 <- beta
  sigmainit.temp1 <- sigmainit
  beta.temp2 <- NULL
  sigmainit.temp2 <- NULL
  clsize <- unlist(lapply(split(id, id), length))
  order <- unlist(lapply(clsize, function(x) 1:x))
  N <- sum(clsize)
  n <- length(clsize)
  Z <- rep(1, N)
  var.MB <- var.ZLCF <- var.ZLMB <- var.sHCF <- var.sHMB <- var.js <- var.ISCF <- var.ISMB <- NaN
  for (i in 1:control$maxiter) {
      ## Point Estimation
      if (rankweights %in% c("gehan", "PW", "Prentice-Wilcoxon", "GP") & i == 1) {
          tbeta <- system.time(beta.temp2 <- nlm(Ln, p = beta.temp1, other = list(Y, X, delta, clsize, sigmainit.temp1, weights, Z), fscale = 0.01)$estimate)
          ## tbeta <- system.time(temp <- BBsolve(beta.temp1, uiFun, Y = Y, X = X, delta = delta, clsize = clsize, sigma = sigmainit.temp1, n = n, Z = Z, weights = weights, smooth = TRUE, constant = 0, quiet = TRUE))
          ## beta.temp2 <- temp$par
          ## beta.conv <- temp$convergence
          beta.conv <- 0
      }
      if (rankweights == "logrank") {
          tbeta <- system.time(temp <- BBsolve(beta.temp1, uilogFun, Y = Y, X = X, delta = delta, clsize = clsize, sigma = sigmainit.temp1, n = n, Z = Z, weights = weights, smooth = TRUE, constant = 0, rankweights = rankweights, quiet = TRUE))
          beta.temp2 <- temp$par
          beta.conv <- temp$convergence
      }
      if (rankweights == "PW" | rankweights == "Prentice-Wilcoxon") {
          beta.temp1 <- beta.temp2
          tbeta <- system.time(temp <- BBsolve(beta.temp1, uilogFun, Y = Y, X = X, delta = delta, clsize = clsize, sigma = sigmainit.temp1, n = n, Z = Z, weights = weights, smooth = TRUE, constant = 0, rankweights = rankweights, quiet = TRUE))
          beta.temp2 <- temp$par
          print(beta.temp1)

          ## if (abs(beta.temp1 - beta.temp2) < control$abstol) {break}
          if (max(abs(beta.temp1 - beta.temp2)) >= control$abstol) {
              beta.temp1 <- beta.temp2
              next
          }
          beta.conv <- ifelse(i == control$maxiter, 1, 0)
      }
      if (rankweights == "GP") {
          beta.temp1 <- beta.temp2
          tbeta <- system.time(temp <- BBsolve(beta.temp1, uilogFun, Y = Y, X = X, delta = delta, clsize = clsize, sigma = sigmainit.temp1, n = n, Z = Z, weights = weights, smooth = TRUE, constant = 0, rankweights = rankweights, quiet = TRUE))
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
      if (B == 0) {break}
      ## begin multipler resampling method
      if (sum(variance %in% "MB") > 0 && pass.MB == FALSE) {
          temp <- matrix(0, nrow = B, ncol = p)
          if(B == 0) {
              var.MB <- NULL
              pass.MB <- TRUE
              if (sum(c("ZLCF", "js", "ZLMB", "sHCF", "sHMB", "ISCF", "ISMB") %in% variance) == 0) {break}
          }
          for (i in 1:B) {
              Z <- rep(rexp(length(clsize)), clsize)

              if (rankweights == "gehan") {
                  temp[i,] <- nlm(Ln, p= beta.temp2, other = list(Y, X, delta, clsize, sigmainit.temp1, weights, Z), fscale = 0.01)$estimate
              }
              if (rankweights == "logrank") {
                  temp[i,] <- BBsolve(beta.temp2, uilogFun, Y = Y, X = X, delta = delta, clsize = clsize, sigma = sigmainit.temp1, n = n, Z = Z, weights = weights, smooth = TRUE, constant = 0, rankweights = rankweights, quiet = TRUE)$par
              }
              if (rankweights == "PW" | rankweights == "Prentice-Wilcoxon") {
                  temp[i,] <- BBsolve(beta.temp2, uilogFun, Y = Y, X = X, delta = delta, clsize = clsize, sigma = sigmainit.temp1, n = n, Z = Z, weights = weights, smooth = TRUE, constant = 0, rankweights = rankweights, quiet = TRUE)$par
              }
              if (rankweights == "GP") {
                  temp[i,] <- BBsolve(beta.temp2, uilogFun, Y = Y, X = X, delta = delta, clsize = clsize, sigma = sigmainit.temp1, n = n, Z = Z, weights = weights, smooth = TRUE, constant = 0, rankweights = rankweights, quiet = TRUE)$par
              }
          }
          var.MB <- var(temp)
          pass.MB <- TRUE
          if (sum(c("ZLCF", "js", "ZLMB", "sHCF", "sHMB", "ISCF", "ISMB") %in% variance) == 0) {break}
      } ## end MBipler resampling method

      ## begin Zeng and Lin's approach, close Si
      if (sum(variance %in% "ZLCF") > 0 && pass.ZLCF == FALSE) {
          var.ZLCF <- zlFun(beta.temp2, Y, delta, X, id, weights, B, vClose = TRUE, rankweights, stratify = stratify, sigma = sigmainit.temp1)
          ## if (var.ZLCF$An.inv == 0) {print (var.ZLCF$An.msg)}
          ZLCF.An.inv <- var.ZLCF$An.inv
          var.ZLCF <- var.ZLCF$covmat
          pass.ZLCF <- TRUE
          if (sum(c("js", "ZLMB", "sHCF", "sHMB", "ISCF", "ISMB") %in% variance) == 0) {break}
      } ## end Zeng and Lin's approach, close Si

      ## begin Zeng and Lin's approach, emprical Si
      if (sum(variance %in% "ZLMB") > 0 && pass.ZLMB == FALSE) {
          var.ZLMB <- zlFun(beta.temp2, Y, delta, X, id, weights, B, vClose = FALSE, rankweights, sigma = sigmainit.temp1)
          ## if(var.ZLMB$An.inv) {print(var.ZLMB$An.msg)}
          ZLMB.An.inv <- var.ZLMB$An.inv
          var.ZLMB <- var.ZLMB$covmat
          pass.ZLMB <- TRUE
          if (sum(c("js", "sHCF", "sHMB", "ISCF", "ISMB") %in% variance) == 0) {break}
      } ## end Zeng and Lin's approach, emprical Si

    ## begin Huang's approach, close Si
    if (sum(variance %in% "sHCF") > 0 && pass.sHCF == FALSE){
      var.sHCF <- huangFun(beta.temp2, Y, delta, X, id, weights, B, vClose = TRUE, rankweights, sigma = sigmainit.temp1, stratify = stratify)$covmat
      pass.sHCF <- TRUE
      if (sum(c("js", "sHMB", "ISCF", "ISMB") %in% variance) == 0) {break}
    } ## end Huang's approach, close Si

    ## begin Huang's approach, emprical Si
    if (sum(variance %in% "sHMB") > 0 && pass.sHMB == FALSE){
      var.sHMB <- huangFun(beta.temp2, Y, delta, X, id, weights, B, vClose = FALSE, rankweights, sigma = sigmainit.temp1)$covmat
      pass.sHMB <- TRUE
      if (sum(c("js", "ISCF", "ISMB") %in% variance) == 0) {break}
    } ## end Huang's approach, emprical Si

    ## begin IS's appraoch, close Si
    if (sum(variance %in% "ISCF") > 0 && pass.ISCF == FALSE) {
        var.ISCF <- isFun(beta.temp2, Y, delta, X, id, weights, sigmainit.temp1, B, vClose = TRUE, rankweights, stratify = stratify)
        ## if(var.ISCF$An.inv == 0){print(var.ISCF$An.msg)}
        ISCF.An.inv <- var.ISCF$An.inv
        var.ISCF <- var.ISCF$covmat
        if (sum(c("js", "ISMB") %in% variance) == 0) {break}
        pass.ISCF <- TRUE
    }
    ## end IS's approach, close Si

    ## begin IS's appraoch, emprical Si
    if (sum(variance %in% "ISMB") > 0 && pass.ISMB == FALSE) {
        var.ISMB <- isFun(beta.temp2, Y, delta, X, id, weights, sigmainit.temp1, B, vClose = FALSE, rankweights)
        ## if(var.ISMB$An.inv == 0) {print(var.ISMB$An.msg)}
        ISMB.An.inv <- var.ISMB$An.inv
        var.ISMB <- var.ISMB$covmat
        if (!("js" %in% variance)) {break}
        pass.ISMB <- TRUE
    }
    ## end IS's approach, emprical Si

    ## begin JS's iterative approach
    if (sum(variance %in% "js") > 0) {
      var.js <- isFun(beta.temp2, Y, delta, X, id, weights, sigmainit.temp1, omega = TRUE, rankweights = rankweights, stratify = stratify)
      ## if(var.js$An.inv == 0) {print(var.js$An.msg)}
      js.An.inv <- var.js$An.inv
      sigmainit.temp2 <- var.js$covmat
      e.sigma <- abs(sigmainit.temp1-sigmainit.temp2)
      e.beta <- abs(beta.temp1 - beta.temp2)
      e.rel <- max(max(e.sigma/abs(sigmainit.temp2)), max(e.beta/abs(beta.temp2)))
      e.abs <- max(e.sigma, e.beta)
      if (abs(e.rel) < control$reltol) break
      if (abs(e.abs) < control$abstol) break
      ## mat add an absolute tolerance check
      beta.temp1 <- beta.temp2
      sigmainit.temp1 <- sigmainit.temp2
      ## var.js <- sigmainit.temp2 / n
    } ## end js
  }
  covmat <- list(MB = var.MB, ZLCF = var.ZLCF, ZLMB = var.ZLMB, sHCF = var.sHCF, sHMB = var.sHMB, ISCF = var.ISCF, ISMB = var.ISMB, js = sigmainit.temp2 / n)
  convergence <- if ( i == control$maxiter) 1 else 0
  out <- list(beta = beta.temp2, covmat = covmat, convergence = convergence, tbeta = tbeta, beta.conv = beta.conv, beta.conv.at = i, var.meth = variance,
              ZLMB.An.inv = ZLMB.An.inv, ZLCF.An.inv = ZLCF.An.inv, ISMB.An.inv = ISMB.An.inv, ISCF.An.inv = ISCF.An.inv, js.An.inv = js.An.inv)
  ## if (ZLMB.An.inv == 0 | ZLCF.An.inv == 0 | ISMB.An.inv == 0 | ISCF.An.inv == 0 | js.An.inv == 0) {
  ##      print("An is singular")
  ## }
  class(out) <- "smoothrr.fit"
  return(out)
}



