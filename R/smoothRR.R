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

abarlogfun <- function(beta, Y, X, delta, clsize, sigma, weights, pw = rep(1, nrow(X))) {
    p <- ncol(X)
    N <- nrow(X)
    n <- length(clsize)
    a <- vector("double", p * p)
    matrix(.C("abarlogfun", as.double(beta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize),
              as.double(pw), as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(weights),
              out = as.double(a), PACKAGE = "aftgee")$out, nrow = p)
}

abarpwfun <- function(beta, Y, X, delta, clsize, sigma, weights, pw) {
    p <- ncol(X)
    N <- nrow(X)
    n <- length(clsize)
    a <- vector("double", p * p)
    pt1 <- matrix(.C("abarpwfun", as.double(beta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize), as.double(pw$fhat), as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(weights), out = as.double(a), PACKAGE = "aftgee")$out, p)
    ## pt1 <- uilogFun(beta, Y, X, delta, clsize, sigma, n, Z = rep(1, nrow(X)), weights, smooth = TRUE, constant = 0, s = 0, pw = pw$fhat)
    pt2 <- abarlogfun(beta, Y, X, delta, clsize, sigma, weights, pw$Shat)
    ## rep(1, p) %o% pt1 + pt2
    ## diag(pt1) + pt2
    pt1 + pt2
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
  ans <- .C("ufun", as.double(beta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize), as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(Z), as.double(weights), out = as.double(sn), PACKAGE = "aftgee")$out
  ans <- ans - constant
}

uilogFun <- function(beta, Y, X, delta, clsize, sigma, n, Z, weights, smooth = TRUE, constant = 0, gpweight, pw, rankweights) {
  N <- sum(clsize)
  p <- ncol(X)
  sn <- vector("double", p)
  ans <- numeric(p)
  if (gpweight > 0) {
    pw <- get_suv(Y, X, beta, N, delta, weights)$Shat ^ gpweight
  }
  if (gpweight > 0 & rankweights == "eGP") {
      pw <- get_suv(Y, X, beta, N, delta, weights)
      pw <- ((pw$Shat - pw$Shatlast) / (1 - pw$Shatlast))^ gpweight
  }
  if (smooth == TRUE) {
      ans <- .C("ulogfun", as.double(beta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize),
                as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(Z), as.double(weights),
                as.double(pw), out = as.double(sn), PACKAGE = "aftgee")$out
  }
  if (smooth != TRUE) {
      ans <- .C("ulognsfun", as.double(beta), as.double(Y), as.double(X), as.double(delta), as.integer(clsize),
                as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(Z), as.double(weights),
                as.double(pw), out = as.double(sn), PACKAGE = "aftgee")$out
  }
  ans - constant
}

viEmp <- function(beta, Y, delta, X, id, weights = rep(1, nrow(X)), B = 500, mb = TRUE, zbeta = FALSE, smooth = TRUE, rankweights = "gehan", sigma = diag(ncol(X)), gpweight = 1){
  p <- ncol(X)
  clsize <- unlist(lapply(split(id, id), length))
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
    ## if (rankweights %in% c("nsPW", "nsGP")) {
    ##     smooth <- FALSE
    ## }
    gpweight <- ifelse(rankweights %in% c("nsGP", "GP"), gpweight, 1)
    if (smooth == TRUE) {
        if (rankweights == "gehan") {
            UnV[,i] <- as.vector(.C("ufun", as.double(newbeta), as.double(Y), as.double(X), as.double(delta),
                                    as.integer(clsize), as.double(sigma), as.integer(length(clsize)), as.integer(p),
                                    as.integer(sum(clsize)), as.double(Z), as.double(weights), out = as.double(sn),
                                    PACKAGE = "aftgee")$out) # / n
        }
        if (rankweights %in% c("nslogrank", "logrank")) {
            UnV[,i] <- uilogFun(newbeta, Y, X, delta, clsize, sigma, length(clsize), Z, weights, smooth, constant = 0,
                                gpweight = 0, pw = rep(1, sum(clsize)), rankweights) # / n
        }
        if (rankweights %in% c("PW", "Prentice-Wilcoxon", "GP", "eGP")) {
            UnV[,i] <- uilogFun(beta = newbeta, Y = Y, X = X, delta = delta, clsize = clsize, sigma = sigma,
                                n = length(clsize), Z = Z, weights = weights, smooth = TRUE, constant = 0,
                                gpweight = gpweight, pw = rep(1, sum(clsize)), rankweights = rankweights) # / n
        }
    }
    if (smooth == FALSE) {
        if (rankweights == "gehan") {
            UnV[,i] <- as.vector(.C("unsfun", as.double(newbeta), as.double(Y), as.double(X), as.double(delta),
                                    as.integer(clsize), as.double(sigma), as.integer(length(clsize)), as.integer(p),
                                    as.integer(sum(clsize)), as.double(Z), as.double(weights), out = as.double(sn),
                                    PACKAGE = "aftgee")$out) # / n ## n and N?
        }
        if (rankweights %in% c("logrank", "nslogrank")) {
            UnV[,i] <- uilogFun(newbeta, Y, X, delta, clsize, sigma, length(clsize), Z, weights, smooth, constant = 0,
                                gpweight = 0, pw = rep(1, sum(clsize)), rankweights) # / n
        }
        if (rankweights %in% c("nsPW", "nsGP", "PW", "GP", "eGP")) {
            UnV[,i] <- uilogFun(beta = newbeta, Y = Y, X = X, delta = delta, clsize = clsize, sigma = sigma,
                                n = length(clsize), Z = Z, weights = weights, smooth = FALSE, constant = 0,
                                gpweight = gpweight, pw = rep(1, sum(clsize)), rankweights = rankweights) # / n
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
    xdif <- X[ik,] - X[jl,]
    edif <- en[ik] - en[jl]
    ind <- ifelse(edif <= 0, 1, 0)
    minEn <- ifelse(edif <= 0, ik, jl)
    si <- s <- NULL
    if (rankweights == "gehan") {
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
        s <- rowsum(s, ik)
        }
    if (rankweights == "logrank") {
        haz <- -1 * log(Shat)
        haz <- ifelse(haz == Inf, 0, haz)
        haz <- ifelse(haz == -Inf, 0, haz)
        haz <- ifelse(is.na(haz) == TRUE, 0, haz)
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
    s
}

get_vi <- function(s, id, delta, weights, n) {
    clweights <- as.numeric(unlist(lapply(split(weights, id), unique)))
    s1 <- rowsum(s, group = id)
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
    list(vi = vi, v1 = v1, v2 = v2)
}



viClo <- function(beta, Y, delta, X, id, weights = rep(1, nrow(X)), B = 500, rankweights = "gehan", stratify = TRUE) {
    s <- v1 <- vi <- v2i <- NULL
    n <- sum(unlist(lapply(split(weights, id), unique)))
    p <- ncol(X)
    clsize <- unlist(lapply(split(id, id), length))
    clid <- unlist(lapply(clsize, function(x) 1:x))
    stra <- match(weights, unique(weights))
    dim <- unique(clsize)
    s <- get_si(beta, Y, delta, X, id, weights, B, rankweights)
    if (stratify == TRUE) {
        v1 <- get_vi(s, id, delta, weights, n)$v1
        vi <- v1
        v2 <- matrix(0, ncol = p, nrow = p)
        if (length(unique(stra)) > 1) {
            for (i in 1:length(unique(stra))) {
                ns <- sum(unlist(lapply(split(weights[stra == i], id[stra == i]), unique)))
                ## weights at cluster level
                v2i <- get_vi(s = s[stra == i, ], id = id[stra == i], delta = delta[stra == i],
                              weights = weights[stra == i],
                              n = ns
                              )$vi
                strPr <- ns / n
                v2 <- v2 + v2i * strPr
            }
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
  if (rankweights %in% c("PW", "GP", "Prentice-Wilcoxon")) {
      pw <- get_smooth_suv(Y, X, beta, n, delta, weights)
      An <- abarpwfun(beta, Y, X, delta, clsize, sigma, weights, pw)
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

zlFun <- function(beta, Y, delta, X, id, weights = rep(1, nrow(X)), B = 500, vClose = FALSE, rankweights = "gehan", stratify = TRUE, sigma = diag(ncol(X)), gpweight = 1, variance) {
  p <- ncol(X)
  clsize <- unlist(lapply(split(id, id), length))
  n <- sum(unlist(lapply(split(weights, id), unique)))
  # n <- length(clsize)
  # N <- sum(clsize)
  smooth <- ifelse(variance %in% c("sZLMB"), TRUE, FALSE)
  UnMat <- zmat <- ahat <- unTime <- NULL
  An.inv <- 1
  UnV <- viEmp(beta, Y, delta, X, id, weights, B, mb = FALSE, zbeta = TRUE, smooth = smooth, rankweights = rankweights, sigma = sigma, gpweight = gpweight)
  zmat <- UnV$zmat
  UnV <- UnV$UnV
  if (vClose == TRUE) {
    vi <- viClo(beta, Y, delta, X, id, weights, B, rankweights, stratify = stratify)$vi * n
  }
  if (vClose != TRUE) {
    vi <- viEmp(beta, Y, delta, X, id, weights, B, mb = TRUE, zbeta = FALSE, smooth = smooth, rankweights = rankweights, gpweight = gpweight)$vi ## smooth
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
          bb <- BBsolve(beta, uiFun, Y = Y, X = X, delta = delta, clsize = clsize, sigma = sigma, n = length(clsize), Z = Z,
                        weights = weights, smooth = TRUE, constant = qq[,i] * n ^ -0.5, quiet = TRUE)
          newBeta[, i] <- bb$par
      }
      if (rankweights == "logrank") {
          bb <- BBsolve(beta, uilogFun, Y = Y, X = X, delta = delta, clsize = clsize, sigma = sigma, n = length(clsize), Z = Z,
                        weights = weights, smooth = TRUE, constant = qq[,i] * n ^ -0.5, gpweight = 0, pw = rep(1, sum(clsize)), rankweights = rankweights, quiet = TRUE)
          newBeta[, i] <- bb$par
      }
      if (rankweights %in% c("PW", "Prentice-Wilcoxon", "GP")) {
          gpweight <- ifelse(rankweights == "GP", ncol(X), 1)
          bb <- BBsolve(beta, uilogFun, Y = Y, X = X, delta = delta, clsize = clsize, sigma = sigma, n = length(clsize), Z = Z,
                        weights = weights, smooth = TRUE, constant = qq[,i] * n ^ -0.5, gpweight, pw = rep(1, sum(clsize)), rankweights = rankweights, quiet = TRUE)
          newBeta[, i] <- bb$par
      }
  }
  dd <- newBeta - beta
  covmat <- n * t(dd) %*% dd
##  covmat <- covmat / n
  list(covmat = covmat)
}

solveGehan <- function(beta, Y, X, delta, clsize, sigma, weights, Z) {
    tbeta <- system.time(temp <- nlm(Ln, p = beta, other = list(Y, X, delta, clsize, sigma, weights, Z), fscale = 0.01))
    conv <- ifelse(temp$code == 1, 0, 1)
    out <- list(beta2 = temp$estimate, tbeta = tbeta, conv = conv)
}

solveLogRank <- function(beta, Y, X, delta, clsize, sigma, weights, Z, n, smooth, constant, rankweights){
    tbeta <- system.time(temp <- BBsolve(beta, uilogFun, Y = Y, X = X, delta = delta, clsize = clsize,
                                         sigma = sigma, n = n, Z = Z, weights = weights, smooth = smooth,
                                         constant = 0, gpweight = 0, pw = rep(1, sum(clsize)), rankweights = rankweights, quiet = TRUE))
    out <- list(beta2 = temp$par, tbeta = tbeta, conv = temp$convergence)
}

solvePW <- function(beta, Y, X, delta, clsize, sigma, weights, Z, n, smooth, rankweights, constant, gpweight, control = aftgee.control()){
    gpweight <- ifelse(rankweights %in% c("GP", "nsGP"), gpweight, 1)
    temp <- BBsolve(beta, uilogFun, Y = Y, X = X, delta = delta, clsize = clsize, sigma = sigma, n = n, Z = Z,
                    weights = weights, smooth = smooth, constant = 0, gpweight = gpweight, pw = rep(1, sum(clsize)), rankweights = rankweights, quiet = TRUE)
    out <- list(beta2 = temp$par, conv = temp$convergence)
}


solvePWiter <- function(beta, Y, X, delta, clsize, sigma, weights, Z, n, smooth, rankweights, constant, gpweight, control = aftgee.control()){
    gpweight <- ifelse(rankweights %in% c("GP", "nsGP"), gpweight, 1)
    convStep <- 0
    N <- sum(clsize)
    betahist <- NULL
    bbhist <- NULL

    for(i in 1:control$maxiter) {
        pw <- rep(1, N)
        ## if (phi == "gehan") {
        ##     pw <- get_gehan(Y, X, beta, N, delta, clsize, sigma, weights)
        ## }
        pw <- get_suv(Y, X, beta, N, delta, weights)
        if (rankweights %in% c("nsPW", "PW")) {
            pw <- pw$Shat
        }
        if (rankweights %in% c("nsGP", "GP")) {
            pw <- pw$Shat ^ gpweight
        }
        if (rankweights == "eGP") {
            pw <- ((pw$Shat - pw$Shatlast) / (1 - pw$Shatlast)) ^ gpweight
        }
        temp <- BBsolve(beta, uilogFun, Y = Y, X = X, delta = delta, clsize = clsize,
                        sigma = sigma, n = n, Z = Z, weights = weights, smooth = smooth,
                        constant = 0, gpweight = 0, pw = pw, rankweights = rankweights, quiet = TRUE)
        beta2 <- temp$par
        betahist <- rbind(betahist, beta2)
        bbhist <- c(bbhist, temp$convergence)
        ## print(beta2)
        ## print(temp$convergence)
        if (max(abs(beta2 - beta)) <= control$abstol | max(abs(beta2 - beta) / abs(beta2)) <= control$reltol) {
            beta <- beta2
            convStep <- convStep + 1
            conv <- 0
            break
        }
            beta <- beta2
            convStep <- convStep + 1
            conv <- 1
    }
    out <- list(beta2 = beta2, conv = conv, step = convStep, betahist = betahist, bbhist = bbhist)
}


smoothrr <- function(formula, data, subset, contrasts = NULL, id,
                     weights = NULL, rankweights = "gehan",
                     binit = "lm", sigmainit = NULL,
                     variance = "ISMB",
                     B = 100, gpweight = 1,
                     strataid = NULL, iter = TRUE,
                     control = aftgee.control()) {
  if (sum(!(variance %in% c("MB", "ZLCF", "sZLMB", "ZLMB", "sHCF", "sHMB", "ISCF", "ISMB", "js"))) == 9) stop("Invaid variance estimates.")
  if (!(rankweights %in% c("gehan", "nslogrank", "logrank", "PW", "Prentice-Wilcoxon", "GP", "nsGP", "nsPW", "eGP"))) stop("Invaid rankweights weights.")
  scall <- match.call()
  mnames <- c("", "formula", "data", "weights", "subset", "na.nation", "id", "strataid")
  cnames <- names(scall)
  cnames <- cnames[match(mnames, cnames, 0)]
  mcall <- scall[cnames]
  mcall[[1]] <- as.name("model.frame")
  m <- eval(mcall, parent.frame())
  id <- model.extract(m, id)
  if (is.null(id)) {
      id <- 1:nrow(y)
  }
  y <- model.extract(m, "response")
  N <- NROW(y)
  mterms <- attr(m, "terms")
  x <- model.matrix(mterms, m, contrasts) ## this x has interception
  weights <- model.extract(m, weights)
  strataid <- model.extract(m, strataid)
  if (is.null(weights)) weights <- rep(1, N)
  stratify <- TRUE
  if (is.null(strataid)) stratify <- FALSE
  xnames <- colnames(x)
  if(is.null(sigmainit)) {sigmainit = diag(ncol(x))}
  out <- smoothrr.fit(Y = log(y[,1]), delta = y[,2], X = as.matrix(x), id = id,
                      weights = weights, binit = binit, sigmainit = sigmainit,
                      variance = variance, B = B, rankweights = rankweights,
                      stratify = stratify, control = control, iter = iter, gpweight = gpweight)
  out$call <- scall
  out$vari.name <- xnames
  class(out) <- "smoothrr"
  return(out)
}


smoothrr.fit <- function(Y, delta, X, binit = "lm",
                         sigmainit, id, weights = rep(1, nrow(X)),
                         variance = c("MB", "ZLCF", "ZLMB", "sZLMB", "sHCF", "sHMB", "ISCF", "ISMB", "js"),
                         B = 100, rankweights = "gehan", stratify = TRUE, iter, gpweight = gpweight,
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
    betahist <- bbhist <- NULL
    beta.temp1 <- beta
    sigmainit.temp1 <- sigmainit
    beta.temp2 <- NULL
    sigmainit.temp2 <- NULL
    clsize <- unlist(lapply(split(id, id), length))
    order <- unlist(lapply(clsize, function(x) 1:x))
    N <- sum(clsize)
    n <- length(clsize)
    Z <- rep(1, N)
    var.MB <- var.ZLCF <- var.ZLMB <- var.sHCF <- var.sHMB <- var.js <- var.ISCF <- var.ISMB <- beta.step <- NaN
    for (i in 1:control$maxiter) {
        ## Point Estimation
        if (rankweights %in% c("gehan", "PW", "Prentice-Wilcoxon", "GP", "nsGP", "nsPW", "eGP") & i == 1) {
            ans <- solveGehan(beta.temp1, Y, X, delta, clsize, sigmainit.temp1, weights, Z)
            tbeta <- ans$tbeta
            beta.temp2 <- ans$beta2
            beta.conv <- ans$conv
        }
        if (rankweights == "nslogrank") {
            ans <- solveLogRank(beta.temp1, Y, X, delta, clsize, sigmainit.temp1, weights, Z, n, FALSE, 0, rankweights)
            tbeta <- ans$tbeta
           beta.temp2 <- ans$beta2
            beta.conv <- ans$conv
        }
        if (rankweights == "logrank") {
            ans <- solveLogRank(beta.temp1, Y, X, delta, clsize, sigmainit.temp1, weights, Z, n, TRUE, 0, rankweights)
            tbeta <- ans$tbeta
           beta.temp2 <- ans$beta2
            beta.conv <- ans$conv
        }
        if (rankweights %in% c("PW", "Prentice-Wilcoxon", "GP", "nsPW", "nsGP", "eGP")) {
            beta.temp1 <- beta.temp2
            if (iter == FALSE) {
                tbeta <- system.time(ans <- solvePW(beta = beta.temp1, Y = Y, X = X, delta = delta, clsize = clsize,
                                                    sigma = sigmainit.temp1, weights = weights, Z = Z, n = n,
                                                    smooth = TRUE, rankweights = rankweights, constant = 0, gpweight = gpweight))
                beta.temp2 <- ans$beta2
                beta.conv <- ans$conv
                beta.step <- 0
            }
            if (iter == TRUE) {
                tbeta <- system.time(ans <- solvePWiter(beta = beta.temp1, Y = Y, X = X, delta = delta, clsize = clsize,
                                                        sigma = sigmainit.temp1, weights = weights, Z = Z, n = n,
                                                        smooth = TRUE, rankweights = rankweights, constant = 0, gpweight = gpweight))
                betahist <- ans$betahist
                bbhist <- ans$bbhist
                beta.temp2 <- ans$beta2
                beta.conv <- ans$conv
                beta.step <- ans$step
            }
        }
        if (rankweights %in% c("nsPW", "nsGP")) {
            beta.temp1 <- beta.temp2
            if (iter == FALSE) {
                tbeta <- system.time(ans <- solvePW(beta = beta.temp1, Y = Y, X = X, delta = delta, clsize = clsize,
                                                    sigma = sigmainit.temp1, weights = weights, Z = Z, n = n,
                                                    smooth = FALSE, rankweights = rankweights, constant = 0, gpweight = gpweight))
                beta.temp2 <- ans$beta2
                beta.conv <- ans$conv
                beta.step <- 0
            }
            if (iter == TRUE) {
                tbeta <- system.time(ans <- solvePWiter(beta = beta.temp1, Y = Y, X = X, delta = delta, clsize = clsize,
                                                        sigma = sigmainit.temp1, weights = weights, Z = Z, n = n,
                                                        smooth = FALSE, rankweights = rankweights, constant = 0, gpweight = gpweight))
                betahist <- ans$betahist
                bbhist <- ans$bbhist
                beta.temp2 <- ans$beta2
                beta.conv <- ans$conv
                beta.step <- ans$step
            }
        }
        beta.MB <- 0
        pass.MB <- pass.ZLCF <- pass.ZLMB <- pass.sHCF <- pass.sHMB <- pass.ISCF <- pass.ISMB <- FALSE
        ZLMB.An.inv <- ZLCF.An.inv <- ISMB.An.inv <- ISCF.An.inv <- js.An.inv <- 1
        if (B == 0) {break}
        ## begin MB resampling method
        if (sum(variance %in% "MB") > 0 && pass.MB == FALSE) {
            temp <- matrix(0, nrow = B, ncol = p)
            if(B == 0) {
                var.MB <- NULL
                pass.MB <- TRUE
                if (sum(c("ZLCF", "js", "ZLMB", "sZLMB", "sHCF", "sHMB", "ISCF", "ISMB") %in% variance) == 0) {break}
            }
            for (i in 1:B) {
                Z <- rep(rexp(length(clsize)), clsize)
                if (rankweights == "gehan") {
                    temp[i,] <- solveGehan(beta.temp2, Y, X, delta, clsize, sigmainit.temp1, weights, Z)$beta2
                }
                if (rankweights == "logrank") {
                    temp[i,] <- solveLogRank(beta.temp2, Y, X, delta, clsize, sigmainit.temp1, weights, Z, n, TRUE, 0, rankweights)$beta2
                }
                if (rankweights %in% c("PW", "Prentice-Wilcoxon", "GP")) {
                    if (iter == FALSE) {
                        temp[i,] <-solvePW(beta.temp1, Y, X, delta, clsize, sigmainit.temp1, weights, Z, n, TRUE, rankweights, 0, gpweight)$beta2
                    }
                    if (iter == TRUE) {
                        temp[i,] <-solvePWiter(beta.temp1, Y, X, delta, clsize, sigmainit.temp1, weights, Z, n, TRUE, rankweights, 0, gpweight)$beta2
                    }
                }
            }
            beta.MB <- temp
            var.MB <- var(temp)
            pass.MB <- TRUE
            if (sum(c("ZLCF", "js", "ZLMB", "sZLMB", "sHCF", "sHMB", "ISCF", "ISMB") %in% variance) == 0) {break}
        } ## end MB resampling method

        ## begin Zeng and Lin's approach, close Si
        if (sum(variance %in% "ZLCF") > 0 && pass.ZLCF == FALSE) {
            var.ZLCF <- zlFun(beta.temp2, Y, delta, X, id, weights, B, vClose = TRUE, rankweights, stratify = stratify, sigma = sigmainit.temp1, gpweight = gpweight)
            ZLCF.An.inv <- var.ZLCF$An.inv
            var.ZLCF <- var.ZLCF$covmat
            pass.ZLCF <- TRUE
            if (sum(c("js", "ZLMB", "sHCF", "sHMB", "ISCF", "ISMB") %in% variance) == 0) {break}
        } ## end Zeng and Lin's approach, close Si

        ## begin Zeng and Lin's approach, emprical Si
        if (sum(variance %in% c("sZLMB", "ZLMB")) > 0 && pass.ZLMB == FALSE) {
            var.ZLMB <- zlFun(beta.temp2, Y, delta, X, id, weights, B, vClose = FALSE, rankweights, sigma = sigmainit.temp1, gpweight = gpweight, variance = variance)
            ZLMB.An.inv <- var.ZLMB$An.inv
            var.ZLMB <- var.ZLMB$covmat
            pass.ZLMB <- TRUE
            ## make sZLMB into ZLMB... for now
            variance[which(variance == "sZLMB")] <- "ZLMB"
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
            ISCF.An.inv <- var.ISCF$An.inv
            var.ISCF <- var.ISCF$covmat
            if (sum(c("js", "ISMB") %in% variance) == 0) {break}
            pass.ISCF <- TRUE
        }
        ## end IS's approach, close Si

        ## begin IS's appraoch, emprical Si
        if (sum(variance %in% "ISMB") > 0 && pass.ISMB == FALSE) {
            var.ISMB <- isFun(beta.temp2, Y, delta, X, id, weights, sigmainit.temp1, B, vClose = FALSE, rankweights)
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
            beta.temp1 <- beta.temp2
            sigmainit.temp1 <- sigmainit.temp2
        }
    }
    covmat <- list(MB = var.MB, ZLCF = var.ZLCF, ZLMB = var.ZLMB, sHCF = var.sHCF, sHMB = var.sHMB, ISCF = var.ISCF, ISMB = var.ISMB, js = sigmainit.temp2 / n)
    convergence <- if ( i == control$maxiter) 1 else 0
    out <- list(beta = beta.temp2, covmat = covmat, convergence = convergence, tbeta = tbeta, beta.conv = beta.conv, beta.step = beta.step,
                beta.conv.at = i, betaMB = beta.MB, var.meth = variance,
                ZLMB.An.inv = ZLMB.An.inv, ZLCF.An.inv = ZLCF.An.inv, ISMB.An.inv = ISMB.An.inv, ISCF.An.inv = ISCF.An.inv, js.An.inv = js.An.inv, betahist = betahist, bbhist = bbhist)
    class(out) <- "smoothrr.fit"
    return(out)
}

get_suv <- function(Y, X, beta, N, delta, weights) {
    en <- Y - X %*% beta
    ik <- rep(1:N, each=N)
    jl <- rep(1:N, N)
    Shat <- fhat <- NULL
    dummy <- 1:N
    ord <- order(en)
    ei <- en[ord]
    weightsi <- weights[ord]
    deltai <- delta[ord]
    repeats <- table(ei)
    ## Shati <- survfit(Surv(ei, deltai) ~ 1, weights = weightsi)$surv
    Shati <- exp(-1 * basehaz(coxph(Surv(ei, deltai)~1, weights = weightsi))$hazard)
    ## fhati <- diff(c(Shati[1], Shati, Shati[length(Shati)]), 2) / diff(c(ei[1], ei, ei[length(ei)]), 2)
    Shati <- rep(Shati, repeats)
    Shatlast <- rev(Shati)[1]
    ## Shatlast <- Shati[1]
    ## fhati <- rep(fhati, repeats)
    Shat[dummy[ord]] <- Shati
    ## fhat[dummy[ord]] <- fhati
    ## fhat <- ifelse(fhat < -1, -1, fhat)
    ## mu <- mean(subset(fhat, fhat > -1))
    ## fhat <- ifelse(fhat < -1, mu, fhat)
    ## list(Shat = Shat, fhat = fhat, Shatlast = Shatlast)
    list(Shat = Shat, Shatlast = Shatlast)
}

get_gehan <- function(Y, X, beta, N, delta, clsize, sigma, weights) {
    p <- ncol(X)
    N <- nrow(X)
    n <- length(clsize)
    a <- vector("double", N)
    out <- matrix(.C("getgehan", as.double(beta), as.double(Y), as.double(X), as.integer(clsize),
                     as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(weights),
                     out = as.double(a), PACKAGE = "aftgee")$out, ncol = 1)
    out
}


get_smooth_suv <- function(Y, X, beta, N, delta, weights) {
    en <- Y - X %*% beta
    ik <- rep(1:N, each=N)
    jl <- rep(1:N, N)
    Shat <- fhat <- NULL
    dummy <- 1:N
    ord <- order(en)
    rij <- sqrt(diag((X[ord,] - X[dummy,]) %*% t((X[ord,] - X[dummy,]))))
    ei <- en[ord]
    weightsi <- weights[ord]
    deltai <- delta[ord]
    repeats <- table(ei)
    ## Shati <- survfit(Surv(ei, deltai) ~ 1, weights = weightsi)$surv
    ## assume no ties
    di <- rev(cumsum(rev(rep(1, N))))
    Z <- pnorm((en - ei) / rij)
    Z <- ifelse(is.na(Z) == T, 0, Z)
    hazi <- cumsum(deltai * Z / di)
    Shati <- exp(-1 * hazi)
    fhati <- diff(c(Shati[1], Shati, Shati[length(Shati)]), 2) / diff(c(ei[1], ei, ei[length(ei)]), 2)
    Shati <- rep(Shati, repeats)
    fhati <- rep(fhati, repeats)
    Shat[dummy[ord]] <- Shati
    fhat[dummy[ord]] <- fhati
    ## fhat <- ifelse(fhat < -1, -1, fhat)
    mu <- mean(subset(fhat, fhat > -1))
    fhat <- ifelse(fhat < -1, mu, fhat)
    list(Shat = Shat, fhat = fhat)
}
