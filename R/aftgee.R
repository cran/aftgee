aftgee <- function(formula, data, subset, id, contrasts = NULL,
                   weights = NULL,
                   margin = NULL,
                   corstr="independence",
                   binit = "lm", B = 100,
                   control = aftgee.control()
                   ) {
  scall <- match.call()
  mnames <- c("", "formula", "data", "weights", "margin", "subset", "na.action", "id")
  cnames <- names(scall)
  cnames <- cnames[match(mnames, cnames, 0)]
  mcall <- scall[cnames]
  ##  if (is.null(mcall$id)) mcall$id <- as.name("id")
  mcall[[1]] <- as.name("model.frame")
  m <- eval(mcall, parent.frame())
  y <- model.extract(m, "response")
  N <- NROW(y)
  mterms <- attr(m, "terms")
  x <- model.matrix(mterms, m, contrasts) ## ~x -1
  if (qr(x)$rank != ncol(x)) {
      stop("Design matirx does not have full rank; suggest removing intercept.")
  }
  weights <- model.extract(m, weights)
  if (is.null(weights)) weights <- rep(1, N)
  id <- model.extract(m, id)
  if (is.null(id)) stop("id variable not found.")
  margin <- model.extract(m, margin)
  if (is.null(margin)) margin <- rep(1, N)
  xnames <- colnames(x)[-1]
  out <- aftgee.fit(y = y, x = x, id = id, corstr = corstr, B = B, binit = binit,
                    weights = weights, margin = margin, control = control)
  out$call <- scall
  out$var.name <- xnames
  ## rownames(out$coefficients) <- rep(c("Intercept", xnames), length(unique(margin)))
  colnames(out$coefficients) <- c("binit", "AFTGEE")
  class(out) <- "aftgee"
  out
}

aftgee.fit <- function(y, x, id, corstr="independence",
                       weights = rep(1, nrow(x)),
                       margin = rep(1, nrow(x)),
                       B = 100, binit = "lm",
                       control = aftgee.control()) {
  x <- as.matrix(x)
  n <- length(unique(id))
  rm <- NULL
  rmName <- NULL
  p <- ncol(x)
  X <- x
  include.intercept <- ifelse(sum(x[,1]) == nrow(x), 1 , 0)
  if (include.intercept == 1) {
      p <- ncol(x) - 1
      X <- x[,-1]
  }
  firstBeta <- NULL
  firstSd <- NULL
  firstSdMat <- NULL
  clsize <- unlist(lapply(split(id, id), length))
  N <- sum(clsize)
  if (is.numeric(binit)) {
      if (length(binit) != p) {
          stop("Binit value length does not match with numbers of covariates.")
      }
  }
  if (!(is.numeric(binit))) {
      if (!(binit %in% c("lm", "srrgehan"))) {
          stop("Invaid binit value method.")
      }
  }
  Y <- log(y[,1])
  Y <- ifelse(Y == -Inf, 0, Y)
  delta <- y[,2]
  if (!(is.numeric(binit))) {
      if (binit == "lm") {
          linfit <- summary(lm(Y ~ x))
          first <- list(beta = linfit$coef[,1], sd = linfit$coef[,2])
          firstBeta <- first$beta
          firstSd <- first$sd
          firstconvergence <- first$convergence
      }

      if (binit == "srrgehan") {
          first <- smoothrr(Surv(exp(Y), delta) ~ X - 1, B = 0, variance = "ISMB", weights = weights, control = control, id = id)
          firstBeta <- first$beta
          firstSdMat <- first$covmat$ISMB
          if (is.na(firstSdMat)[1] == FALSE) {
              firstSd <- NaN
          }
          if (is.na(firstSdMat)[1] == FALSE) {
              firstSd <- sqrt(diag(firstSdMat))
          }
          firstconvergence <- first$convergence
          if (include.intercept == 1){
              firstBeta <- c(mean(Y - X %*% first$beta), firstBeta)
              firstSd <- c(NaN, firstSd)
              firstSdMat <- c(NaN, firstSdMat)
          }
      }
  }
  if (is.numeric(binit) & length(binit) == p) {
      firstBeta <- c(mean(Y - X %*% binit), binit)
      firstSd <- rep(NaN, p)
      firstconvergence <- 0
  }
  binitValue <- list(beta = firstBeta, sd = firstSd, sdMat = firstSdMat)
  result <- aftgee.est(Y, x, delta, binitValue$beta, id, corstr, rep(1, length(Y)), margin, weights, control)
  if (B > 0) {
    sample <- matrix(0, nrow = B, ncol = length(result$beta))
    for (i in 1:B){
      Z <- as.vector(rep(rexp(n,1), time = clsize))
      sample[i,] <- aftgee.est(Y, x, delta, result$beta, id, corstr, Z, margin, weights, control)$beta
    }
    vhat <- var(sample)
  }
  if (B == 0) {
    vhat <- NULL
  }
  ini.beta <- c(binitValue$beta)
  ini.sd <- c(binitValue$sd)
  ini.sdMat <- c(binitValue$sdMat)
  fit <- list(coefficients = cbind(ini.beta, result$beta),
              coef.res = result$beta,
              var.res = vhat,
              varMargin = result$gamma,
              alpha = result$alpha,
              coef.init = ini.beta,
              sd.init = ini.sd,
              var.init.mat = ini.sdMat,
              binit = binit,
              conv = result$convergence,
              ini.conv = firstconvergence,
              conv.step = result$convStep)
  class(fit) <- "aftgee.fit"
  fit
}

aftgee.control <- function(maxiter = 30,
                           reltol = 0.0001,
                           abstol = 0.0001,
                           trace = FALSE) {
  list(maxiter = maxiter,
       reltol = reltol,
       abstol = abstol,
       trace = trace)
}

aftgee.prepX <- function(xmat, margin, name = NULL) {
  newXmat <- NULL
  xmat <- as.matrix(xmat)
  p <- ncol(xmat)
  if (length(unique(margin)) == 1) {
      colnames(xmat) <- name
      xmat
  }
  else {
      dummy <- model.matrix(~factor(margin) -1 )
      xName <- NULL
      bName <- NULL
      for ( i in 1:nrow(xmat)) {
          newXmat <- rbind(newXmat, dummy[i,] %x% xmat[i,])
      }
      if (is.null(name) == TRUE & is.null(colnames(xmat)) == TRUE) {
          xmat <- newXmat
      }
      if (is.null(name) != TRUE | is.null(colnames(xmat)) != TRUE) {
          if (is.null(colnames(xmat)) != TRUE) {
              name <- colnames(xmat)
          }
          n.marg <- rep(unique(margin), each = p)
          name.marg <-rep(name, length(unique(margin)))
          new.name <- NULL
          for(i in 1:length(n.marg)){
              new.name <- c(new.name, paste(name.marg[i], n.marg[i], sep = ""))
          }
          xmat <- newXmat
          colnames(xmat) <- new.name
      }
      for( i in 1:length(names(xmat))) {
          colnames(xmat) <- gsub("(", "", colnames(xmat), fixed = TRUE)
          colnames(xmat) <- gsub(")", "", colnames(xmat), fixed = TRUE)
          colnames(xmat) <- gsub(":", "", colnames(xmat), fixed = TRUE)

      }
      xmat
  }
}

aftgee.est <- function(y, x, delta, beta, id, corstr="independence", Z = rep(1, length(y)),
                       margin = rep(1, length(id)), weights = rep(1, length(y)),
                       control = aftgee.control()) {
    xmat <- as.matrix(x) ## intercept included, pre-prepared.
    if (qr(x)$rank < ncol(x) & sum(x[,1]) == nrow(x)) {
        x <- x[,-1]
    }
    nobs <- length(y)
    for (i in 1:control$maxiter) {
        betaprev <- beta
        eres <- NULL
        eres2 <- NULL
        if (sum(margin == 1) == nobs) {
            e <- y - xmat %*% beta
            eres <- lss.eres(e, delta, Z * weights)
            yhat <- delta * y + (1 - delta) * (eres[[1]] + xmat %*% beta)
#            yhatZ <- sqrt(Z * weights) * yhat
#            xmatZ <- sqrt(Z * weights) * xmat
            geefit <- geese.fit(xmat, yhat, id, corstr = corstr)
        }
        if (sum(margin == 1) != nobs) {
            e <- y - xmat %*% beta
            ## if (res == FALSE) {
            ##     eres <- lss.eres(e, delta, Z * weights)
            ##     yhat <- delta * y + (1 - delta) * (eres[[1]] + xmat %*% beta)
            ##     yhatZ <- sqrt(Z * weights) * yhat
            ##     xmatZ <- sqrt(Z * weights) * xmat
            ##     geefit <- geese.fit(xmatZ, yhatZ, id, zsca = model.matrix(~factor(margin) -1 ), corstr = corstr)
            ## }
            ## else { # (res == TRUE) {
            er1 <- NULL
            er2 <- NULL
            for (m in unique(margin)) {
                temp <- lss.eres(e[margin == m], delta[margin == m], Z[margin == m])
                temp[[2]] <- ifelse(delta[margin == m] == 1, e[margin == m]^2, temp[[2]])
                eres2[m] <- mean(temp[[2]], na.rm = TRUE)
                dum <- cumsum(ifelse(margin == m, 1, 0))
                er1temp <- temp[[1]][ifelse(margin == m, dum, NA)]
                er1 <- rbind(er1, er1temp)
            }
            er1 <- as.vector(er1)
            er1 <- er1[!is.na(er1)]
            yhat <- delta * y + (1 - delta) * (er1 + xmat %*% beta)
#            yhatZ <- sqrt(Z * weights) * yhat
#            xmatZ <- sqrt(Z * weights) * xmat
            er2 <- as.matrix(eres2[margin])
            geefit <- geese.fit(xmat, yhat, id, zsca = er2, scale.fix = TRUE, corstr = corstr)
            ## }
        }
        beta <- geefit$beta
        if (control$trace) {
            cat("\n beta:\n")
            ## cat(beta)
            print(as.numeric(beta))
        }
        convStep = i
        if (max(abs(beta - betaprev)) <= control$abstol) break
        if (max(abs(beta - betaprev)/abs(beta)) <= control$reltol) break
    } ## end i for 1:maxiter
    alpha <- geefit$alpha
    ## if (res == TRUE && sum(margin == 1) != nobs) {
    gamma <- eres2
    ## }
    ## else {# (res == FALSE) {
    ##     gamma <- geefit$gamma
    ## }
    convergence <- if (i == control$maxiter) 1 else 0
    out <- list(beta = beta, alpha = alpha, gamma = gamma,
                convergence = convergence, convStep = convStep)
    return(out)
}


lss.eres <- function(e, delta, z = rep(1, length(e)))
{
  nobs <- length(e)
  ord <- order(e)
  ei <- e[ord]
  deltai <- delta[ord]
  zi <- z[ord]
  dummy <- 1:nobs
  repeats <- table(ei)
  Shat <- survfit(Surv(ei, deltai) ~ 1, weights = zi)$surv
  Shat <- rep(Shat, repeats)
  edif <- c(diff(ei), 0)  ## diff(ei) gives 1 less terms
  ehat <- rev(cumsum(rev(edif * Shat)))
  ehat2 <- rev(cumsum(rev(ei * edif * Shat)))
  ehat <- ehat/Shat + ei    ## +ei because there was a diff() in edif
  ehat2 <- 2 * ehat2/Shat + ei^2
  ehat[is.na(ehat)] = ei[is.na(ehat)]
  ehat2[is.na(ehat2)] = ei[is.na(ehat2)]^2
  ehat2[which(ehat2 < 0)] <- NaN
  eres <- ehat
  eres2 <- ehat2
  eres[dummy[ord]] <- ehat  ## puting it back to the original order
  eres2[dummy[ord]] <- ehat2
  return(list(eres, eres2))
}
