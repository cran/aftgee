aftgee <- function(formula, data, subset, id = NULL, contrasts = NULL,
                   weights = NULL, margin = NULL, 
                   corstr="independence",
                   binit = "srrgehan", B = 100,
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
    id <- model.extract(m, id)
    mterms <- attr(m, "terms")
    weights <- model.extract(m, weights) 
    obj <- unclass(m[,1]) 
    if (class(m[[1]]) != "Surv" || ncol(obj) > 2)
        stop("aftsrr only supports Surv object with right censoring.", call. = FALSE)
    if (is.null(id)) id <- 1:nrow(obj)
    if (is.null(weights)) weights <- rep(1, nrow(obj))
    margin <- model.extract(m, margin)
    if (is.null(margin)) margin <- rep(1, nrow(obj))
    formula[[2]] <- NULL
    ## Create DF; the first 2 columns are from Surv with time and status
    ## time, status, id, weights, margin, x1, x2, ...
    if (formula == ~1) DF <- cbind(obj, zero = 0)
    else {
        DF <- cbind(obj, id, weights, margin, model.matrix(mterms, m, contrasts))
        yint <- (sum(colnames(DF) == "(Intercept)") > 0)
        if (sum(colnames(DF) == "(Intercept)") > 0)
            DF <- DF[,-which(colnames(DF) == "(Intercept)")]
    }
    DF <- as.data.frame(DF)
    out <- NULL
    if (sum(DF$status) == nrow(DF) & corstr %in% c("indep", "independence")) {
        warning("Response is uncensored and correlation structure is independence,
ordinary least squares is used", call. = FALSE)
        out <- lm(log(y[,1]) ~ x - 1)
        out$coef.init <- out$coef.res <- out$coefficients
        out$coefficients <- cbind(out$coefficients, out$coefficients)
        out$var.res <- vcov(out)
    }
    else {
        out <- aftgee.fit(DF = DF, corstr = corstr, B = B, binit = binit, control = control, yint = yint)
    } 
    out$y <- DF$time
    out$x <- DF[,-(1:5)]
    rownames(out$coefficients) <- names(out$coef.res) <- names(out$coef.init) <- colnames(model.matrix(mterms, m, contrasts))
    ## out$intercept <- (sum(x[,1]) == nrow(x))
    colnames(out$coefficients) <- c("binit", "AFTGEE")
    out$call <- scall
    class(out) <- "aftgee"
    out
}

aftgee.fit <- function(DF, corstr="independence",
                       B = 100, binit = "lm", yint = TRUE,
                       control = aftgee.control()) {
    x <- as.matrix(DF[,-(1:5)])
    id <- DF$id
    n <- length(unique(id))
    rm <- NULL
    rmName <- NULL
    firstBeta <- firstSd <- firstSdMat <- firstconvergence <- NULL
    clsize <- unlist(lapply(split(id, id), length))
    N <- sum(clsize)
    if (is.numeric(binit)) {
        if (length(binit) != ncol(x) + yint * 1)
            stop("binit value length does not match with numbers of covariates", call. = FALSE)
        firstBeta <- binit      
    }
    if (!(is.numeric(binit))) {
        if (!(binit %in% c("lm", "srrgehan"))) 
            stop("Invalid binit value method", call. = FALSE)
    }
    if (!(is.numeric(binit))) {
        if (binit == "lm") {
            if (yint) linfit <- summary(lm(log(DF$time) ~ x, subset = DF$time > 0))
            else linfit <- summary(lm(log(DF$time) ~ x - 1, subset = DF$time > 0))
            first <- list(beta = linfit$coef[,1], sd = linfit$coef[,2])
            firstBeta <- first$beta
            firstSd <- first$sd
            firstconvergence <- first$convergence
        }
        if (binit == "srrgehan") {
            engine.control <- control[names(control) %in% names(attr(getClass("gehan.is"), "slots"))]
            engine <- do.call("new", c(list(Class = "gehan.is"), engine.control))
            if (engine@b0 == 0) {
                engine@b0 <- as.numeric(coef(lm(DF$time ~ as.matrix(DF[,-(1:5)]))))[-1]
            }
            engine@sigma0 <- diag(length(engine@b0))
            first <- rankFit.gehan.is(DF[,-5], engine, NULL)
            firstBeta <- first$beta
            firstSdMat <- NA
            firstconvergence <- first$convergence
            if (yint) firstBeta <- c(mean(log(DF$time) - x %*% firstBeta), firstBeta) 
        }
    }
    if (yint) x <- cbind(1, x)
    binitValue <- list(beta = firstBeta, sd = firstSd, sdMat = firstSdMat)
    result <- aftgee.est(log(DF$time), x, DF$status, binitValue$beta, id, corstr,
                         rep(1, nrow(DF)), DF$margin, DF$weights, control)
    ## variance estimation
    sample <- zout <- NULL
    if (B > 0) {
        sample <- matrix(0, nrow = B, ncol = length(result$beta))
        for (i in 1:B){
            Z <- as.vector(rep(rexp(n,1), time = clsize))
            zout <- cbind(zout, Z)
            DF0 <- DF
            DF0$weights <- Z
            if (control$seIni) {
                boot.ini <- rankFit.gehan.is(DF0[,-5], engine, NULL)
                if (yint) boot.ini$beta <- c(mean(log(DF$time) - as.matrix(DF[,-(1:5)]) %*% boot.ini$beta), boot.ini$beta)
            } else boot.ini <- result
            sample[i,] <- aftgee.est(log(DF$time), x, DF$status, boot.ini$beta, id, corstr,
                                     Z, DF$margin, DF$weights, control)$beta
        }
        vhat <- var(sample)
    }
###############################################################################################################
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
                bhat = sample,
                zout = zout, 
                conv.step = result$convStep)
    class(fit) <- "aftgee.fit"
    fit
}

## aftgee.se; bootstrap or resampling, true value or aftsrr as initial value

aftgee.control <- function(maxiter = 50, reltol = 0.001, trace = FALSE, seIni = TRUE) {
    list(maxiter = maxiter, reltol = reltol, trace = trace, seIni = seIni)
}

aftgee.est <- function(y, x, delta, beta, id, corstr = "independence", Z = rep(1, length(y)),
                       margin = rep(1, length(id)), weights = rep(1, length(y)),
                       control = aftgee.control()) {
    xmat <- as.matrix(x) 
    nobs <- length(y)
    for (i in 1:control$maxiter) {
        betaprev <- beta
        eres <- NULL
        eres2 <- NULL
        if (sum(margin == 1) == nobs) {
            e <- y - xmat %*% beta
            eres <- eRes(e, delta = delta, z = Z * weights)
            yhat <- delta * y + (1 - delta) * (eres[[1]] + xmat %*% beta)
            yhatZ <- sqrt(Z) * yhat
            xmatZ <- sqrt(Z) * xmat
            geefit <- geese.fit(xmatZ, yhatZ, id, corstr = corstr, weights =  weights)
        }
        if (sum(margin == 1) != nobs) {
            e <- y - xmat %*% beta
            er1 <- NULL
            er2 <- NULL
            for (m in unique(margin)) {
                temp <- eRes(e[margin == m], delta[margin == m], Z[margin == m])
                temp[[2]] <- ifelse(delta[margin == m] == 1, e[margin == m]^2, temp[[2]])
                eres2[m] <- mean(temp[[2]], na.rm = TRUE)
                dum <- cumsum(ifelse(margin == m, 1, 0))
                er1temp <- temp[[1]][ifelse(margin == m, dum, NA)]
                er1 <- rbind(er1, er1temp)
            }
            er1 <- as.vector(er1)
            er1 <- er1[!is.na(er1)]
            yhat <- delta * y + (1 - delta) * (er1 + xmat %*% beta)
            yhatZ <- sqrt(Z * weights) * yhat
            xmatZ <- sqrt(Z * weights) * xmat
            er2 <- as.matrix(eres2[margin])
            ## geefit <- geese.fit(xmat, yhat, id, zsca = er2, scale.fix = TRUE, corstr = corstr, weights = Z * weights)
            geefit <- geese.fit(xmatZ, yhatZ, id, zsca = er2, scale.fix = TRUE, corstr = corstr)
        }
        beta <- geefit$beta
        if (control$trace) {
            cat("\n beta:\n")
            print(as.numeric(beta))
        }
        convStep <- i
        if (max(abs(beta - betaprev)/abs(beta)) <= control$reltol) break
    } ## end i for 1:maxiter
    iniBeta <- geefit$beta
    if ("(Intercept)" %in% colnames(x)) {
        ##    beta <- c(eRes(e = y - as.matrix(x[,-1]) %*% geefit$beta[-1], delta = delta)[[3]],
        ##            geefit$beta[-1])
        beta <- c(mean(yhat - as.matrix(x[,-1]) %*% geefit$beta[-1]), geefit$beta[-1])
    } else {
        beta <- geefit$beta
    }
    alpha <- geefit$alpha
    gamma <- eres2
    convergence <- ifelse(i == control$maxiter, 1, 0)
    out <- list(beta = beta, alpha = alpha, gamma = gamma, iniBeta = iniBeta,
                convergence = convergence, convStep = convStep)
    return(out)
}

eRes <- function(e, delta, z = rep(1, length(e))) {
    nobs <- length(e)
    ord <- order(e)
    ei <- e[ord]
    deltai <- delta[ord]
    zi <- z[ord]
    dummy <- 1:nobs
    tmp <- survfit(Surv(ei, deltai) ~ 1, weights = zi)
    Shat <- with(tmp, approx(time, surv, ei))$y
    edif <- c(diff(ei), 0)  ## diff(ei) gives 1 less terms
    ehat <- rev(cumsum(rev(edif * Shat)))
    inpt <- mean(ehat)
    ehat2 <- rev(cumsum(rev(ei * edif * Shat)))
    ehat <- ehat/Shat + ei    ## +ei because there was a diff() in edif
    ehat2 <- 2 * ehat2/Shat + ei^2
    ehat[is.na(ehat)] <- ei[is.na(ehat)]
    ehat2[is.na(ehat2)] <- ei[is.na(ehat2)]^2
    ehat2[which(ehat2 < 0)] <- NaN
    eres <- ehat
    eres2 <- ehat2
    eres[dummy[ord]] <- ehat  ## puting it back to the original order
    eres2[dummy[ord]] <- ehat2
    return(list(eres, eres2, inpt))
}

## ## Internal function for obtaining the se estimator for an aftgee object 
## aftgee.se <- function(DF, x, B, b0, yint, control = aftgee.control()) {
##     ## use the standard bootstrap when resampling = FALSE
##     clsize <- unlist(lapply(split(DF$id, DF$id), length))
##     n <- length(unique(DF$id))
##     bb <- NULL
##     if (control$parallel) {
##         cl <- makeCluster(stdErr@parCl)
##         clusterExport(cl = cl, varlist=c("DF", "x"), envir = environment())
##         if (control$resampling) {
##             out <- parSapply(cl, 1:B, function(x) {
##                 Z <- as.vector(rep(rexp(n, 1), time = clsize))
##                 DF$weights <- Z
##                 boot.ini <- b0
                
##                 }


##         }
##     }
## }
