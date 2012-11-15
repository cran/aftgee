print.aftgee<- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("Call:\n")
  print(x$coefficients)
}

summary.aftgee <- function(object,...){
  z <- object
  if (class(z) != "aftgee"){
    stop("Most be aftgee class")
  }
  ans <- z["call"]
  TAB.ini <- NULL

  ## aftgee part
  est.gee <- z$coef.res
  se.gee <- sqrt(diag(z$var.res))
  est.temp.gee <- ifelse(se.gee == "NaN", "NaN", est.gee)
  z.val.gee <- as.numeric(est.temp.gee)/as.numeric(se.gee)
  TAB <- cbind(Estimate = est.gee, StdErr = se.gee, z.value = z.val.gee, p.value = 2 * pnorm(-abs(z.val.gee)))
  rownames(TAB) <- rownames(z$coefficients)
  ## binit part
  est.ini <- z$coef.init
  res <- list(call=object$call, coefficients=TAB, binit = z$binit, iniEst = z$iniEst, est.ini = z$coef.init)
  class(res) <- "summary.aftgee"
  res
}

summary.smoothrr <- function(object,...){
  z <- object
  if (class(z) != "smoothrr"){
    stop("Most be smoothrr class")
  }
  ans <- z["call"]
  est.srr <- z$beta
  p <- length(z$beta)
  var.meth <- z$var.meth[z$var.meth %in% c("MB", "ZLCF", "ZLMB", "sHCF", "sHMB", "ISCF", "ISMB", "js")]
  se.count <- length(var.meth)
  se.name <- match(var.meth, names(z$covmat))
  se.covmat <- z$covmat[se.name]
  TAB.srr <- NULL
  for (i in 1:se.count) {
      se.srr <- sqrt(diag(matrix(unlist(se.covmat[i]), p)))
      z.val.srr <- as.numeric(est.srr)/as.numeric(se.srr)
      temp.srr <- cbind(Estimate = est.srr, StdErr = se.srr, z.value = z.val.srr, p.value = 2 * pnorm(-abs(z.val.srr)))
      rownames(temp.srr) <- z$vari.name
      TAB.srr <- append(TAB.srr, list(temp.srr))
  }
  res <- list(call = object$call, coefficients = TAB.srr, var.name = names(z$covmat)[se.name])
  class(res) <- "summary.smoothrr"
  res
}


print.summary.aftgee <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("AFTGEE Estimator")
  cat("\n")
  printCoefmat(as.matrix(x$coefficients), P.values = TRUE, has.Pvalue = TRUE)
  ## if (is.numeric(x$binit) != TRUE) {
  ##     if (x$binit == "lm") {
  ##         cat("Initial Estimator from lm:")
  ##         cat("\n")
  ##         cat(format(round(as.numeric(x$est.ini), digits = 5), nsmall = 5))
  ##     }
  ##     if (x$binit == "srrgehan") {
  ##         cat("Initial Estimator from smoothrr with Gehan's weight:")
  ##         cat("\n")
  ##         cat(format(round(as.numeric(x$est.ini), digits = 5), nsmall = 5))
  ##     }
  ##     cat("\n")
  ##     cat("AFTGEE Estimator:")
  ##     cat("\n")
  ##     printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
  ## }
  ## if (is.numeric(x$binit) == TRUE) {
  ##     cat("Gehan Estimator:")
  ##     cat("\n")
  ##     cat(format(round(as.numeric(x$est.ini), digits = 5), nsmall = 5))
  ##     cat("\n")
  ##     cat("AFTGEE Estimator")
  ##     cat("\n")
  ##     printCoefmat(as.matrix(x$coefficients), P.values = TRUE, has.Pvalue = TRUE)
  ## }
}


print.summary.smoothrr <- function(x, ...){
  se.count <- length(x$var.name)
  cat("Call:\n")
  print(x$call)
  cat("\n")
  for (i in 1:se.count){
      cat("\n")
      cat("Variance Estimator:", as.character(x$var.name[i]))
      cat("\n")
      printCoefmat(as.data.frame(x$coefficients[i]), P.values = TRUE, has.Pvalue = TRUE)
  }
}
