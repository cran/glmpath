"step.length.cox" <- function(corrector, x, d, rslist, wlist, min.lambda, max.arclength, add.newvars, backshoot, approx.Gram, h0=NULL, eps=.Machine$double.eps)
  {
    active <- corrector$active
    lambda <- corrector$lambda - min.lambda
    C <- corrector$corr
    b <- corrector$b[active]
    eta <- corrector$eta
    wsum <- corrector$wsum
    n <- length(d)
    p1 <- length(active)
    p2 <- ncol(x)
    A <- rep(0, n)
    if (!approx.Gram) AA <- matrix(0, n, n)
    rset <- w <- NULL
    for (i in 1:sum(d==1)) {
      if (!is.null(rslist[[i]])) {
        rset0 <- rset
        rset <- c(rset0, rslist[[i]])
      }
      w <- c(rep(1,length(rset0)), wlist[[i]]) * eta[rset]
      w1 <- w/wsum[i]
      A[rset] <- A[rset] + w1 - w1^2
      if (!approx.Gram) {
        k <- length(rset)
        AA[1:k,1:k] <- AA[1:k,1:k] - outer(w1[rset], w1[rset])
      }
    }
    if (approx.Gram) dhdb <- t(x[ ,active,drop=FALSE] * A) %*% x
    else {
      diag(AA) <- A
      dhdb <- t(x[ ,active,drop=FALSE]) %*% AA %*% x
    }
    db <- solve(dhdb[ ,active,drop=FALSE]) %*% sign(C[active])
    newa <- NULL
    if (!backshoot) {
      Cmax <- max(abs(C))
      inactive <- seq(C)[-active]
      ninact <- length(inactive)
      a <- drop(t(db) %*% dhdb[ ,inactive,drop=FALSE])
      gam <- c((Cmax - C[inactive])/(1 - a), (Cmax + C[inactive])/(1 + a))
      ha <- min(gam[gam > eps], lambda)
      hd <- -b/db
      h <- min(hd[hd > eps], ha)
      if (h==ha & h < lambda) {
        ii <- which(gam > eps)
        nii <- length(ii)
        oo <- order(gam[ii])
        oo <- oo[1:min(nii, add.newvars)]
        oo <- ii[oo]
        oo <- unique(ifelse(oo <= ninact, oo, oo-ninact))
        newa <- inactive[oo]
      }
      h <- min(h, max.arclength/sum(abs(db)))
    }
    else {
      hd <- b/db
      ii <- hd > eps & hd < -h0
      if (any(ii)) h <- -max(hd[ii])
      else h = 0
    }
    list(h = -h, db = db, newa = newa)
  }

"predictor.cox" <- function(b, step)
  {
    b - step$db * step$h
  }

"corrector.cox" <- function(x, y, d, rslist, wlist, rept, method, active, tmpa, lambda, b0, a0, M0, bshoot.threshold, relax.lambda, trace, eps = .Machine$double.eps)
  {
    param <- c(b0[tmpa], M0, a0, -b0[tmpa]-a0, 0, b0[tmpa]-a0)
    xa <- x[ , tmpa, drop=FALSE]
    nobs <- nrow(xa)
    p <- ncol(xa)
    xstate <- rep(2, 4*p+2)
    xstate[param==0] <- 0
    lenz <- 10 + (p+5)*nobs 
    zsmall <- rep(0, lenz)
    zsmall[1:3] <- c(nobs, lambda, method)
    zsmall[10 + seq((p+3)*nobs)] <- c(as.vector(xa), y, d, rept)
    sol <- .Fortran("coxsolution",
                    k = as.integer(p),
                    n = as.integer(2*p+1),
                    nb = as.integer(4*p+2),
                    ne = as.integer(5*p+1),
                    hs = as.integer(xstate),
                    xn = as.double(param),
                    zsmall = as.double(zsmall),
                    lenz = as.integer(lenz),
                    inform = as.integer(0))
    b0[tmpa] <- sol$xn[1:p]
    M0 <- sol$xn[p+1]
    a0 <- sol$xn[(p+2):(2*p+1)]
    if (sol$inform != 0) cat("\nconvergence warning in corrector step\n")
    ii <- 10+(p+3)*nobs
    eta <- sol$zsmall[ii+c(1:nobs)]
    ii <- 10+(p+4)*nobs
    wsum <- sol$zsmall[ii+c(1:nobs)][d==1]
    rset <- NULL
    a <- d==1
    for (i in 1:sum(a)) {
      if (!is.null(rslist[[i]])) {
        rset0 <- rset
        rset <- c(rset0, rslist[[i]])
      }
      w <- c(rep(1,length(rset0)), wlist[[i]]) * eta[rset]
      a[rset] <- a[rset] - w/wsum[i]
    }
    corr <- apply(x*a, 2, sum)
    i <- which(abs(corr) >= lambda*(1-relax.lambda))
    newa <- i[!i%in%tmpa]
    newactive <- i[!i%in%active]
    i <- which(abs(b0[active]) < eps)
    inactive <- active[i]
    active <- active[!active%in%inactive]
    active <- c(active, newactive)
    b0[-active] <- 0
    df <- length(active) - length(newactive)
    backshoot <- ifelse(any(abs(b0[newactive]) > bshoot.threshold), TRUE, FALSE)
    list(eta=eta, wsum=wsum, b=b0, a=a0, M=M0, lp=sol$zsmall[4], active=active, newactive=newactive, newa=newa, inactive=inactive, corr=corr, lambda=lambda, df=df, backshoot=backshoot)
  }

"coxpath" <- function(data, method = c("breslow", "efron"), max.steps = NULL, max.norm = 100*ncol(data$x), min.lambda = 0, max.arclength = Inf, add.newvars = 1, bshoot.threshold = 0.1, relax.lambda = 1e-7, approx.Gram = FALSE, eps = .Machine$double.eps, trace = FALSE)
  {
    method <- match.arg(method)
    method <- switch(method, breslow = 1, efron = 2)
    x <- data$x
    time <- data$time
    status <- data$status
    n <- nrow(x)
    m <- ncol(x)
    o <- order(status)
    oo <- o[order(time[o], decreasing=T)]
    x <- x[oo, ]
    time <- time[oo]
    status <- status[oo]
    complete <- which(status==1)
    nnc <- length(complete)
    rept <- rep(0, n)
    for (i in complete) rept[i] <- sum(time[i:n]==time[i] & status[i:n]==1)
    rslist <- wlist <- vector("list", length=nnc)
    for (i in 1:nnc) {
      if (i==1) {
        ii <- time >= time[complete[i]]
        rslist[[i]] <- which(ii)
      }
      else if (rept[complete[i]] >= rept[complete[i]-1]) {
        ii <- (time >= time[complete[i]]) & (time < time[complete[i-1]])
        rslist[[i]] <- which(ii)
      }
      wlist[[i]] <- rep(1, sum(ii))
      if (method == 2) {
        if (rept[complete[i]] > 0) {
          tie <- time[ii]==time[complete[i]] & status[ii]==1
          di <- max(rept[ii][tie])
          wlist[[i]][tie] <- wlist[[i]][tie] - (di - rept[complete[i]])/di
        }
      }  
    }
    if (m >= n) min.lambda <- ifelse(min.lambda == 0, 1e-3, min.lambda)
    if (is.null(max.steps)) max.steps <- 10*min(m, n)
    one <- rep(1, n)
    meanx <- drop(one %*% x)/n
    x <- scale(x, meanx, FALSE)
    sdx <- sqrt(drop(one %*% (x^2))/(n-1))
    ignore <- sdx < eps
    if (any(ignore)) sdx[ignore] <- 1
    x <- scale(x, FALSE, sdx)
    if (is.null(dimnames(x)[[2]])) xnames <- paste("x",seq(m),sep="")
    else xnames <- dimnames(x)[[2]]
    lam.vec <- step.len <- rep(0, max.steps)
    bmat.pred <- bmat.corr <- cmat <- matrix(0, nrow=max.steps, ncol=m)
    lp <- df <- new.df <- rep(0, max.steps)
    actions <- vector("list", length=max.steps)
    backshoot <- FALSE
    b <- rep(0, m)
    rset <- NULL
    wsum <- rep(0, nnc)
    a <- status==1
    for (i in 1:nnc) {
      if (!is.null(rslist[[i]])) {
        rset0 <- rset
        rset <- c(rset0, rslist[[i]])
      }
      w <- c(rep(1,length(rset0)), wlist[[i]])
      wsum[i] <- sum(w)
      a[rset] <- a[rset] - w/wsum[i]
    }
    init.corr <- apply(x*a, 2, sum)
    k <- 1
    lp[k] <- -sum(log(wsum))
    cmat[k, ] <- init.corr
    lam.vec[k] <- lambda <- max(abs(init.corr))
    actions[[k]] <- active <- which(abs(init.corr) == lambda)
    names(actions[[k]]) <- xnames[abs(active)]
    if (trace) {
      cat(paste("Lambda=",lambda,"lets the first factor in.\nStep",k,":"))
      cat(paste("\t",xnames[active],"added"))
    }
    corrector <- list(eta=rep(1, n), wsum=wsum, b=sign(init.corr), active=active, corr=init.corr, lambda=lambda) 
    while(TRUE) {
      if (!backshoot) {
        k <- k + 1
        if (trace) cat(paste("\nStep",k,":"))
        step <- step.length.cox(corrector, x, status, rslist, wlist, min.lambda, max.arclength, add.newvars, backshoot, approx.Gram)
        b[active] <- predictor.cox(b[active], step)
        bmat.pred[k, ] <- b
        step.len[k-1] <- h <- step$h
        lam.vec[k] <- lambda <- lambda + h
        tmpa <- c(active, step$newa)
        a <- abs(b[tmpa])
        M <- sum(a)
      }
      else {
        if (trace) cat(paste("\nStep",k,":"))
        step <- step.length.cox(corrector, x, status, rslist, wlist, min.lambda, Inf, add.newvars, backshoot, approx.Gram, h)
        step.len[k-1] <- h + step$h
        h <- step$h
        lam.vec[k] <- lambda <- lambda + h
        a <- abs(b[tmpa])
        M <- sum(a)
      }
      corrector <- corrector.cox(x, time, status, rslist, wlist, rept, method, active, tmpa, lambda, b, a, M, bshoot.threshold, relax.lambda, trace)
      newa <- corrector$newa
      while(length(newa) > 0) {
        if (trace) cat(paste("\nRepeating step",k,":"))
        tmpa <- c(tmpa, newa)
        a <- abs(b[tmpa])
        M <- sum(a)
        corrector <- corrector.cox(x, time, status, rslist, wlist, rept, method, active, tmpa, lambda, b, a, M, bshoot.threshold, relax.lambda, trace)
        newa <- corrector$newa 
      }
      newaction <- c(corrector$newactive, -corrector$inactive)
      if (length(newaction) > 0 & length(corrector$active) <= n) {
        if (corrector$backshoot & !backshoot) {
          if (trace) cat("\nOvershooting occurred: increasing lambda again")
          backshoot <- TRUE
        }
        else {
          active <- corrector$active
          b <- corrector$b
          actions[[k]] <- newaction
          names(actions[[k]]) <- xnames[abs(newaction)]
          new.df[k] <- corrector$df
          if (trace) {
            na <- newaction[newaction > 0]
            ina <- -newaction[newaction < 0]
            if (length(na) > 0) cat(paste("\t",xnames[na],"added"))
            if (length(ina) > 0) cat(paste("\t",xnames[ina],"dropped"))
          }
          backshoot <- FALSE
        }
      }
      else if (length(corrector$active) <= n) {
        active <- corrector$active
        b <- corrector$b
        backshoot <- FALSE
      }
      if (!backshoot) {
        bmat.corr[k, ] <- b 
        cmat[k, ] <- corrector$corr
        lp[k] <- corrector$lp
        df[k] <- corrector$df
        if (lambda <= min.lambda | k == max.steps | length(corrector$active) > n | corrector$M >= max.norm) {
          if (trace & lambda <= min.lambda) cat("\nLambda=",min.lambda,"\n")
          else if (trace & k == max.steps) cat(paste("\nMaximum steps (",max.steps,") taken.\n"))                     
          else if (length(corrector$active) > n) {
            k <- k - 1
            if (trace) cat("\nNumber of active variables has reached its maximum.\n")
          }
          else if (trace) cat("\n|beta| >=",max.norm,"\n")
          break
        }
      }
    }
    bmat.pred <- bmat.pred[1:k, ]
    bmat.corr <- bmat.corr[1:k, ]
    cmat <- cmat[1:k, ]
    bmat.pred <- scale(bmat.pred, FALSE, sdx)
    bmat.corr <- scale(bmat.corr, FALSE, sdx)
    dimnames(bmat.pred) <- dimnames(bmat.corr) <- dimnames(cmat) <- list(seq(k), xnames)
    df <- df[1:k]
    lp <- lp[1:k]
    aic <- -2*lp + 2*df
    bic <- -2*lp + log(n)*df
    object <- list(lambda=lam.vec[1:k], step.length=abs(step.len[1:(k-1)]), corr = cmat, new.df = new.df[1:k], df = df, loglik = lp, aic = aic, bic = bic, b.predictor = bmat.pred, b.corrector = bmat.corr, actions=actions[1:k], meanx=meanx, sdx=sdx, xnames=xnames, method=ifelse(method==1, "breslow", "efron"))
    class(object) <- "coxpath"
    object
  }

"plot.coxpath" <- function(x, xvar = c("norm", "lambda", "step"), type = c("coefficients", "aic", "bic"), xlimit = NULL, predictor = FALSE, omit.zero = TRUE, breaks = TRUE, mar = NULL, main = NULL, eps = .Machine$double.eps, ...)
  {
    object <- x
    lam <- object$lambda
    xvar <- match.arg(xvar)
    type <- match.arg(type)
    coef.pred <- scale(object$b.predictor, FALSE, 1/object$sdx)
    coef.corr <- scale(object$b.corrector, FALSE, 1/object$sdx)
    xnames <- object$xnames
    if (omit.zero) {
      c1 <- drop(rep(1, nrow(coef.corr)) %*% abs(coef.corr))
      nonzero <- c1 > eps
      xnames <- xnames[nonzero]
      coef.pred <- coef.pred[ ,nonzero]
      coef.corr <- coef.corr[ ,nonzero]
    }
    m <- ncol(coef.pred)
    k <- nrow(coef.pred)
    s <- switch(xvar, norm = apply(abs(coef.corr), 1, sum),
      lambda = lam,
      step = seq(k))
    if (xvar != "lambda") {
      if (is.null(xlimit)) xlimit <- max(s)
      else if (xlimit <= min(s)) stop("Increase xlimit.")
      xi <- s <= xlimit
    }
    else {
      if (is.null(xlimit)) xlimit <- min(s)
      else if (xlimit >= max(s)) stop("Decrease xlimit.")
      xi <- s >= xlimit
    }
    k <- max(which(xi))
    new <- c(T, object$new.df[-1] > 0) & xi
    xname <- switch(xvar, norm = "|beta|", lambda = "lambda", step = "step")
    if (!is.null(mar)) par(mar=mar)
    if (type == "aic") {
      aic <- object$aic[xi]
      plot(s[xi], aic, xlab = xname, ylab = "AIC", type = "b", pch = 16, cex = 0.3, ...)
      if (is.null(main)) title("AIC", line=2.5)
      else title(main, line=2.5)
    }
    else if (type == "bic") {
      bic <- object$bic[xi]
      plot(s[xi], bic, xlab = xname, ylab = "BIC", type = "b", pch = 16, cex = 0.3, ...)
      if (is.null(main)) title("BIC", line=2.5)
      else title(main, line=2.5)
    }
    else {
      matplot(s[xi], coef.corr[xi, ], xlab = xname, ..., type = "b", pch = "*", ylab = "Standardized coefficients", lty=1)
      if (is.null(main)) title("Coefficient path", line=2.5)
      else title(main, line=2.5)
      abline(h = 0, lty = 3)
      axis(4, at = coef.corr[k, ], labels = xnames, cex = 0.8, adj = 0, las = 1)
      if (predictor) {
        for (i in 1:m) segments(s[xi][-k], coef.corr[xi, ][-k,i], s[xi][-1], coef.pred[xi, ][-1,i], lty=2, col=i)
      }
    }
    if (breaks) {
      axis(3, at = s[new], labels = object$new.df[new], cex = 0.8)
      abline(v = s[new])
    }
  }

"predict.coxpath" <- function(object, data, s, type = c("coefficients", "loglik", "lp", "risk", "coxph"), mode = c("step", "norm.fraction", "norm", "lambda.fraction", "lambda"), eps = .Machine$double.eps, ...)
  {
    mode <- match.arg(mode)
    type <- match.arg(type)
    if (missing(data) & type!="coefficient") {
        warning("No data argument; type switched to coefficients")
        type <- "coefficients"
    }    
    if (length(s) > 1 & type=="coxph") {
        warning("Length(s) > 1. Only the first element is used.")
        s <- s[1]
    }
    b <- object$b.corrector
    std.b <- scale(b, FALSE, 1/object$sdx)
    k <- nrow(b)
    steps <- seq(k)
    if (missing(s)) {
        s <- steps
        mode <- "step"
    }
    sb <- switch(mode, step = {
        if (any(s < 1) | any(s > k)) 
            stop("Argument s out of range")
        steps
    }, norm.fraction = {
        if (any(s > 1) | any(s < 0)) 
            stop("Argument s out of range")
        bnorm <- apply(abs(std.b), 1, sum)
        bnorm / bnorm[k]
    }, norm = {
        bnorm <- apply(abs(std.b), 1, sum)
        if (any(s > bnorm[k]) | any(s < bnorm[1])) 
            stop("Argument s out of range")
        bnorm
    }, lambda.fraction = {
        if (any(s > 1) | any(s < 0))
            step("Argument s out of range")
        lam <- object$lambda
        lam[lam < eps] <- eps
        lam <- log(lam)
        (lam - min(lam)) / (max(lam) - min(lam))
    }, lambda = {
        lam <- object$lambda
        if (any(s > lam[1]) | any(s < lam[k]))
          stop("Argument s out of range")
        lam
    })
    sfrac <- (s - sb[1])/(sb[k] - sb[1])
    sb <- (sb - sb[1])/(sb[k] - sb[1])
    usb <- unique(sb)
    useq <- match(usb, sb)
    sb <- sb[useq]
    b <- b[useq, ]
    coord <- approx(sb, seq(sb), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    newb <- ((sb[right] - sfrac) * b[left, , drop = FALSE] + 
        (sfrac - sb[left]) * b[right, , drop = FALSE])/(sb[right] - sb[left])
    newb[left == right, ] <- b[left[left == right], ]    
    coef <- newb
    if (type=="coefficients") {
      fit <- coef
      dimnames(fit) <- list(s, object$xnames)
    }
    else if (type=="loglik") {
      fit <- logplik(data$x, data$time, data$status, t(coef), object$method)
      names(fit) <- s
    }
    else if (type=="lp" | type=="risk") {
      b0 <- coef %*% object$meanx 
      fit <- scale(data$x %*% t(coef), b0, F)
      if (type=="risk") fit <- exp(fit)
      dimnames(fit) <- list(seq(nrow(data$x)), s)      
    }
    else {
      coef <- drop(coef)
      active <- abs(coef) > eps
      coef <- coef[active]
      x <- data$x[ ,active,drop=FALSE]
      time <- data$time
      status <- data$status
      fit <- coxph(Surv(time, status) ~ x, method=object$method)
      junk <- logplik(x, time, status, coef, object$method, TRUE)
      w <- junk$w
      dmat <- junk$dmat
      oo <- junk$oo
      a <- sum(active)
      info <- matrix(0, a, a)
      for (i in 1:sum(status==1)) {
        ind <- dmat[ ,i] > 0
        xr <- x[oo[ind], ,drop=FALSE]
        wr <- w[ind, i]
        v1 <- xr * wr
        v2 <- apply(v1, 2, sum)
        info <- info + t(xr) %*% v1 - outer(v2, v2)
      }
      fit$coefficients <- coef
      fit$var <- solve(info)
      fit$loglik <- c(fit$loglik[1], junk$loglik)
      fit$iter <- fit$residuals <- NULL
      fit$linear.predictors <- junk$eta - sum(coef*object$meanx)
      fit$method <- object$method
      fit$assign <- seq(a)
      fit$wald.test <- sum(coef*(info %*% coef))
    }    
    attr(fit, "s") <- s
    attr(fit, "fraction") <- sfrac
    attr(fit, "mode") <- mode
    return(fit)
  }

"logplik" <- function(x, time, status, b, method = c("breslow", "efron"), return.all=FALSE)
  {
    method <- match.arg(method)
    n <- length(time)
    o <- order(status, decreasing=T)
    oo <- o[order(time[o])]
    time <- time[oo]
    status <- status[oo]
    rept <- rep(0, n)
    for (i in 1:n) rept[i] <- sum(time[i:n]==time[i] & status[i:n]==1)
    complete <- which(status==1)
    nnc <- length(complete)
    dmat <- matrix(0, n, nnc)
    for (i in 1:nnc) {
      dmat[time >= time[complete[i]], i] <- 1
      if (method=="efron") {
        if (rept[complete[i]] > 0) {
          tie <- time==time[complete[i]] & status==1
          di <- max(rept[tie])
          dmat[tie, i] <- dmat[tie, i] - (di - rept[complete[i]])/di
        }
      }
    }
    eta <- x %*% b
    eeta <- exp(eta)
    k <- ncol(eta)
    loglik <- rep(0, k)
    for (i in 1:k) {
      w <- dmat * eeta[oo, i]
      wsum <- apply(w, 2, sum)
      loglik[i] <- sum(eta[oo, i][status==1]) - sum(log(wsum))
    }
    if (return.all)
      return(list(loglik=loglik, w=scale(w, F, wsum), eta=eta, dmat=dmat, oo=oo))
    else return(loglik)    
  }

"cv.coxpath" <- function (data, method = c("breslow", "efron"), nfold = 5, fraction = seq(from=0, to=1, length=100), mode = c("norm","lambda"), plot.it = TRUE, se = TRUE, ...)
  {
    method <- match.arg(method)
    mode <- match.arg(mode)
    x <- data$x
    time <- data$time
    status <- data$status
    n <- length(time)
    all.folds <- split(sample(seq(n)), rep(1:nfold, length=n))
    lpmat <- matrix(0, length(fraction), nfold)
    for (i in seq(nfold)) {
        omit <- all.folds[[i]]
        trdata <- list(x=x[-omit, ], time=time[-omit], status=status[-omit])
        tsdata <- list(x=x[omit, ], time=time[omit], status=status[omit])
        fit <- coxpath(trdata, method, ...)
        pred <- switch(mode, norm = {
          predict(fit, tsdata, fraction, "loglik", "norm.fraction")
        }, lambda = {
          predict(fit, tsdata, fraction, "loglik", "lambda.fraction")
        })
        if (length(omit) == 1) pred <- matrix(pred, nrow = 1)
        lpmat[ ,i] <- pred
        cat("CV Fold", i, "\n")
    }
    cv.error <- apply(lpmat, 1, mean)
    cv.se <- sqrt(apply(lpmat, 1, var)/nfold)
    object <- list(fraction = fraction, cv.error = cv.error, cv.se = cv.se, folds=all.folds)
    if (plot.it) {
      plot(fraction, cv.error, type="l", ylim=range(c(cv.error-cv.se, cv.error+cv.se)), xlab=switch(mode, norm="Norm fraction", lambda="log(lambda) fraction"), ylab="Log-partial-likelihood", main="Cross-validated log-partial-likelihood")
      if (se) segments(fraction, cv.error-cv.se, fraction, cv.error+cv.se)
    }
    invisible(object)
}

"print.coxpath" <- function(x, ...)
  {
    actions <- x$actions
    xn <- x$xnames
    k <- length(actions)
    for (i in 1:k) {
      if (length(actions[[i]]) > 0) cat("Step",i,":",xn[abs(actions[[i]])],"\n")
    }
  }
