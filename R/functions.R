"step.length" <- function(corrector, min.lambda, max.arclength, add.newvars, backshoot, h0=NULL, eps=.Machine$double.eps)
  {
    active <- corrector$active
    lambda <- corrector$lambda - min.lambda
    C <- corrector$corr
    b <- corrector$b[active]
    xw <- corrector$xw
    xwx <- t(xw) %*% xw[ ,active]
    xwx.in <- solve(xwx[active, ])
    db <- drop(xwx.in %*% c(0, sign(C[active[-1]])))
    newa <- NULL
    if (!backshoot) {
      Cmax <- max(abs(C))
      inactive <- seq(C)[-active]
      ninact <- length(inactive)
      a <- drop(xwx[inactive, ] %*% db)
      gam <- c((Cmax - C[inactive])/(1 - a), (Cmax + C[inactive])/(1 + a))
      ha <- min(gam[gam > eps], lambda)
      hd <- -b[-1]/drop(db[-1])
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

"predictor1" <- function(b, step)
  {
    b - step$db * step$h
  }

"corrector1" <- function(x, y, family, weight, k, active, tmpa, lambda, b0, a0, M0, bshoot.threshold, relax.lambda, trace, no.iter = FALSE, eps = .Machine$double.eps)
  {
    variance <- family$variance
    linkinv <- family$linkinv
    mu.eta <- family$mu.eta
    dev.resids <- family$dev.resids
    k <- length(active) - 1
    kk <- length(tmpa) - 1
    param <- c(b0[tmpa], M0, a0, 0, -b0[tmpa[-1]]-a0, 0, b0[tmpa[-1]]-a0)
    xa <- x[ ,tmpa]
    nobs <- nrow(xa)
    p <- ncol(xa)
    xstate <- rep(2, 4*kk+4)
    xstate[param==0] <- 0
    if (!no.iter) {
      dstr <- switch(family$family, binomial=1, poisson=2)
      lenz <- 10+(p+2)*nobs
      zsmall <- rep(0, lenz)
      zsmall[1:2] <- c(lambda, dstr)
      zsmall[10 + seq((p+2)*nobs)] <- c(as.vector(xa), y, weight)
      sol <- .Fortran("solution",
                      k = as.integer(kk),
                      nobs = as.integer(nobs),
                      n = as.integer(2*kk+2),
                      nb = as.integer(4*kk+4),
                      ne = as.integer(5*kk+2),
                      hs = as.integer(xstate),
                      xn = as.double(param),
                      zsmall = as.double(zsmall),
                      lenz = as.integer(lenz),
                      inform = as.integer(0))
      b0[tmpa] <- sol$xn[1:p]
      M0 <- sol$xn[kk+2]
      a0 <- sol$xn[(kk+3):(2*kk+2)]
      if (sol$inform != 0) cat("\nconvergence warning in corrector step\n")
    }
    eta <- drop(x%*%b0)
    mu <- linkinv(eta)
    mu.eta.eta <- mu.eta(eta)
    w <- (weight * mu.eta.eta^2/variance(mu))^0.5
    z <- (y-mu)/mu.eta.eta
    xw <- x * w
    wz <- w * z
    corr <- drop(t(xw) %*% wz)
    i <- which(abs(corr) >= lambda*(1-relax.lambda))
    newa <- i[!i%in%tmpa]
    newactive <- i[!i%in%active]
    i <- which(abs(b0[active[-1]]) < eps)
    inactive <- active[-1][i]
    active <- active[!active%in%inactive]
    active <- c(active, newactive)
    b0[-active] <- 0
    df <- length(active) - length(newactive)
    aic <- sum(dev.resids(y, mu, weight)) + 2*df
    bic <- aic + (log(nobs)-2)*df
    backshoot <- ifelse(any(abs(b0[newactive]) > bshoot.threshold), TRUE, FALSE)
    list(b=b0, a=a0, M=M0, active=active, newactive=newactive, newa=newa, inactive=inactive, corr=corr, lambda=lambda, xw=xw, df=df, aic=aic, bic=bic, backshoot=backshoot)
  }

"glmpath" <- function(x, y, family = binomial, weight = rep(1,length(y)), max.steps = NULL, max.norm = 100*ncol(x), min.lambda = 0, max.arclength = Inf, standardize = TRUE, add.newvars = 1, bshoot.threshold = 0.5, relax.lambda = 1e-9, eps = .Machine$double.eps, trace = FALSE)
  {
    if (is.character(family)) 
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
        family <- family()
    else if (family$family=="gaussian") family <- gaussian()
    else if (family$family=="binomial") family <- binomial()
    else if (family$family=="poisson") family <- poisson()    
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    n <- nrow(x)
    m <- ncol(x)
    if (m >= n) min.lambda <- max(min.lambda, 1e-6)
    if (is.null(max.steps)) max.steps <- 10*min(m, n)
    no.iter <- FALSE
    if (family$family=="gaussian" & family$link=="identity") no.iter <- TRUE
    if (family$family=="binomial") {
      uy <- unique(y)
      if (all(uy==c(1,-1)) | all(uy==c(-1,1))) {
        y[y==-1] <- 0
        cat("y=-1 converted to y=0 \n")
      }
    }
    one <- rep(1, n)
    meanx <- drop(one %*% x)/n
    x <- scale(x, meanx, FALSE)
    sdx <- sqrt(drop(one %*% (x^2))/(n-1))
    ignore <- sdx < eps
    if (any(ignore)) sdx[ignore] <- 1
    if (!standardize) sdx <- rep(1, m)
    x <- scale(x, FALSE, sdx)
    x <- cbind(1, x)
    if (is.null(dimnames(x)[[2]])) xnames <- paste("x",seq(m+1)-1,sep="")
    else xnames <- c("Intercept", dimnames(x)[[2]][-1])
    lam.vec <- step.len <- rep(0, max.steps)
    bmat.pred <- bmat.corr <- cmat <- matrix(0, nrow=max.steps, ncol=(m+1))
    new.df <- aic <- bic <- rep(0, max.steps)
    actions <- vector("list", length=max.steps)
    backshoot <- fullactive <- FALSE
    b <- rep(0, (m+1))
    init.mu <- mean(y)
    b0 <- family$linkfun(init.mu)
    if (trace) cat(paste("Initial intercept is",b0,"\n"))
    init.corr <- t(x)%*%((y-init.mu)*weight)
    lam.vec[1] <- lambda <- max(abs(init.corr[-1]))
    c1 <- which.max(abs(init.corr[-1])) + 1
    active <- c(1, c1)
    k <- 1
    if (trace) cat(paste("Lambda=",lambda,"lets the first factor in","\nStep",k,":","\t",xnames[c1],"added"))
    b[1] <- b0
    bmat.pred[k, ] <- bmat.corr[k, ] <- b
    actions[[k]] <- active
    eta <- drop(x%*%b)
    mu <- family$linkinv(eta)
    mu.eta.eta <- family$mu.eta(eta)
    w <- (weight * mu.eta.eta^2/family$variance(mu))^0.5
    corrector <- list(corr=init.corr, lambda=lambda, active=active, xw=x*w, b=sign(init.corr)) 
    new.df[k] <- 1
    while(TRUE) {
      if (!backshoot) {
        k <- k + 1
        if (trace) cat(paste("\nStep",k,":"))
        step <- step.length(corrector, min.lambda, max.arclength, add.newvars, backshoot)
        predictor <- predictor1(b[active], step)
        b[active] <- predictor
        bmat.pred[k, ] <- b
        step.len[k-1] <- h <- step$h
        lam.vec[k] <- lambda <- lambda + h
        tmpa <- c(active, step$newa)
        a <- abs(b[tmpa[-1]])
        M <- sum(a)
      }
      else {
        if (trace) cat(paste("\nStep",k,":"))
        step <- step.length(corrector, min.lambda, Inf, add.newvars, backshoot, h)
        step.len[k-1] <- h + step$h
        h <- step$h
        lam.vec[k] <- lambda <- lambda + h
        a <- abs(b[tmpa[-1]])
        M <- sum(a)
      }
      corrector <- corrector1(x, y, family, weight, k, active, tmpa, lambda, b, a, M, bshoot.threshold, relax.lambda, trace, no.iter)
      newa <- corrector$newa
      while(length(newa) > 0) {
        if (trace) cat(paste("\nRepeating step",k,":"))
        tmpa <- c(tmpa, newa)
        a <- abs(b[tmpa[-1]])
        M <- sum(a)
        corrector <- corrector1(x, y, family, weight, k, active, tmpa, lambda, b, a, M, bshoot.threshold, relax.lambda, trace, no.iter)
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
          names(actions[[k]]) <- xnames[abs(actions[[k]])]
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
        aic[k] <- corrector$aic
        bic[k] <- corrector$bic
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
    bmat.pred[ ,-1] <- scale(bmat.pred[ ,-1], FALSE, sdx)
    bmat.corr[ ,-1] <- scale(bmat.corr[ ,-1], FALSE, sdx)
    bmat.corr[1,-1] <- 0
    dimnames(bmat.pred) <- dimnames(bmat.corr) <- dimnames(cmat) <- list(seq(k), xnames)
    object <- list(lambda=lam.vec[1:k], step.length=abs(step.len[1:(k-1)]), new.df = new.df[1:k], corr = cmat, aic = aic[1:k], bic = bic[1:k], b.predictor = bmat.pred, b.corrector = bmat.corr, actions=actions[1:k], meanx=meanx, sdx=sdx, xnames=xnames, family=family, weight=weight)
    class(object) <- "glmpath"
    object
  }

"plot.glmpath" <- function(x, xvar = c("norm", "lambda", "step"), type = c("coefficients", "aic", "bic"), xlimit = NULL, predictor = FALSE, omit.zero = TRUE, breaks = TRUE, mar = NULL, eps = .Machine$double.eps, main = NULL, ...)
  {
    object <- x
    lam <- object$lambda
    xvar <- match.arg(xvar)
    type <- match.arg(type)
    coef.pred <- scale(object$b.predictor[ ,-1], FALSE, 1/object$sdx)
    coef.corr <- scale(object$b.corrector[ ,-1], FALSE, 1/object$sdx)
    xnames <- object$xnames[-1]
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
      aic <- object$aic[xi][-1]
      plot(s[xi][-1], aic, xlab = xname, ylab = "AIC", type = "b", pch = 16, cex = 0.3, ...)
      if (is.null(main)) title("AIC", line=2.5)
      else title(main, line=2.5)
    }
    else if (type == "bic") {
      bic <- object$bic[xi][-1]
      plot(s[xi][-1], bic, xlab = xname, ylab = "BIC", type = "b", pch = 16, cex = 0.3, ...)
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

"predict.glmpath" <- function(object, newx, s, type = c("link", "response", "coefficients"), mode = c("step", "norm.fraction", "norm", "lambda.fraction", "lambda"), eps = .Machine$double.eps, ...)
  {
    mode <- match.arg(mode)
    type <- match.arg(type)
    if (missing(newx) & type != "coefficients") {
        warning("Type=fit with no newx argument; type switched to coefficients")
        type <- "coefficients"
    }
    b <- object$b.corrector
    std.b <- scale(b[ ,-1], FALSE, 1/object$sdx)
    k <- nrow(b)
    steps <- seq(k)
    if (missing(s)) {
        s <- steps
        mode <- "step"
    }
    sb <- switch(mode, step = {
        if (any(s < 0) | any(s > k)) 
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
    if (type!="coefficients") {
      if (is.vector(newx)) newx <- matrix(newx, nrow=1)
      fit <- cbind(1,scale(newx, object$meanx, FALSE)) %*% t(newb)      
      if (type=="response") fit <- object$family$linkinv(fit)
      dimnames(fit) <- list(seq(nrow(newx)), s)
    }
    else {
      newb[ ,1] <- newb[ ,1] - newb[ ,-1,drop=FALSE] %*% object$mean
      fit <- newb
      dimnames(fit) <- list(s, object$xnames)
    }
    attr(fit, "s") <- s
    attr(fit, "fraction") <- sfrac
    attr(fit, "mode") <- mode
    fit
}

"cv.glmpath" <- function (x, y, family = binomial, weight = rep(1, length(y)), nfold = 10, fraction = seq(from=0, to=1, length=100), mode = c("norm","lambda"), plot.it = TRUE, se = TRUE, ...)
  {
    mode <- match.arg(mode)
    if (is.character(family)) 
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
        family <- family()
    else if (family$family=="gaussian") family <- gaussian()
    else if (family$family=="binomial") family <- binomial()
    else if (family$family=="poisson") family <- poisson()    
    if (family$family=="binomial") {
      uy <- unique(y)
      if (all(uy==c(1,-1)) | all(uy==c(-1,1))) {
        y[y==-1] <- 0
        cat("y=-1 converted to y=0 \n")
      }
    }    
    n <- length(y)
    all.folds <- split(sample(seq(n)), rep(1:nfold, length=n))
    residmat <- matrix(0, length(fraction), nfold)
    for (i in seq(nfold)) {
        omit <- all.folds[[i]]
        fit <- glmpath(x[-omit, ], y[-omit], family=family, weight=weight[-omit], ...)
        pred <- switch(mode, norm = {
          predict(fit, x[omit, , drop = FALSE], s = fraction,
                  type="response", mode = "norm.fraction")
        }, lambda = {
          predict(fit, x[omit, , drop = FALSE], s = fraction,
                  type="response", mode = "lambda.fraction")
        })
        if (length(omit) == 1) pred <- matrix(pred, nrow = 1)
        if (family$family=="binomial") ifelse(pred > 0.5, 1, 0)
        residmat[, i] <- apply((y[omit] - pred)^2, 2, mean)
        cat("CV Fold", i, "\n")
    }
    cv.error <- apply(residmat, 1, mean)
    cv.se <- sqrt(apply(residmat, 1, var)/nfold)
    object <- list(fraction = fraction, cv.error = cv.error, cv.se = cv.se, folds=all.folds)
    if (plot.it) {
      plot(fraction, cv.error, type="l", ylim=range(c(cv.error-cv.se, cv.error+cv.se)), xlab=switch(mode, norm="Norm fraction", lambda="log(lambda) fraction"), ylab="Cross-validation errors", main="Cross-validation errors")
      if (se) segments(fraction, cv.error-cv.se, fraction, cv.error+cv.se)
    }
    invisible(object)
}

"print.glmpath" <- function(x, ...)
  {
    actions <- x$actions
    xn <- x$xnames
    k <- length(actions)
    for (i in 1:k) {
      if (length(actions[[i]]) > 0) cat("Step",i,":",xn[abs(actions[[i]])],"\n")
    }
  }
