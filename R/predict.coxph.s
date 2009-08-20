# Automatically generated from all.nw using noweb
predict.coxph <- function(object, newdata, 
                       type=c("lp", "risk", "expected", "terms"),
                       se.fit=FALSE, na.action=na.pass,
                       terms=names(object$assign), collapse, ...) {
    if (!inherits(object, 'coxph'))
        stop("Primary argument much be a coxph object")

    type <-match.arg(type)
    n <- object$n
    Terms <-  object$terms

    if (!missing(terms)) {
        if (is.numeric(terms)) {
            if (any(terms != floor(terms) | 
                    terms > length(object$assign) |
                    terms <1)) stop("Invalid terms argument")
            }
        else if (any(is.na(match(terms, names(object$assign)))))
           stop("a name given in the terms argument not found in the model")
        }

    #Do I have strata or a cluster statment?  Don't conflict names with the strata()
    #  function
    strat <- attr(Terms, 'specials')$strata
    dropx <- NULL
    if (length(strat)) {
        stemp <- untangle.specials(Terms, 'strata', 1)
        dropx <- stemp$terms
        }
    if (length(attr(Terms, 'specials')$cluster)) {
        temp <- untangle.specials(Terms, 'cluster', 1)
        dropx <- c(dropx, temp$terms)
        }
    if (length(dropx)) Terms2 <- Terms[-dropx]
    else  Terms2 <- Terms
    if (type != 'expected') Terms2 <- delete.response(Terms2)

    na.action.used <- object$na.action
    has.offset <- !is.null(attr(Terms, 'offset'))
    if (type == 'expected') {
        y <- object[['y']]
        if (is.null(y)) {  # very rare case
            mf <- model.frame(object)
            y <-  model.extract(mf, 'response')
            }
        }

    if (se.fit || !missing(newdata) || type=='terms' || length(strat)) {
        need.x <- TRUE
        if (is.null(object[['x']])) {
            mf <- model.frame(object)
            x <- model.matrix(delete.response(Terms2), mf,
                          contr=object$contrasts)[,-1,drop=FALSE]
            if (length(strat)) {
                if (length(stemp$vars)==1) oldstrat <- mf[[stemp$vars]]
                else oldstrat <- strata(mf[,stemp$vars], shortlabel=TRUE)
                }
            else oldstrat <- rep(0L, nrow(mf))
            weights <- model.weights(mf)
            offset <- model.offset(mf)
            }
        else {
            x <- object[['x']]
            oldstrat <- object$strata
            weights <- object$weights
            offset <-  object$offset
            }
        if (is.null(weights)) weights <- rep(1.0, nrow(x))
        if (is.null(offset))  offset  <- rep(0.0, nrow(x))
        if (length(oldstrat)==0) oldstrat <- rep(0L, nrow(x))
        }
    else {
        oldstrat <- rep(0L, n)
        offset <- 0
        need.x <- FALSE
        }

    if (!missing(newdata)) {
        mf <- model.frame(Terms2, data=newdata, xlev=object$xlevels,
                          na.action=na.action)
        newx <- model.matrix(Terms2, mf,
                             contr=object$contrasts)[,-1,drop=FALSE]
        if (length(strat)) {
            if (length(stemp$vars)==1) newstrat <- mf[[stemp$vars]]
            else newstrat <- strata(mf[,stemp$vars], shortlabel=TRUE)
            if (any(is.na(match(newstrat, oldstrat)))) 
                stop("New data has a strata not found in the original model")
            }
        else newstrat <- rep(0L, nrow(mf))

        newoffset <- model.offset(mf) 
        if (is.null(newoffset)) newoffset <- rep(0.0, nrow(newx))
        na.action.used <- attr(mf, 'na.action')
        if (type== 'expected') {
            newy <- model.response(mf)
            if (attr(newy, 'type') != attr(y, 'type'))
                stop("New data has a different survival type than the model")
            }
        }  
    if (type=='expected') {
        if (missing(newdata))
            pred <- y[,ncol(y)] - object$residuals
        if (!missing(newdata) || se.fit) {
            ustrata <- unique(oldstrat)
            risk <- exp(object$linear.predictors)
            x <- scale(x, center=object$means, scale=FALSE)
            if (se.fit) { 
                se <- double(nrow(mf))
                }
            if (missing(newdata))
                se <- double(nrow(mf))
            if (!missing(newdata)) {
                se <- double(nrow(mf))
                pred <- se
                newx <- scale(newx, center=object$means, scale=FALSE)
                newrisk <- c(exp(newx %*% object$coef))
                }
            survtype= ifelse(fit$method=='efron', 3,2)
            for (i in ustrata) {
                indx <- which(oldstrat == i)
                afit <- agsurv(y[indx,,drop=F], x[indx,,drop=F], 
                                              weights[indx], risk[indx],
                                              survtype, survtype)
                afit.n <- length(afit$time)
                if (missing(newdata)) { 
                    # In this case we need se.fit, nothing else
                    j1 <- approx(afit$time, 1:afit.n, y[indx,1], method='constant',
                                 f=0, yleft=0, yright=afit.n)$y
                    chaz <- c(0, afit$cumhaz)[j1 +1]
                    varh <- c(0, cumsum(afit$varhaz))[j1 +1]
                    xbar <- rbind(0, afit$xbar)[j1+1,,drop=F]
                    if (ncol(y)==2) {
                        dt <- (chaz * x[indx,]) - xbar
                        se[indx] <- sqrt(varh + rowSums((dt %*% object$var) *dt)) *
                            risk[indx]
                        }
                    else {
                        j2 <- approx(afit$time, 1:afit.n, y[indx,2], method='constant',
                                 f=0, yleft=0, yright=afit.n)$y
                        chaz2 <- c(0, afit$cumhaz)[j2 +1]
                        varh2 <- c(0, cumsum(afit$varhaz))[j2 +1]
                        xbar2 <- rbind(0, afit$xbar)[j2+1,,drop=F]
                        dt <- (chaz * x[indx,]) - xbar
                        v1 <- varh +  rowSums((dt %*% object$var) *dt)
                        dt2 <- (chaz2 * x[indx,]) - xbar2
                        v2 <- varh2 + rowSums((dt2 %*% object$var) *dt2)
                        se[indx] <- sqrt(v2-v1)* risk[indx]
                        }
                    }

                else {
                    #there is new data
                    indx2 <- which(newstrat == i)
                    j1 <- approx(afit$time, 1:afit.n, newy[indx2,1], 
                                 method='constant', f=0, yleft=0, yright=afit.n)$y
                    chaz <-c(0, afit$cumhaz)[j1+1]
                    pred[indx2] <- chaz * newrisk[indx2]
                    if (se.fit) {
                        varh <- c(0, cumsum(afit$varhaz))[j1+1]
                        xbar <- rbind(0, afit$xbar)[j1+1,,drop=F]
                        }
                    if (ncol(y)==2) {
                        if (se.fit) {
                            dt <- (chaz * newx[indx2,]) - xbar[indx2,]
                            se[indx2] <- sqrt(varh + rowSums((dt %*% object$var) *dt)) *
                                newrisk[indx2]
                            }
                        }
                    else {
                        j2 <- approx(afit$time, 1:afit.n, newy[indx2,2], 
                                 method='constant', f=0, yleft=0, yright=afit.n)$y
                                    chaz2 <- approx(-afit$time, afit$cumhaz, -newy[indx2,2],
                                   method="constant", rule=2, f=0)$y
                        chaz2 <-c(0, afit$cumhaz)[j2+1]
                        pred[indx2] <- (chaz2 - chaz) * newrisk[indx2]
                    
                        if (se.fit) {
                            varh2 <- c(0, cumsum(afit$varhaz))[j1+1]
                            xbar2 <- rbind(0, afit$xbar)[j1+1,,drop=F]
                            dt <- (chaz * newx[indx2,]) - xbar[indx2,]
                            dt2 <- (chaz2 * newx[indx2,]) - xbar2[indx2,]

                            v2 <- varh2 + rowSums((dt2 %*% object$var) *dt2)
                            v1 <- varh +  rowSums((dt %*% object$var) *dt)
                            se[indx2] <- sqrt(v2-v1)* risk[indx2]
                            }
                        }
                    }
                }
            }
        }
    if (is.null(object$coefficients))
        coef<-numeric(0)
    else {
        # Replace any NA coefs with 0, to stop NA in the linear predictor
        coef <- ifelse(is.na(object$coefficients), 0, object$coefficients)
        }

    if (missing(newdata)) {
        offset <- offset - mean(offset)
        if (length(strat)) {
            for (i in unique(oldstrat)) {
                j <- which(oldstrat==i)
                if (length(j)==1) x[j,] <- 0   # only 1 subject in the strata
                else if (length(j) >1) {
                    xmean <- colSums(x[j,] * weights[j]) / sum(weights[j])
                    x[j,]  <- scale(x[j,], center=xmean, scale=FALSE)
                    }
                }
            newx <- x
            }
        else if (need.x) newx <- scale(x, center=object$means, scale=FALSE)
        }
    else {
        offset <- newoffset - mean(offset)
        if (length(strat)) {
            for (i in unique(oldstrat)) {
                j <- which(newstrat==i)    #no matches is a possibility
                if (length(j) ==1) newx[j,] <- 0
                else if (length(j) >1) { 
                    xmean <- colSums(x[j,,drop=FALSE] * weights[j]) / sum(weights[j])
                    newx[j,]  <- scale(newx[j,], center=xmean, scale=FALSE)
                    }
                }
            }
        else newx <- scale(newx, center=object$means, scale=FALSE)
        }

    if (type=='lp' || type=='risk') {
        if (need.x) pred <- newx %*% coef + offset
        else pred <- object$linear.predictors
        if (se.fit) se <- sqrt(rowSums((newx %*% object$var) *newx))

        if (type=='risk') {
            pred <- exp(pred)
            if (se.fit) se <- se * sqrt(pred)  # standard Taylor series approx
            }
        }
    else if (type=='terms') { 
        asgn <- object$assign
        nterms<-length(asgn)
        pred<-matrix(ncol=nterms,nrow=NROW(newx))
        dimnames(pred) <- list(rownames(newx), names(asgn))
        if (se.fit) se <- pred
        
        for (i in 1:nterms) {
            tt <- asgn[[i]]
            tt <- tt[!is.na(object$coefficients[tt])]
            xtt <- newx[,tt, drop=F]
            pred[,i] <- xtt %*% object$coefficient[tt]
            if (se.fit)
                se[,i] <- sqrt(rowSums((xtt %*% object$var[tt,tt]) *xtt))
            }
        pred <- pred[,terms, drop=F]
        if (se.fit) se <- se[,terms, drop=F]
        
        attr(pred, 'constant') <- sum(object$coefficients*object$means, na.rm=T)
        }
    if (type != 'terms') {
        pred <- drop(pred)
        if (se.fit) se <- drop(se)
        }

    if (!is.null(na.action.used)) {
        pred <- naresid(na.action.used, pred)
        if (is.matrix(pred)) n <- nrow(pred)
        else               n <- length(pred)
        if(se.fit) se <- naresid(na.action.used, se)
        }

    if (!missing(collapse)) {
        if (length(collapse) != n) stop("Collapse vector is the wrong length")
        pred <- rowsum(pred, collapse)  # in R, rowsum is a matrix, always
        if (se.fit) se <- sqrt(rowsum(se^2, collapse))
        if (type != 'terms') {
            pred <- drop(pred)
            if (se.fit) se <- drop(se)
            }
        }

    if (se.fit) list(fit=pred, se.fit=se)
    else pred

    }
