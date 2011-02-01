# Automatically generated from all.nw using noweb
survConcordance <- function(formula, data,
                            weights, subset, na.action) {
    Call <- match.call()  # save a copy of of the call, as documentation

    m <- match.call(expand=FALSE)
    m[[1]] <- as.name("model.frame")
    m$formula <- if(missing(data)) terms(formula, "strata")
                 else              terms(formula, "strata", data=data)
    m <- eval(m, sys.parent())
    Terms <- attr(m, 'terms')

    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) stop("Response must be a survival object")
    n <- nrow(Y)

    wt <- model.extract(m, 'weights')
    if (length(wt) ==0) wt <- rep(1., n)
    offset<- attr(Terms, "offset")
    if (length(offset)>0) stop("Offset terms not allowed")

    stemp <- untangle.specials(Terms, 'strata')
    if (length(stemp$vars)) {
        if (length(stemp$vars)==1) strat <- m[[stemp$vars]]
        else strat <- strata(m[,stemp$vars], shortlabel=TRUE)
        Terms <- Terms[-stemp$terms]
    }
    else strat <- NULL
    
    x <- model.matrix(Terms, m)[,-1, drop=FALSE]  #remove the intercept
    if (ncol(x) > 1) stop("Only one predictor variable allowed")

    btree <- function(n) {
       tfun <- function(n, id, power) {
           if (n==1) id
           else if (n==2) c(2L *id, id)
           else if (n==3) c(2L*id, id, 2L*id +1L)
           else {
               nleft <- if (n== power*2) power  else min(power-1, n-power/2)
               c(tfun(nleft, 2L *id, power/2), id,
                 tfun(n-(nleft+1), 2L*id +1L, power/2))
               }
           }
       tfun(n, 1L, 2^(floor(logb(n-1,2))))
       }

    docount <- function(stime, risk, wts) {
        if (attr(stime, 'type') == 'right') {
            ord <- order(stime[,1], -stime[,2])
            ux <- sort(unique(risk))
            n2 <- length(ux)
            index <- btree(n2)[match(risk[ord], ux)] - 1L
             .Call('concordance1', stime[ord,], wts[ord], index, length(ux))
        }
        else if (attr(stime, 'type') == "counting") {
            sort.stop <- order(-stime[,2], stime[,3])
            sort.start <- order(-stime[,1])
            ux <- sort(unique(risk))
            n2 <- length(ux)
            index <- btree(n2)[match(risk, ux)] - 1L
            .Call('concordance2', stime, wts, index, length(ux),
                           sort.stop-1L, sort.start-1L)
        }
        else stop("Invalid survival type for concordance")
    }
        
    if (is.null(strat)) {
        count <- docount(Y, x, wt)
        names(count) <- c("concordant", "discordant", "tied risk", "tied time")
        concordance <- (count[1] + count[3]/2)/sum(count[1:3])
    }
    else {
        ustrat <- levels(strat)[table(strat) >0]  #some strata may have 0 obs
        count <- matrix(0., nrow=length(ustrat), ncol=4)
        for (i in 1:length(ustrat)) {
            keep <- which(strat == ustrat[i])
            count[i,] <- docount(Y[keep,,drop=F], x[keep], wt[keep])
        }
        dimnames(count) <- list(ustrat,  c("concordant", "discordant",
                                           "tied risk", "tied time"))
        temp <- colSums(count)
        concordance <- (temp[1] + temp[3]/2)/ sum(temp[1:3])
    }
    
    fit <- list(concordance= concordance, stats=count, n=n, call=Call)
    na.action <- attr(m, "na.action")
    if (length(na.action)) fit$na.action <- na.action

    oldClass(fit) <- 'survConcordance'
    fit
}

print.survConcordance <- function(x, ...) {
    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
        cat("\n")
        }
    omit <- x$na.action
    if(length(omit))
        cat("  n=", x$n, " (", naprint(omit), ")\n", sep = "")
    else cat("  n=", x$n, "\n")
    cat("Concordance= ", format(x$concordance), '\n', sep='')
    print(x$stats)

    invisible(x)
    }
