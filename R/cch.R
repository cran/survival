### Suite of programs for case-cohort analysis

### Main program

cch <- function(formula, data=sys.parent(), subcoh, id, cohort.size, 
                method=c("Prentice", "SelfPrentice", "LinYing")){
    call <- match.call()
    
    if (is.data.frame(data)){
        if (inherits(id,"formula"))
            id<-model.frame(id,data)[,1]
        if (inherits(subcoh,"formula"))
            subcoh<-model.frame(subcoh,data)[,1]
    }
    ## Check id, subcoh and cohort.size variables
    if(length(id)!=length(unique(id)))
        stop("Multiple records per id not allowed")
    if (is.logical(subcoh))
        subcoh <- as.numeric(subcoh)
    tt <- table(subcoh)
    if(min(charmatch(names(tt), c("0","1"), 0))==0)
        stop("Permissible values for subcohort indicator are 0/1 or TRUE/FALSE")
    if(length(id)>cohort.size)
        stop("Number of records greater than cohort size")
    nn <- cohort.size	
    ## Evaluate model formula
    m <- match.call(expand.dots=FALSE)
    m$method <- m$cohort.size <- m$id <- m$subcoh <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval(m,sys.parent())
    Terms <- attr(m,"terms")
    Y <- model.extract(m, response)
    if(!inherits(Y, "Surv"))
        stop("Response must be a survival object")
    type <- attr(Y, "type")
    itype<-charmatch(type,c("right","counting"),nomatch=0)
    cens<-switch(itype+1,
                 stop(paste("Cox model doesn't support \"", type, "\" survival data", sep = "")),
                 Y[,2],
                 Y[,3])
    cc<-cens+1-subcoh
    texit<-switch(itype+1, stop(), Y[,1], Y[,2])
    tenter<-switch(itype+1, stop(), rep(0,length(texit)), Y[,1])
    X <- model.matrix(Terms, m)
    X <- X[,2:ncol(X)]
    method<-match.arg(method)
    fitter <- get(method)
    out<-fitter(tenter=tenter, texit=texit, cc=cc, id=id, X=X, ntot=nn)
    ##out <- eval(z, sys.parent())
    out$method <- method
    names(out$coef) <- dimnames(X)[[2]]
    if(!is.null(out$var))
        dimnames(out$var) <- list(dimnames(X)[[2]], dimnames(X)[[2]])
    if(!is.null(out$naive.var))
        dimnames(out$naive.var) <- list(dimnames(X)[[2]], dimnames(X)[[2]])
    out$call <- call
    out$cohort.size <- nn
    out$subcohort.size <- tt[2]
    class(out) <- "cch"
    out
}

### Subprograms

Prentice <- function(tenter, texit, cc,  id, X, ntot){
    eps <- 0.00000001
    cens <- as.numeric(cc>0) # Censorship indicators
    subcoh <- as.numeric(cc<2) # Subcohort indicators

    ## Calculate Prentice estimate
    ent2 <- tenter
    ent2[cc==2] <- texit[cc==2]-eps
    fit1 <- coxph(Surv(ent2,texit,cens)~X,eps=eps,x=TRUE)

    ## Calculate Prentice estimate and variance
    nd <- sum(cens) # Number of failures
    nc <- sum(subcoh) # Number in subcohort
    ncd <- sum(cc==1) #Number of failures in subcohort
    X <- as.matrix(X)
    aent <- c(tenter[cc>0],tenter[cc<2])
    aexit <- c(texit[cc>0],texit[cc<2])
    aX <- rbind(as.matrix(X[cc>0,]),as.matrix(X[cc<2,]))
    aid <- c(id[cc>0],id[cc<2])
    dum <- rep(-100,nd)
    dum <- c(dum,rep(0,nc))
    gp <- rep(1,nd)
    gp <- c(gp,rep(0,nc))
    fit <- coxph(Surv(aent,aexit,gp)~aX+offset(dum)+cluster(aid),eps=eps,x=TRUE,
                 iter.max=35,init=fit1$coefficients)
    db <- resid(fit,type="dfbeta")
    db <- as.matrix(db)
    db <- db[gp==0,]
    fit$naive.var <- fit$naive.var+(1-(nc/ntot))*t(db)%*%(db)
    fit$coef <- fit1$coefficients
    fit
}

SelfPrentice <- function(tenter, texit, cc,  id, X, ntot){
    eps <- 0.00000001
    cens <- as.numeric(cc>0) # Censorship indicators
    subcoh <- as.numeric(cc<2) # Subcohort indicators

    ## Calculate Self-Prentice estimate and variance
    nd <- sum(cens) # Number of failures
    nc <- sum(subcoh) # Number in subcohort
    ncd <- sum(cc==1) #Number of failures in subcohort
    X <- as.matrix(X)
    aent <- c(tenter[cc>0],tenter[cc<2])
    aexit <- c(texit[cc>0],texit[cc<2])
    aX <- rbind(as.matrix(X[cc>0,]),as.matrix(X[cc<2,]))
    aid <- c(id[cc>0],id[cc<2])
    dum <- rep(-100,nd)
    dum <- c(dum,rep(0,nc))
    gp <- rep(1,nd)
    gp <- c(gp,rep(0,nc))
    fit <- coxph(Surv(aent,aexit,gp)~aX+offset(dum)+cluster(aid),eps=eps,x=TRUE)
    db <- resid(fit,type="dfbeta")
    db <- as.matrix(db)
    db <- db[gp==0,,drop=FALSE]
    fit$naive.var <- fit$naive.var+(1-(nc/ntot))*t(db)%*%(db)
    fit
}

LinYing <- function(tenter, texit, cc,  id, X, ntot){
    eps <- 0.00000001
    cens <- as.numeric(cc>0) # Censorship indicators
    subcoh <- as.numeric(cc<2) # Subcohort indicators
    nd <- sum(cens) # Number of failures
    nc <- sum(subcoh) # Number in subcohort
    ncd <- sum(cc==1) #Number of failures in subcohort

    ## Calculate Lin-Ying estimate and variance
    offs <- rep((ntot-nd)/(nc-ncd),length(texit))
    offs[cc>0] <- 1
    loffs <- log(offs)
    fit <- coxph(Surv(tenter, texit, cens)~X+offset(loffs)+cluster(id),
                 eps=eps,x=TRUE)
    db <- resid(fit,type="dfbeta")
    db <- as.matrix(db)
    db <- db[cens==0,,drop=FALSE]
    dbm <- apply(db,2,mean)
    db <- sweep(db,2,dbm)
    fit$naive.var <- fit$naive.var+(1-(nc-ncd)/(ntot-nd))*t(db)%*%db
    fit
}

vcov.cch<-function(object,...) object$var

"print.cch"<- function(x,...)
{
    ## produces summary from an x of the class "cch"
    call<-x$call
    coef <- x$coef
    method <- x$method
    se <- sqrt(diag(x$var))
    Z<- coef/se
    p<- pnorm(Z)
    cohort.size<-x$cohort.size
    subcohort.size<-x$subcohort.size
    coefficients <- matrix(0, nrow = length(coef), ncol = 4)
    dimnames(coefficients) <- list(names(coef), c("Value", 
                                                  "SE", "Z", "p"))
    coefficients[, 1] <- coef
    coefficients[, 2] <- se
    coefficients[, 3] <- Z
    coefficients[, 4] <- 2*(1-p)

    cat("Case-cohort analysis, ",x$method," method,\n with subcohort of",
        subcohort.size,"from cohort of", cohort.size,"\n\n")
    cat("Call: "); print(x$call)
    cat("\nCoefficients:\n")
    print(coefficients)
    invisible(x)
}


"summary.cch"<-function(object,...)
{
    ## produces summary from an object of the class "cch"
    call<-object$call
	coef <- object$coef
	method <- object$method
	se <- sqrt(diag(object$var))
      Z<- coef/se
      p<- pnorm(Z)
      cohort.size<-object$cohort.size
      subcohort.size<-object$subcohort.size
    coefficients <- matrix(0, nrow = length(coef), ncol = 4)
    dimnames(coefficients) <- list(names(coef), c("Value", 
                                                  "SE", "Z", "p"))
    coefficients[, 1] <- coef
    coefficients[, 2] <- se
    coefficients[, 3] <- Z
    coefficients[, 4] <- 2*(1-p)
    structure(list(call=call, method=method, cohort.size=cohort.size,
                   subcohort.size=subcohort.size, coefficients = coefficients), 
              class = "summary.cch")
}

print.summary.cch <- function(x,digits=3,...){

    cat("Case-cohort analysis, ",x$method," method,\n with subcohort of",
        x$subcohort.size,"from cohort of", x$cohort.size,"\n\n")
    cat("Call: "); print(x$call)
    cat("\nCoefficients:\n")
    output<-cbind(Coef=x$coefficients[,1],HR=exp(x$coefficients[,1]),
                  "(95%"=exp(x$coefficients[,1]-2*x$coefficients[,2]),
                  "CI)"=exp(x$coefficients[,1]+2*x$coefficients[,2]),
                  "p"=x$coefficients[,4]
                  )
    print(round(output,3))
    invisible(x)
}

