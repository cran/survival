
vcov.coxph<-function (object, ...) {
    rval<-object$var
    dimnames(rval)<-list(names(coef(object)),names(coef(object)))
    rval
}

vcov.survreg<-function (object, ...) {
        object$var
}

extractAIC.coxph.penal<- function(fit,scale,k=2,...){
    edf<-sum(fit$df)
    loglik <- fit$loglik[length(fit$loglik)]
    c(edf, -2 * loglik + k * edf)
}


 anova.coxph<-function (object, ...,  test = NULL) 
{
    if (length(object$rscore)>0)
        stop("Can't do anova tables with robust variances")
    dotargs <- list(...)
    named <- if (is.null(names(dotargs))) 
        rep(FALSE, length(dotargs))
    else (names(dotargs) != "")
    if (any(named)) 
        warning(paste("The following arguments to anova.coxph(..)", 
            "are invalid and dropped:", paste(deparse(dotargs[named]), 
                collapse = ", ")))
    dotargs <- dotargs[!named]
    is.coxmodel <- unlist(lapply(dotargs, function(x) inherits(x, 
        "coxph")))
    dotargs <- dotargs[is.coxmodel]
    if (length(dotargs) > 0) 
        return(anova.coxphlist(c(list(object), dotargs), test = test))
    varlist <- attr(object$terms, "variables")
    termlist<-attr(object$terms,"term.labels")
    resdev <- resdf <- NULL
    form<-".~."
    fenv<-environment(formula(object))
    nvars<-length(varlist)
    if (nvars > 1) {
        for (i in rev(termlist[-1])) {
            form<-paste(form,i,sep="-")
            fit <-update(object,as.formula(form,env=fenv))
            resdev <- c(resdev, -2*fit$loglik[2])
            resdf <- c(resdf,object$n-sum(!is.na(coef(fit))))
        }
    }
    resdf <- c(object$n, rev(resdf), object$n-sum(!is.na(coef(object))))
    resdev <- c(-2*object$loglik[1], rev(resdev), -2*object$loglik[2])
    table <- data.frame(c(NA, -diff(resdf)), c(NA, pmax(0, -diff(resdev))), 
        resdf, resdev)
    if (nvars == 0) 
        table <- table[1, , drop = FALSE]
    dimnames(table) <- list(c("NULL", attr(object$terms, "term.labels")), 
        c("Df", "Deviance", "Resid. Df", "Resid. Dev"))
    title <- paste("Analysis of Deviance Table\n Cox model: response is ",deparse(object$terms[[2]]),"\nTerms added sequentially (first to last)\n", 
        sep = "")
    df.dispersion <- Inf
    dispersion<-1
    if (!is.null(test)) 
        table <- stat.anova(table = table, test = test, scale = dispersion, 
            df.scale = df.dispersion, n =object$n)
    structure(table, heading = title, class = c("anova", "data.frame"))
}

anova.coxphlist<-function (object, ..., test = NULL) 
{
    responses <- as.character(lapply(object, function(x) {
        deparse(formula(x)[[2]])
    }))
    sameresp <- responses == responses[1]
    if (!all(sameresp)) {
        object <- object[sameresp]
        warning(paste("Models with response", deparse(responses[!sameresp]), 
            "removed because response differs from", "model 1"))
    }
    ns <- sapply(object, function(x) length(x$residuals))
    if (any(ns != ns[1])) 
        stop("models were not all fitted to the same size of dataset")
    nmodels <- length(object)
    if (nmodels == 1) 
        return(anova.glm(object[[1]], dispersion = dispersion, 
            test = test))
    resdf <- as.numeric(lapply(object, function(x) x$n-sum(!is.na(coef(x)))))
    resdev <- as.numeric(lapply(object, function(x) -2*x$loglik[2]))
    table <- data.frame(resdf, resdev, c(NA, -diff(resdf)), c(NA, 
        -diff(resdev)))
    variables <- lapply(object, function(x) paste(deparse(formula(x)), 
        collapse = "\n"))
    dimnames(table) <- list(1:nmodels, c("Resid. Df", "Resid. Dev", 
        "Df", "Deviance"))
    title <- "Analysis of Deviance Table\n"
    topnote <- paste("Model ", format(1:nmodels), ": ", variables, 
        sep = "", collapse = "\n")
    if (!is.null(test)) {
        bigmodel <- object[[order(resdf)[1]]]
        dispersion <-1
        df.dispersion <-  Inf
        table <- stat.anova(table = table, test = test, scale = dispersion, 
            df.scale = df.dispersion, n = length(bigmodel$residuals))
    }
    structure(table, heading = c(title, topnote), class = c("anova", 
        "data.frame"))
}
