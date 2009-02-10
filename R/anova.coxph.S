# $Id$
# The anova function for a coxph object
anova.coxph <- function (object, ...,  test = NULL) {
    # All the ... args need to be coxph fits.  If any of them
    #  have a name attached, e.g., 'charlie=T' we assume a priori
    #  that they are illegal
    #
    dotargs <- list(...)
    named <- if (is.null(names(dotargs))) 
	           rep(FALSE, length(dotargs))
             else (names(dotargs) != "")
    if (any(named)) 
        warning(paste("The following arguments to anova.coxph(..)", 
            "are invalid and dropped:", paste(deparse(dotargs[named]), 
                collapse = ", ")))
    dotargs <- dotargs[!named]

    # Now check that they are all cox models
    is.coxmodel <- unlist(lapply(dotargs, function(x) inherits(x,  "coxph")))
    dotargs <- dotargs[is.coxmodel]
    # If multiple models remain, call a smarter function
    if (length(dotargs) > 0) 
        return(anova.coxphlist(c(list(object), dotargs), test = test))

    #
    # I have one Cox model 
    #
    if (length(object$rscore)>0)
        stop("Can't do anova tables with robust variances")

    varlist <- attr(object$terms, "variables")
    termlist<-attr(object$terms,"term.labels")
    resdev <- resdf <- NULL

    # Fit a series of Cox models, adding terms 1 by 1
    form<-".~."
    if (is.R()) fenv<-environment(formula(object))
    nvars<-length(varlist)
    if (nvars > 1) {
        for (i in rev(termlist[-1])) {
            form<-paste(form,i,sep="-") # add var to the formula
            if (is.R()) fit <-update(object, as.formula(form,env=fenv))
	    else        fit <- update(object, as.formula(form))
            resdev <- c(resdev, -2*fit$loglik[2])
            resdf <- c(resdf, object$n-sum(!is.na(coef(fit))))
        }
    }

    # Add the original model to the list & make a table
    resdf <- c(object$n, rev(resdf), object$n-sum(!is.na(coef(object))))
    resdev <- c(-2*object$loglik[1], rev(resdev), -2*object$loglik[2])
    table <- data.frame(c(NA, -diff(resdf)), c(NA, pmax(0, -diff(resdev))), 
			resdf, resdev)
    if (nvars == 0) 
        table <- table[1, , drop = FALSE]
    dimnames(table) <- list(c("NULL", attr(object$terms, "term.labels")), 
			    c("Df", "Deviance", "Resid. Df", "Resid. Dev"))
    title <- paste("Analysis of Deviance Table\n Cox model: response is ",
		   deparse(object$terms[[2]]),
		   "\nTerms added sequentially (first to last)\n", 
		   sep = "")
    df.dispersion <- Inf
    dispersion<-1
    if (!is.null(test)) 
        table <- stat.anova(table = table, test = test, scale = dispersion, 
            df.scale = df.dispersion, n =object$n)

    if (is.R())
	    structure(table, heading = title, class = c("anova", "data.frame"))
    else    structure(table, heading= title, class='anova')
}
