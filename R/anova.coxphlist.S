#
# This is usually called from anova.coxph, not a user
#  It's first argument must be a list of coxph models
anova.coxphlist <- function (object, test =  NULL,...) {
    if (!is.list(object)) stop("First argument must be a list")
    if (!all(lapply(object, function(x) inherits(x, 'coxph'))))
	     stop("Argument must be a list of coxph models")
    
    responses <- as.character(lapply(object, 
				     function(x) deparse(formula(x)[[2]])))
    sameresp <- (responses == responses[1])
    if (!all(sameresp)) {
        object <- object[sameresp]
        warning(paste("Models with response", deparse(responses[!sameresp]), 
            "removed because response differs from", "model 1"))
    }
    ns <- sapply(object, function(x) length(x$residuals))
    if (any(ns != ns[1])) 
        stop("models were not all fitted to the same size of dataset")
    nmodels <- length(object)
    if (nmodels == 1) # only one model remains
        return(anova.coxph(object[[1]], test = test))
    resdf <- as.numeric(lapply(object, function(x) x$n-sum(!is.na(coef(x)))))
    resdev <- as.numeric(lapply(object, function(x) -2*x$loglik[2]))
    table <- data.frame(resdf, resdev, c(NA, -diff(resdf)), 
			c(NA, -diff(resdev)))
    variables <- lapply(object, 
			function(x) paste(deparse(formula(x)),collapse = "\n"))
    dimnames(table) <- list(1:nmodels, 
			    c("Resid. Df", "Resid. Dev", "Df", "Deviance"))
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
    if (is.R()) structure(table, heading = c(title, topnote), 
			  class = c("anova", "data.frame"))
    else structure(table, heading = c(title, topnote), 
			  class = "anova")	     
} 

