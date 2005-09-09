#SCCS 03/25/97 @(#)summary.coxph.s	4.8
summarycoxph <-
 function(object, table = TRUE, coef = TRUE, conf.int = 0.95, scale = 1,
          digits=getOption("digits"),...)
    {
        cox<-object
    if (!is.null(cl<- cox$call)) {
	cat("Call:\n")
	dput(cl)
	cat("\n")
	}
    if (!is.null(cox$fail)) {
	cat(" Coxreg failed.", cox$fail, "\n")
	return()
	}
    savedig <- options(digits = digits)
    on.exit(options(savedig))

    omit <- cox$na.action
    if (length(omit))
	cat("  n=", cox$n, " (", naprint(omit), ")\n", sep="")
    else cat("  n=", cox$n, "\n")
    if (length(cox$icc))
	cat("  robust variance based on", cox$icc[1],
	    "groups, intra-class correlation =", format(cox$icc[2:3]), "\n")
    if (is.null(cox$coef)) {   # Null model
	cat ("   Null model\n")
	return()
	}

    beta <- cox$coef
    nabeta <- !(is.na(beta))          #non-missing coefs
    beta2 <- beta[nabeta]
    if(is.null(beta) | is.null(cox$var))
        stop("Input is not valid")
    se <- sqrt(diag(cox$var))
    if (!is.null(cox$naive.var)) nse <- sqrt(diag(cox$naive.var))

    if(coef) {
	if (is.null(cox$naive.var)) {
	    tmp <- cbind(beta, exp(beta), se, beta/se,
		   signif(1 - pchisq((beta/ se)^2, 1), digits -1))
	    dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)",
		"se(coef)", "z", "p"))
	    }
	else {
	    tmp <- cbind(beta, exp(beta), nse, se, beta/se,
		   signif(1 - pchisq((beta/ se)^2, 1), digits -1))
	    dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)",
		"se(coef)", "robust se", "z", "p"))
	    }
        cat("\n")
        prmatrix(tmp)
        }
    if(conf.int) {
        z <- qnorm((1 + conf.int)/2, 0, 1)
        beta <- beta * scale
        se <- se * scale
        tmp <- cbind(exp(beta), exp(-beta), exp(beta - z * se),
            exp(beta + z * se))
        dimnames(tmp) <- list(names(beta), c("exp(coef)", "exp(-coef)",
            paste("lower .", round(100 * conf.int, 2), sep = ""),
            paste("upper .", round(100 * conf.int, 2), sep = "")))
        cat("\n")
        prmatrix(tmp)
        }
    logtest <- -2 * (cox$loglik[1] - cox$loglik[2])
    sctest <- cox$score
    df <- length(beta2)
    cat("\n")
    cat("Rsquare=", format(round(1-exp(-logtest/cox$n),3)),
	"  (max possible=", format(round(1-exp(2*cox$loglik[1]/cox$n),3)),
	")\n" )
    cat("Likelihood ratio test= ", format(round(logtest, 2)), "  on ",
	df, " df,", "   p=", format(1 - pchisq(logtest, df)),
	"\n", sep = "")
    cat("Wald test            = ", format(round(cox$wald.test, 2)), "  on ",
	df, " df,", "   p=", format(1 - pchisq(cox$wald.test, df)),
	"\n", sep = "")
    cat("Score (logrank) test = ", format(round(sctest, 2)), "  on ", df,
        " df,", "   p=", format(1 - pchisq(sctest, df)), sep ="") 
    if (is.null(cox$rscore)) cat("\n\n")
    else cat(",   Robust = ", format(round(cox$rscore, 2)), 
	   "  p=", format(1 - pchisq(cox$rscore, df)), "\n\n", sep="")   

    if (!is.null(cox$naive.var))
	cat("  (Note: the likelihood ratio and score tests",
	  "assume independence of\n     observations within a cluster,",
	    "the Wald and robust score tests do not).\n")
    invisible()
    }

summary.coxph <-
 function(object,  coef = TRUE, conf.int = 0.95, scale = 1,
          digits = max(options()$digits - 4, 3),...)
{
     cox<-object
     beta <- cox$coef
     nabeta <- !(is.na(beta))          #non-missing coefs
     beta2 <- beta[nabeta]
     if(is.null(beta) | is.null(cox$var))
         stop("Input is not valid")
     se <- sqrt(diag(cox$var))
     if (!is.null(cox$naive.var)) nse <- sqrt(diag(cox$naive.var))

     rval<-list(call=cox$call,fail=cox$fail, na.action=cox$na.action,
                n=cox$n,icc=cox$icc)
     if(coef) {
         if (is.null(cox$naive.var)) {
	    tmp <- cbind(beta, exp(beta), se, beta/se,
		   signif(1 - pchisq((beta/ se)^2, 1), digits -1))
	    dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)",
		"se(coef)", "z", "p"))
	    }
	else {
	    tmp <- cbind(beta, exp(beta), nse, se, beta/se,
		   signif(1 - pchisq((beta/ se)^2, 1), digits -1))
	    dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)",
		"se(coef)", "robust se", "z", "p"))
	    }
        rval$coef<-tmp
        }
    if(conf.int) {
        z <- qnorm((1 + conf.int)/2, 0, 1)
        beta <- beta * scale
        se <- se * scale
        tmp <- cbind(exp(beta), exp(-beta), exp(beta - z * se),
            exp(beta + z * se))
        dimnames(tmp) <- list(names(beta), c("exp(coef)", "exp(-coef)",
            paste("lower .", round(100 * conf.int, 2), sep = ""),
            paste("upper .", round(100 * conf.int, 2), sep = "")))
        rval$conf.int<-tmp
    }
    df <- length(beta2)
     logtest <- -2 * (cox$loglik[1] - cox$loglik[2])
    rval$logtest <- c(test=logtest,
                      df=df,
                      pvalue=1 - pchisq(logtest, df))
    rval$sctest <- c(test=cox$score,
                     df=df,
                     pvalue=1 - pchisq(cox$score, df))
     rval$rsq<-c(rsq=1-exp(-logtest/cox$n),
                 maxrsq=1-exp(2*cox$loglik[1]/cox$n))
     rval$waldtest<-c(test=as.vector(round(cox$wald.test, 2)),
                      df=df,
                      pvalue=1 - pchisq(as.vector(cox$wald.test), df))
     if (!is.null(cox$rscore))
         rval$robscore<-c(test=cox$rscore,
                          df=df,
                          pvalue=1 - pchisq(cox$rscore, df))
     rval$used.robust<-!is.null(cox$naive.var)
     class(rval)<-"summary.coxph"
     rval
 }


print.summary.coxph <-
 function(x, digits = max(options()$digits - 4, 3), ...)
    {
    if (!is.null(x$call)) {
	cat("Call:\n")
	dput(x$call)
	cat("\n")
    }
    if (!is.null(x$fail)) {
	cat(" Coxreg failed.", x$fail, "\n")
	return()
    }
    savedig <- options(digits = digits)
    on.exit(options(savedig))

    omit <- x$na.action
    if (length(omit))
	cat("  n=", x$n, " (", naprint(omit), ")\n", sep="")
    else cat("  n=", x$n, "\n")
    if (length(x$icc))
	cat("  robust variance based on", x$icc[1],
	    "groups\n")
    if (nrow(x$coef)==0) {   # Null model
	cat ("   Null model\n")
	return()
    }


    if(!is.null(x$coef)) {
        prmatrix(x$coef)
        }
    if(!is.null(x$conf.int)) {
        cat("\n")
        prmatrix(x$conf.int)
        }
    cat("\n")
    cat("Rsquare=", format(round(x$rsq["rsq"],3)),
	"  (max possible=", format(round(x$rsq["maxrsq"],3)),
        ")\n" )
    cat("Likelihood ratio test= ", format(round(x$logtest["test"], 2)), "  on ",
	x$logtest["df"], " df,", "   p=", format(x$logtest["pvalue"]),
        "\n", sep = "")
    cat("Wald test            = ", format(round(x$waldtest["test"], 2)), "  on ",
	x$waldtest["df"], " df,", "   p=", format(x$waldtest["pvalue"]),
	"\n", sep = "")
    cat("Score (logrank) test = ", format(round(x$sctest["test"], 2)), "  on ",
        x$sctest["df"]," df,", "   p=", format(x$sctest["pvalue"]), sep ="")
    if (is.null(x$robscore))
        cat("\n\n")
    else cat(",   Robust = ", format(round(x$robscore["test"], 2)), 
             "  p=", format(x$robscore["pvalue"]), "\n\n", sep="")   

    if (x$used.robust)
	cat("  (Note: the likelihood ratio and score tests",
            "assume independence of\n     observations within a cluster,",
	    "the Wald and robust score tests do not).\n")
    invisible()
}
