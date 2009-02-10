#
# This function is supplied for SAS converts who don't read
#  the documentation, and can't find the baseline hazard
#
basehaz <- function(fit, centered=TRUE) {
    if(!inherits(fit,"coxph"))
        stop("must be a coxph object")

    sfit<-survfit(fit)

    H<- -log(sfit$surv)

    strata<-sfit$strata
    if (!is.null(strata))
        strata<-rep(names(strata),strata)
    
    if (!centered){
        z0<-fit$means
        bz0<-sum(z0*coef(fit))
        H<- H*exp(-bz0)
    }

    if (is.null(strata))
      return(data.frame(hazard=H,time=sfit$time))
    else
      return(data.frame(hazard=H,time=sfit$time,strata=strata))

    }
