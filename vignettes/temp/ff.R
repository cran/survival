# 
# this is a wild problem pointed out by a user
#

dummy <- function(object, ...) {
    UseMethod('dummy', object)
}

dummy.coxph <-function(object, newdata, times) {
    if (is.null(object[['y']])) stop("outcome data not in coxph object")
    Y  <- object[['y']]
    ny <- ncol(Y)
    
    if (missing(times)) times <- max(Y[, ny-1])
    if (missing(newdata)) {
        newdata <- data.frame(Y= Y)
        kmfit <- survfit(Y ~ 1, newdata)
    }       else {stop("not done")}
    
    ps <- pseudo(kmfit, times=times)
    1- ps[,1]
}

test <- coxph(Surv(time, status) ~ x, aml)
