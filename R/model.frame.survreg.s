# $Id: model.frame.survreg.S 11059 2008-10-23 12:32:50Z therneau $
model.frame.survreg <- function(formula, ...) {
    Call <- formula$call
    Call[[1]] <- as.name("model.frame")
    Call <- Call[c(1, match(c("formula", "data", "weights", "subset",
			   "na.action"), names(Call), 0))]
    if (is.R()) {
	dots <- list(...)
	nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0)]
	Call[names(nargs)] <- nargs
	env<-environment(formula$terms)
	if (is.null(env)) env<-parent.frame()
	eval(Call, env)
	}
    else eval(Call)
    }
