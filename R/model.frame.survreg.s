# SCCS @(#)model.frame.survreg.s	1.1 11/25/98
model.frame.survreg <- function(object, ...) {
    Call <- object$call
    Call[[1]] <- as.name("model.frame")
    Call <- Call[match(c("", "formula", "data", "weights", "subset",
			   "na.action"), names(Call), 0)]
    env<-environment(Call$formula)
    if (is.null(env)) env<-parent.frame()
    eval(Call, env)
    }
