#  SCCS  @(#)model.frame.coxph.s	4.4 02/21/99
model.frame.coxph <- function(formula, ...) {
    Call <- formula$call
    Call[[1]] <- as.name("model.frame")
    Call <- Call[match(c("", "formula", "data", "weights", "subset",
			   "na.action"), names(Call), 0)]
    env<-environment(formula(formula))
    if (is.null(env)) env<-parent.frame()
    eval(Call, env)

    }
