#  SCCS  @(#)model.frame.coxph.s	4.4 02/21/99
model.frame.coxph <- function(formula, ...) {
    Call <- formula$call
    Call[[1]] <- as.name("model.frame")
    Call <- Call[match(c("", "formula", "data", "weights", "subset",
			   "na.action"), names(Call), 0)]
    dots <- list(...)
    nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0)]
    Call[names(nargs)] <- nargs
    # coxph has a 'formula' component so this is OK
    env <- environment(formula(formula))
    if (is.null(env)) env <- parent.frame()
    eval(Call, env)
    }
