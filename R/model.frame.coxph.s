#  $Id: model.frame.coxph.S 11059 2008-10-23 12:32:50Z therneau $
# This is here because the white book (Models in S) says we should have
#  one, but I am not sure that anything ever calls this routine.
model.frame.coxph <- function(formula, ...) {
    Call <- formula$call
    Call[[1]] <- as.name("model.frame")
    Call <- Call[c(1, match(c("formula", "data", "weights", "subset",
			   "na.action"), names(Call), 0))]

    if (is.R()) {
	# Thomas added these lines, greek to me (TMT)
	dots <- list(...)
	nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0)]
	Call[names(nargs)] <- nargs
	# coxph has a 'formula' component so this is OK
	env <- environment(formula(formula))
	if (is.null(env)) env <- parent.frame()
	eval(Call, env)
	}
    else eval(Call)
    }
