# SCCS @(#)labels.survreg.s	1.1 01/06/99
#labels.survreg <- labels.lm
labels.survreg<-function(x,...) attr(x$terms,"term.labels")
