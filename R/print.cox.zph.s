# SCCS @(#)print.cox.zph.s	4.5 09/27/96
print.cox.zph <- function(x, digits = max(options()$digits - 4, 3),...)
    invisible(print(x$table, digits=digits))
