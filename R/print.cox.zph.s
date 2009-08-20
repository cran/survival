# $Id: print.cox.zph.S 11059 2008-10-23 12:32:50Z therneau $
print.cox.zph <- function(x, digits = max(options()$digits - 4, 3),...)
    invisible(print(x$table, digits=digits))
