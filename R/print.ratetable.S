## $Id: print.ratetable.S 10788 2008-09-18 00:48:23Z therneau $

print.ratetable <- function(x, ...) {
  cat ("Rate table with dimension(s):", attr(x, 'dimid'), "\n")
  attributes(x) <- attributes(x)[c("dim", "dimnames")]
  NextMethod()
}
