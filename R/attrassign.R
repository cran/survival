
##
## redoes attr(modelmatrix,"assign") in the nice S-PLUS 3.4 format
##
attrassign<-function (object, ...) UseMethod("attrassign")

attrassign.lm<-function(lmobj){
	attrassign(model.matrix(lmobj),terms(lmobj))}

attrassign.default<-function(mmat,tt){
        if (!inherits(tt,"terms"))
                stop("need terms object")
        aa<-attr(mmat,"assign")
        if (is.null(aa))
                stop("argument is not really a model matrix")
        ll<-attr(tt,"term.labels")
        if (attr(tt,"intercept")>0)
                ll<-c("(Intercept)",ll)
        aaa<-factor(aa,labels=ll)
        split(order(aa),aaa)
	}

