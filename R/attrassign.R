
##
## redoes attr(modelmatrix,"assign") in the nice S-PLUS 3.4 format
##
attrassign<-function (object, ...) UseMethod("attrassign")

attrassign.lm<-function(object,...){
	attrassign(model.matrix(object),terms(object))}

attrassign.default<-function(object,tt,...){
        if (!inherits(tt,"terms"))
                stop("need terms object")
        aa<-attr(object,"assign")
        if (is.null(aa))
                stop("argument is not really a model matrix")
        ll<-attr(tt,"term.labels")
        if (attr(tt,"intercept")>0)
                ll<-c("(Intercept)",ll)
        aaa<-factor(aa,labels=ll)
        split(order(aa),aaa)
	}

