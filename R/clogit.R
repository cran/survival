## conditional logistic regression
##
## case~exposure+strata(matching)
##

clogit<-function(formula,data,method=c("exact","approximate"),
                 na.action=getOption("na.action"),subset=NULL){
    mf<-match.call()
    mf[[1]]<-as.name("model.frame")
    mf$method<-NULL
    mf<-eval(mf,sys.frame(sys.parent()))
    
    coxcall<-match.call()
    coxcall[[1]]<-as.name("coxph")
    newformula<-formula
    newformula[[2]]<-substitute(Surv(rep(1,nn),case),list(case=formula[[2]],nn=NROW(mf)))
    environment(newformula)<-environment(formula)
    coxcall$formula<-newformula
    coxcall$method<-switch(match.arg(method),exact="exact","breslow")

    coxcall<-eval(coxcall,sys.frame(sys.parent()))
    coxcall$call<-sys.call()
    
    class(coxcall)<-c("clogit","coxph")
    coxcall
}
