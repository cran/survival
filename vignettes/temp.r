#
# test of multiple tt calls
#
tdata <- lung
tdata$age2 <- lung$age

fit1 <- coxph(Surv(time, status) ~ sex + tt(age) + tt(age2), tdata,
             tt=list(function(x, t, ...) x+ t/365, 
                     function(x, t, ...) (x+t)^2))

fit2 <- coxph(Surv(time, status) ~ sex + tt(age) + tt(age), tdata,
             tt=list(function(x, t, ...) x+ t/365, 
                     function(x, t, ...) (x+t)^2))
