#
# Good initial values are key to this data set
#   It killed v4 of survreg; 
#   data courtesy of Deborah Donnell, Fred Hutchinson Cancer Center
#

donnell <- scan("data.donnell", what=list(time1=0, time2=0, status=0))
donnell <- data.frame(donnell)

dfit <- survreg(Surv(time1, time2, status, type='interval') ~1, donnell)
summary(dfit)

