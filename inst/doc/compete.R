### R code from vignette source 'compete.Rnw'

###################################################
### code chunk number 1: compete.Rnw:23-29
###################################################
options(continue="  ", width=60)
options(SweaveHooks=list(fig=function() par(mar=c(4.1, 4.1, .3, 1.1))))
pdf.options(pointsize=10) #text in graph about the same as regular text
options(contrasts=c("contr.treatment", "contr.poly")) #ensure default

require("survival")


###################################################
### code chunk number 2: sfig1
###################################################
getOption("SweaveHooks")[["fig"]]()
# A note to readers of this code: drawing multi-state figures in this
#  way via polygon and arrows statements is a major PITA.  Don't mimic
#  the code below, instead do yourself a favor and use a package
#  designed for the task such as diagram, DiagrammeR, shape or Gmisc.
# Survival is a recommended package that is included by lots of others so
#  I try to limit dependencies in the survival vignettes.
#
par(mar=c(.1, .1, .1, .1))
frame()
oldpar <- par(usr=c(0,100,0,100))
# first figure
xx <- c(0, 10, 10, 0)
yy <- c(0, 0, 10, 10)
polygon(xx +10, yy+70)
polygon(xx +30, yy+70)
arrows( 22, 75, 28, 75, length=.1)
text(c(15, 35), c(75,75), c("Alive", "Dead"))

# second figure
polygon(xx +60, yy+70)  
for (j in c(55, 70, 85)) {
    polygon(xx +80, yy+j)
    arrows(72, (5*75 +j+5)/6, 78, (100+j*5)/6, length=.1)
}
text(c(65, 85,85,85), c(70,55,70,85)+5, c("A", "D3", "D2", "D1")) 

# third figure
polygon(xx+10, yy+25)
for (j in c(15,35)) {
    polygon(xx +30, yy+j)
    arrows(22, (5*30 +j+4)/6, 28, (54+j*5)/6, length=.1)
}
arrows(28, 2+(65 + 35*5)/6, 22, 2+ (160 + 40)/6, length=.1)
arrows(35, 33, 35, 27, length=.1)
text(c(15, 35,35), c(30, 20, 40), c("Health", "Death", "Illness"))

# fourth
for (i in c(50, 68)) polygon(xx+i, yy+25)
arrows(62, 30, 67, 30, length=.1)
arrows(80, 30, 84, 30, length=.1)
text(90, 30, "...", cex=2)
text(c(55, 73), c(30, 30), c("0", "1"))
par(oldpar)


###################################################
### code chunk number 3: crfig2
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mar=c(.1, .1, .1, .1))
frame()
oldpar <- par(usr=c(0,100,0,100))
# first figure
xx <- c(0, 10, 10, 0)
yy <- c(0, 0, 10, 10)
polygon(xx +10, yy+70) 
temp <- c(60, 80) 
for (j in 1:2) {
    polygon(xx + 30, yy+ temp[j])
    arrows(22, 70 + 3*j, 28, temp[j] +5, length=.1)
}
text(c(15, 35, 35), c(75, 65, 85),c("Entry", "Death", "PCM")) 
text(25, 55, "Competing Risk")

# Second figure
polygon(xx +60, yy+70)  
for (j in 1:2) {
    polygon(xx + 80, yy+ temp[j])
    arrows(72, 70+ 3*j, 78, temp[j] +5, length=.1)
}
text(50+ c(15, 35, 35), c(75, 65, 85),c("Entry", "Death", "PCM")) 
arrows(85, 78, 85, 72, length=.1)
text(75, 55, "Multi-state 1")

# third figure
polygon(xx+10, yy+25)
temp <- c(15, 35) 
for (j in 1:2) {
    polygon(2*xx +30, yy + temp[j])
    arrows(22, 25 + 3*j, 28, temp[j] +5, length=.1)
}
text(c(15, 40, 40), c(30, 20, 40),c("Entry", "Death w/o PCM", "PCM")) 
polygon(2*xx + 60, yy+temp[2])
arrows(52, 40, 58, 40, length=.1)
text(70, 40, "Death after PCM")
text(40, 10, "Multi-state 2")


###################################################
### code chunk number 4: mgus1
###################################################
getOption("SweaveHooks")[["fig"]]()
oldpar <- par(mfrow=c(1,2))
hist(mgus2$age, nclass=30, main='', xlab="Age")
with(mgus2, tapply(age, sex, mean))

mfit1 <- survfit(Surv(futime, death) ~ sex, data=mgus2)
mfit1
plot(mfit1, col=c(1,2), xscale=12, mark.time=FALSE, lwd=2,
     xlab="Years post diagnosis", ylab="Survival")
legend("topright", c("female", "male"), col=1:2, lwd=2, bty='n')
par(oldpar)


###################################################
### code chunk number 5: mgus2
###################################################
getOption("SweaveHooks")[["fig"]]()
etime <- with(mgus2, ifelse(pstat==0, futime, ptime))
event <- with(mgus2, ifelse(pstat==0, 2*death, 1))
event <- factor(event, 0:2, labels=c("censor", "pcm", "death"))
table(event)

mfit2 <- survfit(Surv(etime, event) ~ sex, data=mgus2)
print(mfit2, rmean=240, scale=12)
mfit2$transitions

plot(mfit2, col=c(1,2,1,2), lty=c(2,2,1,1),
     mark.time=FALSE, lwd=2,  xscale=12,
     xlab="Years post diagnosis", ylab="Probability in State")
legend(240, .6, c("death:female", "death:male", "pcm:female", "pcm:male"), 
       col=c(1,2,1,2), lty=c(1,1,2,2), lwd=2, bty='n')


###################################################
### code chunk number 6: mgus3
###################################################
getOption("SweaveHooks")[["fig"]]()
pcmbad <- survfit(Surv(etime, pstat) ~ sex, data=mgus2)
plot(pcmbad[2], mark.time=FALSE, lwd=2, fun="event", conf.int=FALSE, xscale=12,
     xlab="Years post diagnosis", ylab="Fraction with PCM")
lines(mfit2[2,1], lty=2, lwd=2, mark.time=FALSE, conf.int=FALSE)
legend(0, .25, c("Males, PCM, incorrect curve", "Males, PCM, competing risk"),
       col=1, lwd=2, lty=c(1,2), bty='n')


###################################################
### code chunk number 7: mgus4
###################################################
ptemp <- with(mgus2, ifelse(ptime==futime & pstat==1, ptime-.1, ptime))
newdata <- tmerge(mgus2, mgus2,  id=id, death=event(futime, death),
                  pcm = event(ptemp, pstat))
newdata <- tmerge(newdata, newdata, id, enum=cumtdc(tstart))
with(newdata, table(death, pcm))


###################################################
### code chunk number 8: mgus4g
###################################################
getOption("SweaveHooks")[["fig"]]()
temp <- with(newdata, ifelse(death==1, 2, pcm))
newdata$event <- factor(temp, 0:2, labels=c("censor", "pcm", "death"))  
mfit3 <- survfit(Surv(tstart, tstop, event) ~ sex, data=newdata, id=id)
print(mfit3, rmean=240, digits=2)
mfit3$transitions
plot(mfit3[,1], mark.time=FALSE, col=1:2, lty=1:2, lwd=2,
     xscale=12,
     xlab="Years post MGUS diagnosis", ylab="Prevalence of PCM")
legend(48, .04, c("female", "male"), lty=1:2, col=1:2, lwd=2, bty='n') 


###################################################
### code chunk number 9: mgus5
###################################################
getOption("SweaveHooks")[["fig"]]()
# Death after PCM will correspond to data rows with
#  enum = 2 and event = death 
d2 <- with(newdata, ifelse(enum==2 & event=='death', 4, as.numeric(event)))
e2 <- factor(d2, labels=c("censor", "pcm", "death w/o pcm", 
                          "death after pcm"))
mfit4 <- survfit(Surv(tstart, tstop, e2) ~ sex, data=newdata, id=id)
plot(mfit2[2,], lty=c(1,2),
     xscale=12, mark.time=FALSE, lwd=2, 
     xlab="Years post diagnosis", ylab="Prevalence")
lines(mfit4[2,3], mark.time=FALSE, col=2, lty=1, lwd=2,
      conf.int=FALSE)

legend(200, .5, c("Death w/o PCM", "ever PCM", 
                  "Death after PCM"), col=c(1,1,2), lty=c(2,1,1), 
             lwd=2, bty='n', cex=.82)


###################################################
### code chunk number 10: cfit1
###################################################
options(show.signif.stars = FALSE)  # display intelligence
cfit2 <- coxph(Surv(etime, event=="death") ~ age + sex + mspike, mgus2)
summary(cfit2, scale=c(10, 1, 1))   # scale age in decades 



###################################################
### code chunk number 11: cfit2
###################################################
cfit1 <- coxph(Surv(etime, event=="pcm") ~ age + sex + mspike, mgus2)
cfit1
quantile(mgus2$mspike, na.rm=TRUE)


###################################################
### code chunk number 12: mpyears
###################################################
pfit1 <- pyears(Surv(ptime, pstat) ~ sex, mgus2, scale=12)
round(100* pfit1$event/pfit1$pyears, 1)  # PCM rate per year

temp <- summary(mfit1, rmean="common")  #print the mean survival time
round(temp$table[,1:6], 1)


###################################################
### code chunk number 13: PCMcurve
###################################################
getOption("SweaveHooks")[["fig"]]()
newdata <- expand.grid(sex=c("F", "M"), age=c(60, 80), mspike=1.2)
newdata

temp <- matrix(list(), 3,3)
dimnames(temp) <- list(from=c("Entry", "PCM", "Death"),
                       to  =c("Entry", "PCM", "Death"))
temp[1,2] <- list(survfit(cfit1, newdata, std.err=FALSE))
temp[1,3] <- list(survfit(cfit2, newdata, std.err=FALSE))
csurv  <- survfit(temp, p0 =c(1,0,0))
plot(csurv[,2], xmax=25*12, xscale=12, 
     xlab="Years after MGUS diagnosis", ylab="PCM",
     col=1:2, lty=c(1,1,2,2), lwd=2)
legend(10, .14, outer(c("female", "male   "), 
                     c("diagnosis at age 60", "diagnosis at age 80"), 
                      paste, sep=", "),
       col=1:2, lty=c(1,1,2,2), bty='n', lwd=2)


###################################################
### code chunk number 14: year20
###################################################
# Print out a M/F results at 20 years
temp <- summary(csurv, time=20*12)$prev
cbind(newdata, PCM= round(100*temp[,2], 1))


###################################################
### code chunk number 15: fg1
###################################################
# (first three lines are identical to an earlier section)
etime <- with(mgus2, ifelse(pstat==0, futime, ptime))
event <- with(mgus2, ifelse(pstat==0, 2*death, 1))
event <- factor(event, 0:2, labels=c("censor", "pcm", "death"))

pcmdat <- finegray(Surv(etime, event) ~ ., data=mgus2,
                   etype="pcm")
pcmdat[1:4, c(1:3, 11:14)]

deathdat <- finegray(Surv(etime, event) ~ ., data=mgus2,
                     etype="death")
dim(pcmdat)
dim(deathdat)
dim(mgus2)


###################################################
### code chunk number 16: pfit2
###################################################
# The PCM curves of the multi-state model
pfit2 <- survfit(Surv(fgstart, fgstop, fgstatus) ~ sex,
                data=pcmdat, weight=fgwt)
# The death curves of the multi-state model
dfit2 <- survfit(Surv(fgstart, fgstop, fgstatus) ~ sex, 
                  data=deathdat, weight=fgwt)


###################################################
### code chunk number 17: fg2
###################################################
getOption("SweaveHooks")[["fig"]]()
fgfit1 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ sex, data=pcmdat,
               weight= fgwt)
summary(fgfit1)
fgfit2 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ sex, data=deathdat,
               weight= fgwt)
fgfit2

mfit2 <- survfit(Surv(etime, event) ~ sex, data=mgus2) #reprise the AJ
plot(mfit2[,1], col=1:2,
     lwd=2,  xscale=12,
     conf.times=c(5, 15, 25)*12,
     xlab="Years post diagnosis", ylab="Fraction with PCM")
ndata <- data.frame(sex=c("F", "M"))
fgsurv1 <- survfit(fgfit1, ndata)
lines(fgsurv1, fun="event", lty=2, lwd=2, col=1:2)
legend("topleft", c("Female, Aalen-Johansen", "Male, Aalen-Johansen",
                 "Female, Fine-Gray", "Male, Fine-Gray"),
       col=1:2, lty=c(1,1,2,2), bty='n')

# rate models with only sex
pfitr <- coxph(Surv(etime, event=="pcm") ~ sex, mgus2)
dfitr <- coxph(Surv(etime, event=="death") ~ sex, mgus2)
temp <- matrix(list(), 3,3)
temp[1,2] <- list(survfit(pfitr, ndata, std.err=FALSE))
temp[1,3] <- list(survfit(dfitr, ndata, std.err=FALSE))
rcurve <- survfit(temp, p0=c(entry=1, pcm=0, death=0))


###################################################
### code chunk number 18: fg3
###################################################
fgfit2a <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age + sex + mspike,
                 data=pcmdat, weight=fgwt)

fgfit2b <-  coxph(Surv(fgstart, fgstop, fgstatus) ~ age + sex + mspike,
                 data=deathdat, weight=fgwt)
round(rbind(PCM= coef(fgfit2a), death=coef(fgfit2b)), 3)


###################################################
### code chunk number 19: finegray2
###################################################
getOption("SweaveHooks")[["fig"]]()
oldpar <- par(mfrow=c(1,2))
newdata <- expand.grid(sex= c("F", "M"), age=c(60, 80), mspike=1.2)
fsurv1 <- survfit(fgfit2a, newdata)  # time to progression curves
plot(fsurv1, xscale=12, col=1:2, lty=c(1,1,2,2), lwd=2, fun='event',
     xlab="Years", ylab="Fine-Gray predicted", 
     xmax=12*25, ylim=c(0, .15))
legend(1, .15, c("Female, 60", "Male, 60","Female: 80", "Male, 80"),
       col=c(1,2,1,2), lty=c(1,1,2,2), lwd=2, bty='n')

plot(csurv[,2], xscale=12, col=1:2, lty=c(1,1,2,2), lwd=2,
     xlab="Years", ylab="Multi-state predicted", 
     xmax=12*25, ylim=c(0, .15))
legend(1, .15, c("Female, 60", "Male, 60","Female: 80", "Male, 80"),
       col=c(1,2,1,2), lty=c(1,1,2,2), lwd=2, bty='n')
par(oldpar)


###################################################
### code chunk number 20: finegray-check
###################################################
getOption("SweaveHooks")[["fig"]]()
zph.fgfit2a <- cox.zph(fgfit2a)
zph.fgfit2a
plot(zph.fgfit2a[1])
abline(h=coef(fgfit2a)[1], lty=2, col=2)


###################################################
### code chunk number 21: finegray3
###################################################
getOption("SweaveHooks")[["fig"]]()
fsurv2 <- survfit(fgfit2b, newdata)  # time to progression curves
xtime <- 0:(30*12)  #30 years
y1a <- 1 - summary(fsurv1, times=xtime)$surv  #predicted pcm
y1b <- 1 - summary(fsurv2, times=xtime)$surv #predicted deaths before pcm
y1  <- (y1a + y1b)  #either

matplot(xtime/12, y1, col=1:2, lty=c(1,1,2,2), type='l',
        xlab="Years post diagnosis", ylab="FG: either endpoint")
abline(h=1, col=3)
legend("bottomright", c("Female, 60", "Male, 60","Female: 80", "Male, 80"),
       col=c(1,2,1,2), lty=c(1,1,2,2), lwd=2, bty='n')


###################################################
### code chunk number 22: pcmstack
###################################################
temp1 <- data.frame(mgus2, time=etime, status=(event=="pcm"), group='pcm')
temp2 <- data.frame(mgus2, time=etime, status=(event=="death"), group="death")
stacked <- rbind(temp1, temp2)
allfit <- coxph(Surv(time, status) ~ hgb + (age + sex)*strata(group),
                 data=stacked)


###################################################
### code chunk number 23: compete.Rnw:1091-1094
###################################################
if (file.exists("mstate.rda")) load("mstate.rda") else {
    stop("no local copy of mstate data was found")
}


###################################################
### code chunk number 24: aids1
###################################################
getOption("SweaveHooks")[["fig"]]()
adata <- aidssi
adata$event <- factor(adata$status, 0:2, c("censored", "AIDS", "SI"))

# KM curves that censor the other endpoint (a bad idea)
bad1 <- survfit(Surv(time, event=="AIDS") ~ 1, adata)
bad2 <- survfit(Surv(time, event=="SI") ~1, adata)

# The correct Aalen-Johansen curves
ajfit <- survfit(Surv(time, event) ~1, adata)
ajfit$transitions
plot(ajfit, xmax=13, col=1:2, lwd=2,
     xlab="Years from HIV infection", ylab="Probability")
legend(8, .2, c("AIDS", "SI"), lty=1, lwd=2, col=1:2, bty='n')


###################################################
### code chunk number 25: fitT2
###################################################
getOption("SweaveHooks")[["fig"]]()
# re-create figure T2
plot(bad1, conf.int=FALSE, xmax=13, 
     xlab="Years from HIV infection", ylab="Probability")
lines(bad2, conf.int=FALSE, fun='event', xmax=13)
text(c(8,8), c(.8, .22), c("AIDS", "SI"))


###################################################
### code chunk number 26: figT3
###################################################
getOption("SweaveHooks")[["fig"]]()
# re-create figure T3
plot(bad1, conf.int=FALSE, xmax=13, col="lightgray",
     xlab="Years from HIV infection", ylab="Probability")
lines(bad2, conf.int=FALSE, fun='event', xmax=13, col='lightgray')
text(c(8,8), c(.8, .22), c("AIDS", "SI"))

lines(ajfit[,2], conf.int=FALSE, lwd=2, xmax=13)
lines(ajfit[,1], conf.int=FALSE, lwd=2, ,xmax=13, fun=function(x) x)


###################################################
### code chunk number 27: figT4
###################################################
getOption("SweaveHooks")[["fig"]]()
# re-create figure T4
temp <- ajfit
temp$prev <- t(apply(temp$prev, 1, cumsum))  # apply() transposes
plot(temp, xmax=13, lwd=2, col=1, ylim=c(0,1), 
        xlab="Years from HIV infection", ylab="Probability")
text(c(11, 11, 11), c(.2, .55, .9), c("AIDS", "SI", "Event free"))


###################################################
### code chunk number 28: cfit
###################################################
cfit1 <- coxph(Surv(time, event=="AIDS") ~ ccr5, adata)
cfit1

cfit2 <- coxph(Surv(time, event=="SI") ~ ccr5, adata)
cfit2


###################################################
### code chunk number 29: stack
###################################################
temp <- subset(adata, select=c(time, ccr5))
temp1 <- data.frame(temp, status= 1*(adata$event=="AIDS"), cause="AIDS")
temp2 <- data.frame(temp, status= 1*(adata$event=="SI"),   cause="SI")
stack <- rbind(temp1, temp2)

cfit3 <- coxph(Surv(time, status) ~ ccr5 * strata(cause), data=stack)
cfit3
sum(coef(cfit3))


###################################################
### code chunk number 30: stack2
###################################################
stack$ccr5.1 <- (stack$ccr5=="WM") * (stack$cause == "AIDS")
stack$ccr5.2 <- (stack$ccr5=="WM") * (stack$cause == "SI")
coxph(Surv(time, status) ~ ccr5.1 + ccr5.2 + strata(cause), stack)


###################################################
### code chunk number 31: cfit4
###################################################
cfit4 <- coxph(Surv(time, status) ~ ccr5.1 + ccr5.2 + cause, stack)
cfit4
summary(cfit4)$coefficients
cfit4b <-  coxph(Surv(time, status) ~ ccr5*cause, stack) #same result


###################################################
### code chunk number 32: figT5
###################################################
getOption("SweaveHooks")[["fig"]]()
# re-create figure T5 in a single panel
tdata <- expand.grid(ccr5=c("WW","WM"), cause=c("AIDS", "SI"))
tdata

tsurv <- survfit(cfit4b, newdata=tdata)
smat <- matrix(list(), 3, 3,
               dimnames=list(from= c("AIDS", "SI", "entry"),
                             to =  c("AIDS", "SI", "entry")))
smat[3,1] <- list(tsurv[1:2])
smat[3,2] <- list(tsurv[3:4])
smat   #did we put things in the right place?

csurv1 <- survfit(smat, p0=c(0,0,1))
plot(csurv1[,1:2], col=1:2, lty=c(1,1,2,2), xmax=13, lwd=2, 
     xlab="Years from HIV infection", ylab="Probability")
legend(0, .4, c("AIDS, WW", "AIDS, WM", "SI, WW", "SI, WM"),
       col=1:2, lty=c(1,1,2,2), lwd=2, bty='n')


###################################################
### code chunk number 33: compete.Rnw:1242-1251
###################################################
fdata1 <- finegray(Surv(time, event) ~ ., adata, etype='AIDS')
fgfit1 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ ccr5, fdata1,
                weight = fgwt)
fgfit1

fdata2 <- finegray(Surv(time, event) ~., adata, etype="SI")
fgfit2 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ ccr5, fdata2,
                weight = fgwt)
fgfit2


###################################################
### code chunk number 34: figT8
###################################################
getOption("SweaveHooks")[["fig"]]()
# re-create figure T8: Fine-Gray curves
fgsurv1<-survfit(fgfit1,newdata=tdata)
fgsurv2<-survfit(fgfit2,newdata=tdata)

oldpar <- par(mfrow=c(1,2), mar=c(4.1, 3.1, 3.1, 1)) #leave room for title
plot(fgsurv1, col=1:2, lty=c(1,1,2,2), lwd=2, xmax=13,
     ylim=c(0, .5),fun='event',
     xlab="Years from HIV infection", ylab="Probability")
title("AIDS")
plot(fgsurv2, col=1:2, lty=c(1,1,2,2), lwd=2, xmax=13,
     ylim=c(0, .5), fun='event',
     xlab="Years from HIV infection", ylab="Probability")
title("SI appearance")     
par(oldpar)


###################################################
### code chunk number 35: figT9
###################################################
getOption("SweaveHooks")[["fig"]]()
# re-create figure T9: curves by CCR type
aj2 <- survfit(Surv(time, event) ~ ccr5, adata)
oldpar <- par(mfrow=c(1,2))
plot(aj2[,1], xmax=13, col=1:2, lwd=2, ylim=c(0, .5),
     xlab="Years from HIV infection", ylab="Probability of AIDS")
text(c(10, 10), c(.35, .07), c("WW", "WM"))

plot(aj2[,2], xmax=13, col=1:2, lwd=2, ylim=c(0, .5), 
     xlab="Years from HIV infection", ylab="Probability of SI")
text(c(8, 8), c(.34, .18), c("WW", "WM"))
par(oldpar)


###################################################
### code chunk number 36: tableT2
###################################################
table(ebmt3$dissub)
table(ebmt3$drmatch)
table(ebmt3$tcd)
table(ebmt3$age)


###################################################
### code chunk number 37: data1
###################################################
getOption("SweaveHooks")[["fig"]]()
temp <- subset(ebmt3, select = -c(prtime, prstat, rfstime, rfsstat))
edata <- tmerge(temp, ebmt3, id, 
                rstat = event(rfstime, rfsstat),
                pstat = event(prtime, prstat),
                enum  = tdc(prtime))
print(edata[15:20,-(3:5)])

# Check that no one has recovery and death on the same day
with(edata, table(rstat, pstat))

# Create the factor outcome
edata$event <- with(edata, factor(pstat + 2*rstat, 0:2,
                           labels = c("censor", "PR", "RelDeath")))
levels(edata$drmatch) <- c("Match", "Mismatch")

surv1 <- survfit(Surv(tstart, tstop, event) ~ 1, edata, id=id)
surv1$transitions   # matches the Frequencies on page C5
plot(surv1, col=1:2, xscale=365.25, lwd=2, 
     xlab="Years since transplant", ylab="Fraction in state")
legend(1000, .2, c("Platelet recovery", "Death or Relapse"), 
       lty=1, col=1:2, lwd=2, bty='n')


###################################################
### code chunk number 38: data2
###################################################
temp1 <- with(edata, data.frame(edata[enum==0,], status=pstat[enum==0], 
                 trans="1->2", from=1, to=2))  # baseline to PR
temp2 <- with(edata, data.frame(edata[enum==0,], status=rstat[enum==0], 
                 trans="1->3", from=1, to=3)) # baseline to relapse/death
temp3 <- with(edata, data.frame(edata[enum==1,], status=rstat[enum==1], 
                 trans="2->3", from=2, to=3)) # PR to replase/death
edata2 <- rbind(temp1, temp2, temp3)  # the stacked data set


###################################################
### code chunk number 39: efit1
###################################################
efit1.2 <- coxph(Surv(tstop, event=='PR') ~ 
                    dissub + age + drmatch + tcd, ties='breslow',
                    data=edata, subset = (enum==0))
efit1.3 <- coxph(Surv(tstop, event=='RelDeath') ~ 
                    dissub + age + drmatch + tcd, ties= 'breslow',
                    data=edata, subset = (enum==0))
efit2.3 <- coxph(Surv(tstart, tstop, event=='RelDeath') ~ 
                    dissub + age + drmatch + tcd, ties='breslow',
                    data=edata, subset = (enum==1))
round(cbind('1->2'= coef(efit1.2) , '1->3'= coef(efit1.3), 
            '2->3' =coef(efit2.3)), 3)


###################################################
### code chunk number 40: compete.Rnw:1390-1405
###################################################
efit1b <- coxph(Surv(tstart, tstop, status) ~ (dissub + age + drmatch + tcd)
                           *strata(trans), data=edata2, ties="breslow")
efit1c <- coxph(Surv(tstart, tstop, status) ~ 
               strata(trans)/(dissub + age + drmatch + tcd),
               data=edata2, ties="breslow")
# Rearrange the order so as to match the paper
mycoef <- summary(efit1c)$coefficients[, c(1,3)]  # coef and se
index <- as.vector(matrix(1:18, ncol=3, byrow=T))
round(mycoef[index,], 3)
    
# Compare the log-likelihoods of the fits          
matrix(c(efit1.2$loglik + efit1.3$loglik + efit2.3$loglik, 
         efit1b$loglik, efit1c$loglik), nrow=2,
       dimnames=list(c("Initial LL", "Final LL"), 
                     c("fit a", "fit b", "fit c")))


###################################################
### code chunk number 41: efit2
###################################################
efit2a <- coxph(Surv(tstart, tstop, status) ~ 
                 factor(from)*(dissub + age + drmatch + tcd),
                data=edata2, ties="breslow", subset= (to==3))
efit2b <- coxph(Surv(tstart, tstop, status) ~ 
                 factor(from)/(dissub + age + drmatch + tcd),
                data=edata2, ties="breslow", subset= (to==3))
mycoef <- summary(efit2a)$coefficients[, c(1,3)]
round(mycoef, 3)[c(2,4,6,8,10,12, 3,5,7,9,11,13,1),]

efit2c <- coxph(Surv(tstart, tstop, status) ~ strata(to) +
               trans/(dissub + age + drmatch + tcd),
               data=edata2, ties="breslow")
round(summary(efit2b)$coefficients[, c(1,3)], 3)


###################################################
### code chunk number 42: zph
###################################################
getOption("SweaveHooks")[["fig"]]()
z2 <- cox.zph(efit2a, transform=function(x) sqrt(x))
plot(z2[7], resid=FALSE)
abline(h=coef(efit2a)[7], col=2)
z2[7]


###################################################
### code chunk number 43: tableS4
###################################################
# The extra variable for column 3 of table 3
efit3 <- coxph(Surv(tstart, tstop, status) ~ I(tstart/365.25) +
                factor(from)/(dissub + age + drmatch + tcd),
                data=edata2, ties="breslow", subset= (to==3)) 
coef(efit3)["I(tstart/365.25)"]


###################################################
### code chunk number 44: tableT4
###################################################
efit4 <- coxph(Surv(tstop - tstart, status) ~ 
               strata(trans)/(dissub + age + drmatch + tcd),
               data=edata2, ties="breslow")
mycoef4 <- summary(efit4)$coefficients[, c(1,3)]  # coef and se
index4 <- index[7:18]  # leave off transition 1->2
round(mycoef4[index4,], 3)  


###################################################
### code chunk number 45: compete.Rnw:1518-1538
###################################################
getOption("SweaveHooks")[["fig"]]()
newdata1 <- expand.grid(age="<=20", dissub="AML", drmatch="Mismatch",
                       tcd=c("No TCD", "TCD"))
newdata2 <- cbind(newdata1[c(1,2,1,2),], from=c(1,1,2,2))
newdata2

tcurve1 <- survfit(efit1.2, newdata1, se.fit=FALSE)
tcurve2 <- survfit(efit2a,  newdata2, se.fit=FALSE)

tmat <- matrix(list(), 3, 3,
               dimnames=list(from=c("Tx", "PR", "R/D"),
                             to  =c("Tx", "PR", "R/D")))
tmat[1,2] <- list(tcurve1)
tmat[1,3] <- list(tcurve2[1:2] )
tmat[2,3] <- list(tcurve2[3:4] ) 
ecurve <- survfit(tmat, p0=c(1,0,0))
plot(ecurve, col=c(1,1,2,2,3,3), lty=1:2, lwd=2, xscale=365.25,
     xlab="Years since transplant", ylab="Predicted probabilities")
legend(700, .9, c("Currently alive in remission, no PR", "Currently in PR",
               "Relapse or death"), col=1:3, lwd=2, bty='n')
text(700, .95, "Solid= No TCD, dashed = TCD", adj=0)


###################################################
### code chunk number 46: fourstate
###################################################
getOption("SweaveHooks")[["fig"]]()
dtemp <- c("PR", "R/D after PR", "R/D without PR", "Tx")
xmat <- matrix(list(), 4, 4, dimnames=list(from=dtemp, to=dtemp))
xmat[4,1] <- list(tcurve1)
xmat[4,3] <- list(tcurve2[1:2] )
xmat[1,2] <- list(tcurve2[3:4] ) 
ecurve2 <- survfit(xmat, p0=c(0,0,0,1))
dim(ecurve2)   # rows are the two "subjects", cols are the states
plot(ecurve2[, 1:3], lwd=2, lty=1:2, col=c(2,2,3,3,4,4,1,1), 
     xscale=365.25,
     xlab="Years since transplant", ylab="Predicted probabilities")
text(c(1500, 1500, 1500), c(.55, .3, .05), 
     c("PR", "R/D after PR", "R/D w/o PR"), col=2:4)
text(1550, .66, "solid= no SCD, dashed = SCD", col='gray')


###################################################
### code chunk number 47: figT14
###################################################
getOption("SweaveHooks")[["fig"]]()
tcurve <- ecurve2
tcurve$prev <- t(apply(tcurve$prev, 1, cumsum))
oldpar <- par(mfrow=c(1,2), mar=c(4.1, 3.1, 3.1, .1))
plot(tcurve[1,1:4], col=1, xscale=365.25, ylim=c(0,1), 
      xlab="Years since transplant", ylab="Predicted probabilities")
text(rep(4*365, 4), c(.4, .55, .7, .95), cex=.7, 
     c("Alive in remission, PR", "Relapse or death after PR",
       "Relapse or death without PR", "Alive in remission, no PR"))
title("No TCD")
plot(tcurve[2,1:4], col=1, xscale=365.25, ylim=c(0,1), 
      xlab="Years since transplant", ylab="Predicted probabilities")
text(rep(4*365, 4), c(.45, .7, .82, .95), cex=.7, 
     c("Alive in remission, PR", "Relapse or death after PR",
       "Relapse or death without PR", "Alive in remission, no PR"))
title("TCD")
par(oldpar)


