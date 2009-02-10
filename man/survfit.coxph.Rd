\name{survfit.coxph}
\alias{survfit.coxph}
\title{
  Compute a Survival Curve from a Cox model
}
\description{
Computes the predicted survivor function for a Cox proportional 
hazards model. 
}
\usage{
\method{survfit}{coxph}(formula, newdata, 
        se.fit=TRUE, conf.int=.95,
        individual=FALSE, 
        type,vartype,
        conf.type=c("log","log-log","plain","none"),...)
}
\arguments{
  \item{formula}{
    A \code{coxph} object. 
  }
  \item{newdata}{
    a data frame with the same variable names as those that appear 
    in the \code{coxph} formula.  
    
    The curve(s) produced will be representative of a cohort whose 
    covariates correspond to the values in \code{newdata}. 
    Default is the mean of the covariates used in the 
    \code{coxph} fit. 
  }
  \item{individual}{
    a logical value indicating whether each row of \code{newdata} represents
    a distinct individual (FALSE, the default), 
    or if each row of the data frame represents different 
    time epochs for only one individual (TRUE).  
    In the former case the result will have one curve for each row in
    \code{newdata}, in the latter only a single curve will be produced.
  }
  \item{conf.int}{
    the level for a two-sided confidence interval on the survival curve(s). 
    Default is 0.95. 
  }
  \item{se.fit}{
    a logical value indicating whether standard errors should be 
    computed.  Default is \code{TRUE}. 
  }
  \item{type,vartype}{
    a character string specifying the type of survival curve. 
    Possible values are  \code{"aalen"} or \code{"kaplan-meier"} 
    (only the first two characters are necessary). 
    The default is \code{"aalen"}  
    Earlier versions of \code{survfit} 
    used \code{type="tsiatis"} to get the \code{"aalen"} estimator.  For backward 
    compatibility, this is still allowed. 
  }
 
  \item{conf.type}{
    One of \code{"none"}, \code{"plain"}, \code{"log"} (the default),
    or \code{"log-log"}.  Only
    enough of the string to uniquely identify it is necessary.
    The first option causes confidence intervals not to be
    generated.  The second causes the standard intervals
    \code{curve +- k *se(curve)}, where k is determined from
    \code{conf.int}.  The log option calculates intervals based on the
    cumulative hazard or log(survival). The last option bases
    intervals on the log hazard or log(-log(survival)). 
  }
  \item{\dots}{for future methods}
}
\value{
an object of class \code{"survfit"}.  
See \code{survfit.object} for 
details. Methods defined for survfit objects are  
\code{print}, \code{plot}, 
\code{lines}, and \code{points}. 
}
\details{
See \code{\link{survfit}} for more Details.
  
The estimates used are the Kalbfleisch-Prentice 
(Kalbfleisch and Prentice, 1980, p.86) and the Tsiatis/Link/Breslow, 
which reduce to the Kaplan-Meier and Fleming-Harrington estimates, 
respectively, when the weights are unity.  
When the parent data is a Cox model, there is an extra
term in the variance of the curve, due to the variance of the coefficients and 
hence variance in the computed weights.  Details of the Aalen 
estimator and its variance are found in Tsiatis (1981). 

The two variants of the F-H estimate have to do with how ties are handled.
If there were 3 deaths out of 10 at risk, then the first
increments the hazard by 3/10 and the second
by 1/10 + 1/9 + 1/8.
For the first method S(t) = exp(H), where H is 
the Nelson-Aalen cumulative hazard estimate,
whereas the \code{fh2} method will 
give results S(t) results closer to the Kaplan-Meier. 
For curves created 
after a Cox model these correspond to the Breslow and Efron estimates, 
respectively, and the proper choice is made automatically. 

}
\section{References}{
Dorey, F. J. and Korn, E. L. (1987).  Effective sample sizes for confidence 
intervals for survival probabilities.  \emph{Statistics in Medicine} 
\bold{6}, 679-87. 

Fleming, T. H. and Harrington, D. P. (1984).  Nonparametric estimation of the 
survival distribution in censored data.  \emph{Comm. in Statistics}  
\bold{13}, 2469-86. 

Kablfleisch, J. D. and Prentice, R. L. (1980).  
\emph{The Statistical Analysis of Failure Time Data.}
New York:Wiley. 

Link, C. L. (1984). Confidence intervals for the survival 
function using Cox's proportional hazards model with  
covariates.  \emph{Biometrics}  
\bold{40}, 601-610.

Tsiatis, A. (1981). A large sample study of the estimate 
for the integrated hazard function in Cox's regression 
model for survival data. \emph{Annals of Statistics}  
\bold{9}, 93-108. 
}
\seealso{
\code{\link{print}},  
\code{\link{plot}},  
\code{\link{lines}},   
\code{\link{coxph}},  
\code{\link{Surv}},  
\code{\link{strata}}.   
}
\examples{
#fit a Kaplan-Meier and plot it 
fit <- survfit(Surv(time, status) ~ x, data = aml) 
plot(fit, lty = 2:3) 
legend(100, .8, c("Maintained", "Nonmaintained"), lty = 2:3) 

#fit a Cox proportional hazards model and plot the  
#predicted survival for a 60 year old 
fit <- coxph(Surv(futime, fustat) ~ age, data = ovarian) 
plot(survfit(fit, newdata=data.frame(age=60)),
     xscale=365.25, xlab = "Years", ylab="Survival") 

# Here is the data set from Turnbull
#  There are no interval censored subjects, only left-censored (status=3),
#  right-censored (status 0) and observed events (status 1)
#
#                             Time
#                         1    2   3   4
# Type of observation
#           death        12    6   2   3
#          losses         3    2   0   3
#      late entry         2    4   2   5
#
tdata <- data.frame(time  =c(1,1,1,2,2,2,3,3,3,4,4,4),
                    status=rep(c(1,0,2),4),
                    n     =c(12,3,2,6,2,4,2,0,2,3,3,5))
fit  <- survfit(Surv(time, time, status, type='interval') ~1, 
              data=tdata, weight=n)

#
# Time to progression/death for patients with monoclonal gammopathy
#  Competing risk curves (cumulative incidence)
fit1 <- survfit(Surv(stop, event=='progression') ~1, data=mgus1,
                    subset=(start==0))
fit2 <- survfit(Surv(stop, status) ~1, data=mgus1,
                    subset=(start==0), etype=event) #competing risks
# CI curves are always plotted from 0 upwards, rather than 1 down
plot(fit2, fun='event', xscale=365.25, xmax=7300, mark.time=FALSE,
            col=2:3, xlab="Years post diagnosis of MGUS")
lines(fit1, fun='event', xscale=365.25, xmax=7300, mark.time=FALSE,
            conf.int=FALSE)
text(10, .4, "Competing Risk: death", col=3)
text(16, .15,"Competing Risk: progression", col=2)
text(15, .30,"KM:prog")
}
\keyword{survival}
