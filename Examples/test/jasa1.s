expect <- survexp(futime ~ ratetable(age=(accept.dt - birth.dt), sex=1,
		year=accept.dt), jasa, cohort=F, ratetable=survexp.uswhite)

survdiff(Surv(jasa$futime, jasa$fustat) ~ offset(expect))
