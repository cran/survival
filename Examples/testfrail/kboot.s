# Bootstrap the kidney data
#    Use the fact that the ids are 1,1 2,2,3,3,4,.....

nboot <- 100
ngrp <- length(unique(kidney$id))

ksig <- double(nboot)
idlist <- unique(kidney$id)
for (i in 1:nboot) {
    idx <- sample(idlist, size=ngrp, replace=T)
    keep <- c(idx*2 -1, idx*2)
    idnew <- rep(1:ngrp, 2)

    fit <- coxph(Surv(time, status) ~ age + sex + frailty(idnew),
		 kidney[keep,])
    ksig[i] <- fit$history[[1]]$theta
cat(i," ")
    }


