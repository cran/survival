library(survival)
#
# Why does coxph lose label attributes?
#
adata <- aml
attr(adata$x, "label") <- "Group"
test <- na.omit(adata)
attributes(test$x)  # here's the problem, na.action loses labels

# The copyMostAttribNoTs function was added to xtra.R, called by coxph
fit <- coxph(Surv(time, status) ~ x, adata, model=TRUE)
attr(fit$model$x, "label")== "Group"

# Now try a Surv2 model
bdata <- rbind(data.frame(id=1:23, time=0, status=0, x=aml$x),
               data.frame(id=1:23, time=aml$time, status=aml$status, x=NA))
bdata <- bdata[order(bdata$id, bdata$time),]
attr(bdata$x, "label") <- "Group"

# also called by surv2counting
test2 <- fromtimeline(Surv2(time, status) ~., bdata, id=id)
attr(test2$x, "label") == "Group"

fit2 <- coxph(Surv2(time, status) ~x, bdata, id=id, model=TRUE)
attr(fit2$model$x, "label") == "Group"
