#temp <- sas.get("../../../../data/moertel/sasdata", "anal")
#colon <- temp[temp$study==1,]
#rm(temp)
#colon$rx <- factor(colon$rx, levels=1:3, labels=c("Obs", "Lev", "Lev+5FU"))

data.restore('data.colon')
