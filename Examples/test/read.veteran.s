#       Treatment  1=standard,  2=test
#       Celltype   1=squamous,  2=smallcell,  3=adeno,  4=large
#       Survival in days
#       Status     1=dead, 0=censored
#       Karnofsky score
#       Months from Diagnosis
#       Age in years
#       Prior therapy  0=no, 10=yes

veteran <- read.table('veteran.data', col.names=c("rx", "celltype",
			'futime', 'fustat', 'karno', 'months', 'age',
			'prior'))

veteran$celltype <- factor(veteran$celltype, levels=1:4,
		labels=c('squamous', 'smallcell', 'adeno', 'large'))
veteran$prior <- veteran$prior/10
