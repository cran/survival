bladder <- read.table('data.bladder4',
	     col.names=c('id', 'rx', 'number', 'size', 'start', 'stop',
				'event', 'enum'))

bladder2 <- bladder[bladder$start< bladder$stop, ]
bladder$start <- NULL
bladder <- bladder[bladder$enum<5,]
