temp _ scan("data.jasa", what=list(id=0, b.mo=0, b.d=0, b.y=0,
				       a.mo=0, a.d=0, a.y=0,
				       t.mo=0, t.d=0, t.y=0,
				       f.mo=0, f.d=0, f.y=0,
				       fu.stat=0, surg=0, mismatch=0,
				       hla.a2=0, mscore=0, reject=0))

temp3 _ mdy.date(temp$b.mo, temp$b.d, temp$b.y)
temp4 _ mdy.date(temp$a.mo, temp$a.d, temp$a.y)
temp5 _ mdy.date(temp$t.mo, temp$t.d, temp$t.y)
temp6 _ mdy.date(temp$f.mo, temp$f.d, temp$f.y)

# Make sure that a particular idiocy is turned off: turning logicals into 
#	factors!
as.data.frame.logical <- as.data.frame.vector
jasa <- data.frame( birth.dt=temp3, accept.dt=temp4, tx.date=temp5,
		    fu.date=temp6, fustat=temp$fu.stat,
		    surgery = temp$surg,  age= temp4-temp3,
		    futime = 1+temp6 - temp4, wait.time= 1+temp5 - temp4,
		    transplant = c(!is.na(temp5)),
		    mismatch=temp$mismatch, hla.a2=temp$hla.a2,
		    mscore = temp$mscore, reject=temp$reject)

row.names(jasa) <- temp$id

# The "1+" above causes us to match the analysis in Kalbfleisch and Prentice.
#  Someone accepted and transplanted on the same day is assumed to have one
#  day under the non-transplant risk.

rm(temp, temp3, temp4, temp5, temp6)




