#
# Do the full shooting match of tests
#
 cat setup.s > temp
 cat simple.s  >> temp
 cat rfit.s >> temp
 cat cancer.s >> temp
 cat ovarian.s  ovarian2.s >> temp
 cat kfit.s kfit2.s >> temp
 cat readcolon.s colon.s >> temp
 cat resid.s resid.null.s >> temp
 cat rat2.s >> temp
 echo 'q()' >> temp
 rm .Data/*
 Splus BATCH temp testall.out
