#
# Do further tests of parametric regression models
#
 cat setup.s > temp
 cat capacitor.s >> temp
 cat donnell.s donnell2.s  >> temp
 cat lung.s lung2.s >> temp
 cat peterson.s >> temp
 cat ptest.s >> temp
# cat resid.s >> temp
 cat stanford.s stanford2.s >> temp
 cat testscale.s >> temp
 cat strata.s >> temp
 cat testrat.s >> temp
 cat mpip.s    >> temp
 echo 'q()' >> temp
 rm -rf .Data
 mkdir .Data
 Splus BATCH temp testall.out
