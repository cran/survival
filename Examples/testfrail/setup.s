#
# Set up for the test
#
dyn.load("../loadmod.o")
attach("../.Data")
options(na.action="na.omit", contrasts='contr.treatment')
