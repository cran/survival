attach("../.Data")
dyn.load('../loadmod.o')
postscript(file='testall.ps')
options(na.action='na.omit', contrasts='contr.treatment')
