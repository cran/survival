# A subset of some local data on lung cancer patients.  A useful internal test
#  because it has multiple strata and lots of missing values.
#

# inst = enrolling institution
# sex  1=male  2=female
# ph.ecog  physician's estimate of the ECOG performace score.  0=fully active,
#              4=bedridden
# ph.karno  physician's estimate of the Karnofsky score, a competitor to the
#              ECOG ps.
# pat.karno  patient's assesment of his/her Karnofsky score
# meal.cal   # calories consumed at meals (exclude beverages and snacks)
# wt.loss    weight loss in the last 6 months


# 12/98: this is now part of the Splus distribution, but they call it
#  "lung"
cancer <- lung
