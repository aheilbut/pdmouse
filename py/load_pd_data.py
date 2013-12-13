# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import pandas
import pd_locals

# <codecell>

# load data 
pd_all = pandas.read_table(pd_locals.datadir + "/2012_10_29/PD_arraydata.tab")
pd_covar = pandas.read_table(pd_locals.datadir + "/2012_10_29/pd.covar.tab")

mo430symbol = pandas.read_table(pd_locals.datadir + "/2013_july_9/mo4302symbols.tab")
mo430symbol.index = mo430symbol.probe_id

mo430names = pandas.read_table(pd_locals.datadir + "/2013_july_9/mo4302genenames.tab")
mo430names.index =  mo430names.probe_id

mo430info = mo430symbol.merge(mo430names)
mo430info.index = mo430info.probe_id

# <codecell>

f

