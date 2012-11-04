# -*- coding: utf-8 -*-
"""
Created on Sat Nov  3 22:58:57 2012

@author: adrian
"""

import pd_analysis as pda

import pandas
import cPickle

# load data 
pd_all = pandas.read_table("/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/Oct29/PD_arraydata.tab")
pd_covar = pandas.read_table("/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/Oct29/pd.covar.tab")

mo430symbol = pandas.read_table("/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/Oct29/mo4302symbols.tab")
mo430symbol.index = mo430symbol.probe_id

mo430names =pandas.read_table("/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/Oct29/mo4302genenames.tab")
mo430names.index =  mo430names.probe_id

mo430info = mo430symbol.merge(mo430names)
mo430info.index = mo430info.probe_id

# define subsets
ss_cp73_chronicHigh = pd_covar.select(lambda x: (pd_covar.ix[x, "MouseType"] == "CP73" 
                                                and pd_covar.ix[x, "LesionType"] == "6-OHDA" 
                                                and pd_covar.ix[x, "DrugTreat"] == "Chronic high levodopa"))

ss_cp73_chronicLow = pd_covar.select(lambda x: (pd_covar.ix[x, "MouseType"] == "CP73" 
                                                and pd_covar.ix[x, "LesionType"] == "6-OHDA" 
                                                and pd_covar.ix[x, "DrugTreat"] == "Chronic low levodopa"))

ss_cp101_chronicHigh = pd_covar.select(lambda x: (pd_covar.ix[x, "MouseType"] == "CP101" 
                                                and pd_covar.ix[x, "LesionType"] == "6-OHDA" 
                                                and pd_covar.ix[x, "DrugTreat"] == "Chronic high levodopa")) 

ss_cp101_chronicLow = pd_covar.select(lambda x: (pd_covar.ix[x, "MouseType"] == "CP101" 
                                                and pd_covar.ix[x, "LesionType"] == "6-OHDA" 
                                                and pd_covar.ix[x, "DrugTreat"] == "Chronic low levodopa" 
                                                and pd_covar.ix[x, "MouseID"] != 1343)) 


cp73_chronic_saline = pd_covar.select(lambda x: pd_covar.ix[x, "MouseType"] == "CP73" and pd_covar.ix[x, "LesionType"] == "6-OHDA" and pd_covar.ix[x, "DrugTreat"] == "Chronic saline")
cp101_chronic_saline = pd_covar.select(lambda x: pd_covar.ix[x, "MouseType"] == "CP101" and pd_covar.ix[x, "LesionType"] == "6-OHDA" and pd_covar.ix[x, "DrugTreat"] == "Chronic saline")

aim_model_set = ['aim_models_cp73_chronicHigh', 'aim_models_cp73_chronicLow', 
                 'aim_models_cp101_chronicHigh', 'aim_models_cp101_chronicLow']

print "save pickle models"
#for m in aim_model_set:
#    cPickle.dump(eval(m), open("/data/adrian/data/temp/%s.pickle" % m,"w"), protocol=-1)

print "load pickled models"
for m in aim_model_set:
    locals()[m] = cPickle.load(open("/data/adrian/data/temp/%s.pickle" % m))

print "calc change stats"
cp73_highVsSaline_changeStats = pda.calcChangeStats(pd_all, ss_cp73_chronicHigh, cp73_chronic_saline)
cp101_highVsSaline_changeStats = pda.calcChangeStats(pd_all, ss_cp101_chronicHigh, cp101_chronic_saline) 
cp73_lowVsSaline_changeStats = pda.calcChangeStats(pd_all, ss_cp73_chronicLow, cp73_chronic_saline) 
cp101_lowVsSaline_changeStats = pda.calcChangeStats(pd_all, ss_cp101_chronicLow, cp73_chronic_saline) 

print "calc t tests"
cp73_high_ttest = pda.calc_ttest(pd_all, ss_cp73_chronicHigh, cp73_chronic_saline)
cp73_low_ttest = pda.calc_ttest(pd_all, ss_cp73_chronicLow, cp73_chronic_saline)
cp101_high_ttest = pda.calc_ttest(pd_all, ss_cp101_chronicHigh, cp101_chronic_saline)
cp101_low_ttest = pda.calc_ttest(pd_all, ss_cp101_chronicLow, cp101_chronic_saline)

print "merge results"
cp73_highVsSaline_changeStats = cp73_highVsSaline_changeStats.merge(cp73_high_ttest, left_index=True, right_index=True)
cp73_lowVsSaline_changeStats = cp73_lowVsSaline_changeStats.merge(cp73_low_ttest, left_index=True, right_index=True)
cp101_highVsSaline_changeStats = cp101_highVsSaline_changeStats.merge(cp101_high_ttest, left_index=True, right_index=True)
cp101_lowVsSaline_changeStats = cp101_lowVsSaline_changeStats.merge(cp101_low_ttest, left_index=True, right_index=True)


cp101_High_2xFC_AIMcorr = pda.changeFilter(aim_models_cp101_chronicHigh,
                                       cp101_highVsSaline_changeStats, mo430info)

cp101_Low_2xFC_AIMcorr = pda.changeFilter(aim_models_cp101_chronicLow, 
                                      cp101_lowVsSaline_changeStats, mo430info)

cp73_High_2xFC_AIMcorr = pda.changeFilter(aim_models_cp73_chronicHigh, 
                                      cp73_highVsSaline_changeStats, mo430info)
cp73_Low_2xFC_AIMcorr = pda.changeFilter(aim_models_cp73_chronicLow, 
                                     cp73_lowVsSaline_changeStats, mo430info)

cp101_High_AIMcorr = pda.changeFilter(aim_models_cp101_chronicHigh, cp101_highVsSaline_changeStats, mo430info,
                                  aim_fit_max_pval=0.01, 
                                  fc_min=0,
                                  diff_max_bh_pval=1)
cp101_Low_AIMcorr = pda.changeFilter(aim_models_cp101_chronicLow, cp101_lowVsSaline_changeStats, mo430info, aim_fit_max_pval=0.01, fc_min=0, diff_max_bh_pval=1)

cp73_High_AIMcorr = pda.changeFilter(aim_models_cp73_chronicHigh, cp73_highVsSaline_changeStats, mo430info, aim_fit_max_pval=0.01, fc_min=0, diff_max_bh_pval=1)
cp73_Low_AIMcorr = pda.changeFilter(aim_models_cp73_chronicLow, cp73_lowVsSaline_changeStats, mo430info, aim_fit_max_pval=0.01, fc_min=0, diff_max_bh_pval=1)



#makePlots(cp101_Low_2xFC_AIMcorr, "/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/cp101_Low_2xFC_AIMcorr.pdf")
#makePlots(cp101_High_2xFC_AIMcorr, "/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/cp101_High_2xFC_AIMcorr.pdf")
#makePlots(cp73_Low_2xFC_AIMcorr, "/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/cp73_Low_2xFC_AIMcorr.pdf")
#makePlots(cp73_High_2xFC_AIMcorr, "/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/cp73_High_2xFC_AIMcorr.pdf")

#makePlots(cp101_Low_AIMcorr, "/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/cp101_Low_AIMcorr.pdf")
#makePlots(cp101_High_AIMcorr, "/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/cp101_High_AIMcorr.pdf")
#makePlots(cp73_Low_AIMcorr, "/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/cp73_Low_AIMcorr.pdf")
makePlots(pd_all, pd_covar, cp73_High_AIMcorr, "/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/cp73_High_AIMcorr.pdf")

pda.makePlots(pd_all, pd_covar, cp101_High_2xFC_AIMcorr, "/data/adrian/test.pdf")
