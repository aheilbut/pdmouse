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

mo430names = pandas.read_table("/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/Oct29/mo4302genenames.tab")
mo430names.index =  mo430names.probe_id

mo430info = mo430symbol.merge(mo430names)
mo430info.index = mo430info.probe_id

# define subsets

def define_subsets():
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

    ss_cp101_allchronic = pd_covar.select(lambda x: (pd_covar.ix[x, "MouseType"] == "CP101" 
                                                    and pd_covar.ix[x, "LesionType"] == "6-OHDA" 
                                                    and (pd_covar.ix[x, "DrugTreat"] == "Chronic low levodopa" 
                                                         or pd_covar.ix[x, "DrugTreat"] == "Chronic high levodopa")
                                                    and pd_covar.ix[x, "MouseID"] != 1343)) 

    ss_cp101_allchronic = pd_covar.select(lambda x: (pd_covar.ix[x, "MouseType"] == "CP101" and pd_covar.ix[x, "LesionType"] == "6-OHDA" and (pd_covar.ix[x, "DrugTreat"] == "Chronic low levodopa" or pd_covar.ix[x, "DrugTreat"] == "Chronic high levodopa") and pd_covar.ix[x, "MouseID"] != 1343)) 


    ss_cp73_allchronic = pd_covar.select(lambda x: (pd_covar.ix[x, "MouseType"] == "CP73" 
                                                    and pd_covar.ix[x, "LesionType"] == "6-OHDA" 
                                                    and (pd_covar.ix[x, "DrugTreat"] == "Chronic low levodopa" \
                                                         or pd_covar.ix[x, "DrugTreat"] == "Chronic high levodopa")
                                                    and pd_covar.ix[x, "MouseID"] != 1343)) 


    ss_cp73_allChronic_LDOPA = pd_covar.select(lambda x: (pd_covar.ix[x, "MouseType"] == "CP73" 
                                                    and pd_covar.ix[x, "LesionType"] == "6-OHDA" 
                                                    and (pd_covar.ix[x, "DrugTreat"] == "Chronic high levodopa" or pd_covar.ix[x, "DrugTreat"] == "Chronic low levodopa")))
                                              
    ss_cp101_allChronic_LDOPA = pd_covar.select(lambda x: (pd_covar.ix[x, "MouseType"] == "CP101" 
                                                    and pd_covar.ix[x, "LesionType"] == "6-OHDA" 
                                                    and pd_covar.ix[x, "MouseID"] != 1343
                                                    and (pd_covar.ix[x, "DrugTreat"] == "Chronic high levodopa" or pd_covar.ix[x, "DrugTreat"] == "Chronic low levodopa")))
                                          

    cp73_chronic_saline = pd_covar.select(lambda x: pd_covar.ix[x, "MouseType"] == "CP73" 
                                                    and pd_covar.ix[x, "LesionType"] == "6-OHDA" 
                                                    and pd_covar.ix[x, "DrugTreat"] == "Chronic saline")
    
    cp101_chronic_saline = pd_covar.select(lambda x: pd_covar.ix[x, "MouseType"] == "CP101" 
                                                    and pd_covar.ix[x, "LesionType"] == "6-OHDA"
                                                    and pd_covar.ix[x, "DrugTreat"] == "Chronic saline")


    return locals()




aim_model_set = ['aim_models_cp73_chronicHigh', 'aim_models_cp73_chronicLow', 
                 'aim_models_cp101_chronicHigh', 'aim_models_cp101_chronicLow']

print "save pickle models"
#for m in aim_model_set:
#    cPickle.dump(eval(m), open("/data/adrian/data/temp/%s.pickle" % m,"w"), protocol=-1)

aim_models_cp101_allchronic = pda.fitModels(pd_all.index, ss_cp101_allchronic, pd_all)
aim_models_cp73_allchronic = pda.fitModels(pd_all.index, ss_cp73_allchronic, pd_all)

aim_models_cp73_chronic = pda.fitModels(pd_all.index, ss_cp73_allChronic_LDOPA)
aim_models_cp101_chronic = pda.fitModels(pd_all.index, ss_cp101_allChronic_LDOPA)


print "load pickled models"
for m in aim_model_set:
    locals()[m] = cPickle.load(open("/data/adrian/data/temp/%s.pickle" % m))

print "calc change stats"
cp73_highVsSaline_changeStats = pda.calcChangeStats(pd_all, ss_cp73_chronicHigh, cp73_chronic_saline)
cp101_highVsSaline_changeStats = pda.calcChangeStats(pd_all, ss_cp101_chronicHigh, cp101_chronic_saline) 
cp73_lowVsSaline_changeStats = pda.calcChangeStats(pd_all, ss_cp73_chronicLow, cp73_chronic_saline) 
cp101_lowVsSaline_changeStats = pda.calcChangeStats(pd_all, ss_cp101_chronicLow, cp101_chronic_saline) 

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


cp73_highVsLow_changeStats = pda.calcChangeStats(pd_all, ss_cp73_chronicHigh, ss_cp73_chronicLow)
cp101_highVsLow_changeStats = pda.calcChangeStats(pd_all, ss_cp101_chronicHigh, ss_cp101_chronicLow)

cp73_highVsLow_ttest = pda.calc_ttest(pd_all, ss_cp73_chronicHigh, ss_cp73_chronicLow)
cp101_highVsLow_ttest = pda.calc_ttest(pd_all, ss_cp101_chronicHigh, ss_cp101_chronicLow)

cp73_highVsLow_changeStats = cp73_highVsLow_changeStats.merge(cp73_highVsLow_ttest, left_index=True, right_index=True)
cp101_highVsLow_changeStats = cp101_highVsLow_changeStats.merge(cp101_highVsLow_ttest, left_index=True, right_index=True)


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


cp73_wide = cp73_lowVsSaline_changeStats.merge( cp73_HighVsSaline_changeStats, left_index=True, right_index=True )


cp73_p = pandas.Panel(data = { "highVsSaline" : cp73_highVsSaline_changeStats, 
                               "lowVsSaline" : cp73_lowVsSaline_changeStats, 
                               "highVsLow" : cp73_highVsLow_changeStats,
                               "bothVsAIM" : cp73_High_AIMcorr } )

def selector(x):
    return ( cp73_p.ix["highVsSaline"]["fc_exp"][x] > 1
        and cp73_p.ix["highVsLow"]["fc_exp"][x] > 1
    )

cp73_p.select( axis=1, crit = selector)

#makePlots(cp101_Low_2xFC_AIMcorr, "/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/cp101_Low_2xFC_AIMcorr.pdf")
#makePlots(cp101_High_2xFC_AIMcorr, "/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/cp101_High_2xFC_AIMcorr.pdf")
#makePlots(cp73_Low_2xFC_AIMcorr, "/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/cp73_Low_2xFC_AIMcorr.pdf")
#makePlots(cp73_High_2xFC_AIMcorr, "/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/cp73_High_2xFC_AIMcorr.pdf")

#makePlots(cp101_Low_AIMcorr, "/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/cp101_Low_AIMcorr.pdf")
#makePlots(cp101_High_AIMcorr, "/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/cp101_High_AIMcorr.pdf")
#makePlots(cp73_Low_AIMcorr, "/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/cp73_Low_AIMcorr.pdf")
makePlots(pd_all, pd_covar, cp73_High_AIMcorr, "/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/cp73_High_AIMcorr.pdf")

pda.makePlots(pd_all, pd_covar, cp101_High_2xFC_AIMcorr, "/data/adrian/test.pdf")


cp73_highVsSaline_changeStats[["fc_exp", "cv_exp", "pval", "bh"]]
    .merge( cp73_lowVsSaline_changeStats["fc_exp", "cv_exp", "pval", "bh"], 
            left_index=True, right_index=True, suffixes=("_cp73_hi", "_cp73_lo") )

selectdim = ["fc_exp", "cv_exp", "bh"]

o_cp73_chronicHigh_aim = aim_models_cp73_chronicHigh[["rsq_adj","pval"]]
o_cp73_chronicHigh_aim.columns = ["aim_rsq_cp73_hi", "aim_pval_cp73_hi"]

o_cp73_chronicLow_aim = aim_models_cp73_chronicLow[["rsq_adj","pval"]]
o_cp73_chronicLow_aim.columns = ["aim_rsq_cp73_lo", "aim_pval_cp73_lo"]

o_cp101_chronicHigh_aim = aim_models_cp73_chronicHigh[["rsq_adj","pval"]]
o_cp101_chronicHigh_aim.columns = ["aim_rsq_cp101_hi", "aim_pval_cp101_hi"]

o_cp101_chronicLow_aim = aim_models_cp73_chronicLow[["rsq_adj","pval"]]
o_cp101_chronicLow_aim.columns = ["aim_rsq_cp101_lo", "aim_pval_cp101_lo"]


o_cp73 = cp73_highVsSaline_changeStats.ix[:,selectdim].merge( cp73_lowVsSaline_changeStats.ix[:,selectdim], left_index=True, right_index=True, suffixes=("_cp73_t_hi", "_cp73_t_lo") )
 
o_cp101 = cp101_highVsSaline_changeStats.ix[:,selectdim].merge( cp101_lowVsSaline_changeStats.ix[:,selectdim], left_index=True, right_index=True, suffixes=("_cp101_t_hi", "_cp101_t_lo") ) 

o_m_cp73 = o_cp73.merge( o_cp73_chronicHigh_aim, left_index=True, right_index=True ).merge(o_cp73_chronicLow_aim, left_index=True, right_index=True )
o_m_cp101 = o_cp101.merge( o_cp101_chronicHigh_aim, left_index=True, right_index=True ).merge(o_cp101_chronicLow_aim, left_index=True, right_index=True )

o = o_m_cp73.merge(o_m_cp101, left_index=True, right_index=True)


m = []
for k in alib.wikipath.hs_wp_symbols.keys():
    matches = [a in set(alib.wikipath.hs_wp_symbols[k]) for a in [str(s).upper() for s in ss_cp73_filtered_allDoseAIMcorr.symbol]]
    m.append( (k, np.sum(matches), len(alib.wikipath.hs_wp_symbols[k]), len(matches)) )
 
"""

   
fits_expr_dose = cPickle.load(open("/data/adrian/data/temp/fits_exprDoseModel.pickle", "rb"))
fits_expr = cPickle.load(open("/data/adrian/data/temp/fits_expr.pickle", "rb"))
    

    model = "AIM ~ dose * expression"
    model = "AIM ~ dose + expression"
    model = "AIM ~ dose"
    model = "AIM ~ expression"


@dview.parallel
    def gomodel(x):
        return pda.fitAIM(x, ss_cp101_allchronic, pd_all)
       
start = time.time()
%px y = pda.fitModels(a, ss_cp101_allchronic, pd_all)
end = time.time()
print end-start     

start = time.time()
z = gomodel(list(pd_all.index))
end = time.time()
print end-start     


start = time.time()
y = pda.fitModels(a, ss_cp101_allchronic, pd_all)
end = time.time()
print end-start
"""