#import rpy2
import psycopg2
import numpy as np
import networkx as nx
import pandas
import matplotlib as mp
import matplotlib.pyplot as plt
import statsmodels.api as sm
import random
from statsmodels.formula.api import ols 
from matplotlib.backends.backend_pdf import PdfPages
import cPickle
import datetime
import time, datetime
import scipy.stats as st
import statsmodels.sandbox.stats.multicomp
import pd_locals
import json

mo430symbol = pandas.read_table(pd_locals.datadir +  "/2012_10_29/mo4302symbols.tab")
mo430symbol.index = mo430symbol.probe_id

mo430names = pandas.read_table(pd_locals.datadir + "/2012_10_29/mo4302genenames.tab")
mo430names.index =  mo430names.probe_id

mo430info = mo430symbol.merge(mo430names)
mo430info.index = mo430info.probe_id

pd_all = pandas.read_table(pd_locals.datadir + "/2012_10_29/PD_arraydata.tab")
pd_covar = pandas.read_table(pd_locals.datadir + "/2012_10_29/pd.covar.tab")


class TagMan():
    def __init__(self):
        pass
    
    # encode
    def e(self, tagGroup):
        # sort tag tree recursively 
        return json.dumps( sorted( tagGroup ) )
        
    # decode
    def d(self, tagString):
        return json.loads(tagString)
    
    # match
    # return all the strings that match the specified tags
    # if only a tag type is specified, then any string containing that tag is ok, but must have that tag
    # if a value is also specified, then value must match too
    def m(self, searchTagList, encodedStrings):
        resultStrings = []
        for s in encodedStrings:
            try:
                td = dict(self.d(s))
            except:
                print "not an encoded tag list"
            current_match = True
            for searchTag in searchTagList:
                if isinstance(searchTag, (list, tuple)):
                    if td.get(searchTag[0], None) == searchTag[1]:
                        pass
                    else:
                        current_match = False
                        
                else:
                    if td.has_key(searchTag):
                        pass
                    else:
                        current_match = False
                    
                    
            if current_match is True:
                resultStrings.append(s)
                
        return resultStrings
    

tm = TagMan()
    
    # filter

    
#def tagDF(df, df_col_type, tags):
#    """ return the given dataframe with columns relabeled to have the tags specified, 
#    and existing column names changed to tags of type df_col_type """
#    df.columns = [ tags + ((df_col_type, cn),) for cn in df.columns]
#    return df
    
    
def uniqueSym(symbolSeries):
    symbolSeries = symbolSeries.fillna("-")
    r = set(symbolSeries)
    r.discard(np.nan)
    r.discard("-")
    return(r)
    

def fitAIM(probeset, covar_subset, dataset, AIM_dimension="totalAIM"):
    model = "AIM ~ expression"
    cur_gene = pandas.DataFrame( { "expression" : dataset.ix[probeset,covar_subset.filenames], 
                                  "AIM" : list(covar_subset[AIM_dimension].fillna(0)) } )
    test_model = ols(model, cur_gene).fit()
    return test_model

def fit_dose_expr_model(probeset, covar_subset, dataset, AIM_dimension="totalAIM"):
    model = "AIM ~ expression*dose"
    cur_gene = pandas.DataFrame( { "expression" : dataset.ix[probeset, covar_subset.filenames],
                                          "AIM" : list(covar_subset[AIM_dimension].fillna(0)),
                                         "dose" : [int(i) for i in covar_subset["DrugTreat"] == "Chronic high levodopa"]
                                })
    test_model = ols(model, cur_gene).fit()
    return test_model

def fit_model(model, probeset, covar_subset, dataset, AIM_dimension="totalAIM"):
    cur_gene = pandas.DataFrame( { "expression" : dataset.ix[probeset, covar_subset.filenames],
                                          "AIM" : list(covar_subset[AIM_dimension].fillna(0)),
                                         "dose" : [int(i) for i in covar_subset["DrugTreat"] == "Chronic high levodopa"]
                                })
    test_model = ols(model, cur_gene).fit()
    return test_model



def fit_model_to_all_probes(model, probeset_list, covar_subset, dataset):
    fits = pandas.DataFrame(index = probeset, 
                            data = [fit_model(model, p, covar_subset, dataset) for p in probeset], 
                            columns=["model"])
    fits["pval"] = fits["model"].map( lambda x: x.pvalues["expression"] )
    fits["expr_coef"] = fits["model"].map( lambda x: x.params["expression"] )
    fits["bh"] = statsmodels.sandbox.stats.multicomp.multipletests(fits.pval, method="fdr_bh")[1]
    fits["bonf"] = statsmodels.sandbox.stats.multicomp.multipletests(fits.pval, method="bonferroni")[1]
    fits["sign"] = np.sign(fits["expr_coef"])
    fits["rsq_adj"] = fits["model"].map( lambda x: x.rsquared_adj )
    
    return fits                                  

def fitModels(probeset, covar_subset, dataset):
    fits = pandas.DataFrame( index = probeset, 
                            data = [fitAIM(p, covar_subset, dataset) for p in probeset], 
                            columns=["model"])
    fits["pval"] = fits["model"].map( lambda x: x.pvalues["expression"] )
    fits["expr_coef"] = fits["model"].map( lambda x: x.params["expression"] )
    fits["bh"] = statsmodels.sandbox.stats.multicomp.multipletests(fits.pval, method="fdr_bh")[1]
    fits["bonf"] = statsmodels.sandbox.stats.multicomp.multipletests(fits.pval, method="bonferroni")[1]
    fits["sign"] = np.sign(fits["expr_coef"])
    fits["rsq_adj"] = fits["model"].map( lambda x: x.rsquared_adj )
    
    return fits

def dofitModels():
    aim_models_cp73_chronicHigh = fitModels(pd_all.index, ss_cp73_chronicHigh)
    aim_models_cp73_chronicLow = fitModels(pd_all.index, ss_cp73_chronicLow)
    
    aim_models_cp101_chronicHigh = fitModels(pd_all.index, ss_cp101_chronicHigh) 
    aim_models_cp101_chronicLow = fitModels(pd_all.index, ss_cp101_chronicLow)

    return (aim_models_cp73_chronicHigh, aim_models_cp73_chronicLow, aim_models_cp101_chronicHigh, aim_models_cp101_chronicLow)


def plotProbe(d, pd_covar, probeset, mousetype, ax, legend=False, xlim=None):   
    try:
        symbol = mo430info.ix[probeset,"symbol"]
        genename = mo430info.ix[probeset,"gene_name"]
    except:
        symbol = "--"
        genename = "--"

    #ax.patch.set_facecolor((0.9, 0.9, 0.9))
    ax.grid(which='both')        
        
    ax.set_title(mousetype + " : " + probeset + " : " + symbol + "\n" + genename)
    ax.set_xlabel("expression")
    ax.set_ylabel("sum of AIM")
    ax.set_ylim((-5, 160))
    
    ss = pd_covar.select(lambda x: pd_covar.ix[x, "MouseType"] == mousetype 
                         and pd_covar.ix[x, "LesionType"] == "Ascorbate" 
                         and pd_covar.ix[x, "DrugTreat"] == "Chronic saline")  
    cur_gene = pandas.DataFrame( { "expression" : d.ix[probeset,ss.filenames], "AIM" : list(ss.totalAIM) } )
    ax.plot(cur_gene["expression"], cur_gene["AIM"], "y.", ms=10, label="Ascorbate / Saline", alpha=0.6 )    
    
    ss = pd_covar.select(lambda x: pd_covar.ix[x, "MouseType"] == mousetype 
                         and pd_covar.ix[x, "LesionType"] == "6-OHDA" 
                         and pd_covar.ix[x, "DrugTreat"] == "Chronic high levodopa")  
    cur_gene = pandas.DataFrame( { "expression" : d.ix[probeset,ss.filenames], "AIM" : list(ss.totalAIM) } )
    ax.plot(cur_gene["expression"], cur_gene["AIM"], "b.", ms=10, label="Chronic High L-DOPA", alpha=0.6 )
    
    ss = pd_covar.select(lambda x: pd_covar.ix[x, "MouseType"] == mousetype 
                         and pd_covar.ix[x, "LesionType"] == "6-OHDA" 
                         and pd_covar.ix[x, "DrugTreat"] == "Chronic low levodopa")  
    cur_gene = pandas.DataFrame( { "expression" : d.ix[probeset,ss.filenames], "AIM" : list(ss.totalAIM) } )
    ax.plot(cur_gene["expression"], cur_gene["AIM"], "r.", ms=10, label="Chronic low L-DOPA", alpha=0.6 )    
    
    ss = pd_covar.select(lambda x: pd_covar.ix[x, "MouseType"] == mousetype and pd_covar.ix[x, "LesionType"] == "6-OHDA" and pd_covar.ix[x, "DrugTreat"] == "Acute high levodopa")  
    cur_gene = pandas.DataFrame( { "expression" : d.ix[probeset,ss.filenames], "AIM" : list(ss.totalAIM) } )
    ax.plot(cur_gene["expression"], cur_gene["AIM"], "m.", ms=10, label="Acute L-DOPA", alpha=0.6)    

    ss = pd_covar.select(lambda x: pd_covar.ix[x, "MouseType"] == mousetype 
                         and pd_covar.ix[x, "LesionType"] == "6-OHDA" 
                         and pd_covar.ix[x, "DrugTreat"] == "Chronic saline")  
    cur_gene = pandas.DataFrame( { "expression" : d.ix[probeset,ss.filenames], "AIM" : list(ss.totalAIM) } )
    ax.plot(cur_gene["expression"], cur_gene["AIM"], "g.", ms=10, label="6-OHDA / Saline", alpha=0.6)        
    
    if xlim is not None:
        ax.set_xlim(xlim)
    
    if (legend == True):
        #ax.legend()
        ax.legend( bbox_to_anchor=(1.2, 0), loc='lower left', borderaxespad=0., title="groups")
           
    return True # plt.gca().get_legend_handles_labels()



def probe_boxPlots(probeset, d, groups, mousetype, ax):
    mp.rcParams.update({'font.size': 10 })
    l = []
    for (group_name, group_colnames) in groups:
        l.append( d.ix[probeset, group_colnames] )
        
    p = ax.boxplot(l)
    for (i, d) in enumerate(l):
        z = [(i+1)] * len(d)
        ax.plot(z + np.random.rand(len(d))*0.2 - 0.1, d, '.', ms=5)
    ax.xaxis.set_ticks(np.arange(0.5, len(l) + 1))
    ax.xaxis.set_ticklabels( [g[0] for g in groups], rotation=45 )
                      

    try:
        symbol = mo430info.ix[probeset,"symbol"]
        genename = mo430info.ix[probeset,"gene_name"]
    except:
        symbol = "--"
        genename = "--"

    #ax.patch.set_facecolor((0.9, 0.9, 0.9))
        
    ax.set_title(mousetype + " : " + probeset + " : " + symbol + "\n" + genename + "\n expression vs. lesion / drug treatment\n")
    ax.set_ylabel("log2 expression")
    ax.set_xlabel("treatment")
    ax.grid(which='both')        

    return True
    
def plotBoth(d, covar, probeset):
    f = plt.figure(figsize=(18,5))
    plt.rcParams.update({'font.size': 8})
#    subplot2grid((1,5), (0,0), colspan=2)
    plt.subplot(1,3,1)
    plotProbe(d, covar, probeset, "CP73")
#    subplot2grid((1,5), (0,2), colspan=2)
    plt.subplot(1,3,2)
    (handles, labels) = plotProbe(d, covar, probeset, "CP101")
#    subplot2grid((1,5), (0,4), colspan=1)
    ax_legend = plt.subplot(1,3,3)
    ax_legend.legend(handles, labels, numpoints=1, loc='upper left', fancybox=True, shadow=True)
    ax_legend.set_frame_on(False)
    ax_legend.set_xticks([])
    ax_legend.set_yticks([])
    return f
 
    
def calcChangeStats(d, exp_set, control_set, tags=()):
    return pandas.DataFrame(
        { tm.e( tags + (("st", "median"),("cg", "exp")) )      : d[list(exp_set.filenames)].median(axis=1),
          tm.e( tags + (("st", "mean"), ("cg", "exp"))  )      : d[list(exp_set.filenames)].mean(axis=1),
          tm.e( tags + (("st", "fc_medians"),)  )              : d[list(exp_set.filenames)].median(axis=1) - d[list(control_set.filenames)].median(axis=1),
          tm.e( tags + (("st", "fc_means"),)    )              : d[list(exp_set.filenames)].mean(axis=1) - d[list(control_set.filenames)].mean(axis=1),
          tm.e( tags + (("st", "stddev"), ("cg", "exp")) )     : d[list(exp_set.filenames)].std(axis=1),
          tm.e( tags + (("st", "stddev"), ("cg", "ctrl")) )    : d[list(exp_set.filenames)].std(axis=1),
          tm.e( tags + (("st", "cv"), ("cg", "exp")) )         : d[list(exp_set.filenames)].std(axis=1) / d[list(exp_set.filenames)].mean(axis=1)
        })


def calc_ttest(data, exp_set, control_set, tags=()):
    d = [ st.ttest_ind( data.ix[probeset, list(exp_set.filenames)], 
                             data.ix[probeset, list(control_set.filenames)], equal_var=False) for probeset in data.index]
    rs = pandas.DataFrame( index=data.index, data=d, columns=[ tm.e( tags+(("st", "t"),("tt", "welch ttest"))), tm.e( tags + (("st", "pval"), ("tt", "welch ttest"), ("mc", "nominal") ))])    
    rs[tm.e( tags + (("tt", "welch ttest"), ("st", "pval"), ("mc", "bonf")))] = statsmodels.sandbox.stats.multicomp.multipletests(rs.ix[:, tm.e( tags + (("st", "pval"), ("tt", "welch ttest"), ("mc", "nominal"))) ], method="bonferroni")[1]
    rs[tm.e( tags + (("tt", "welch ttest"), ("st", "pval"), ("mc", "bh")))] = statsmodels.sandbox.stats.multicomp.multipletests(rs.ix[:, tm.e( tags + (("st", "pval"), ("tt", "welch ttest"), ("mc", "nominal")))], method="fdr_bh")[1] 

    d = [ st.ttest_ind( data.ix[probeset, list(exp_set.filenames)], 
                             data.ix[probeset, list(control_set.filenames)], equal_var=True) for probeset in data.index]

    rs[tm.e( tags+(("st", "t"),("tt", "student ttest")))] = [v[0] for v in d]
    rs[tm.e( tags + (("st", "pval"), ("tt", "student ttest"), ("mc", "nominal") ))] = [v[1] for v in d]
    
    rs[tm.e( tags + (("st", "pval"), ("tt", "student ttest"), ("mc", "bonf")))] = statsmodels.sandbox.stats.multicomp.multipletests(rs.ix[:, tm.e( tags + (("st", "pval"), ("tt", "student ttest"), ("mc", "nominal"))) ], method="bonferroni")[1]
    rs[tm.e( tags + (("st", "pval"), ("tt", "student ttest"), ("mc", "bh")))] = statsmodels.sandbox.stats.multicomp.multipletests(rs.ix[:, tm.e( tags + (("st", "pval"), ("tt", "student ttest"), ("mc", "nominal")))], method="fdr_bh")[1] 


    # do diagnostic tests for heteroskedasticity
    d = [st.levene( data.ix[probeset, list(exp_set.filenames)], data.ix[probeset, list(control_set.filenames)]) for probeset in data.index ]
    rs[ tm.e( tags + (("tt", "levene"), ("st", "pval")))] = [z[1] for z in d]

    # omnibus test for normality
#    d = [st.normaltest( data.ix[probeset, list(exp_set.filenames)]) for probeset in data.index ]
#    rs[ tm.e( tags + (("tt", "d-p omnibus"), ("st", "pval"), ("cg", "exp") ))] = [z[1] for z in d]

#    d = [st.normaltest( data.ix[probeset, list(control_set.filenames)]) for probeset in data.index ]
#    rs[ tm.e( tags + (("tt", "d-p omnibus"), ("st", "pval"), ("cg", "ctrl") ))] = [z[1] for z in d]

    return rs

def changeFilter(model_stats, change_stats, probeinfo, fc_min=1.0, diff_max_bh_pval=0.20, aim_fit_max_pval=0.05, cv_min=0.01):
    return (model_stats
    .merge(change_stats, suffixes=("_vs_AIM", "_vsSaline"), left_index=True, right_index=True)
    .select(lambda x: change_stats.ix[x, "bh"] < diff_max_bh_pval and abs(change_stats.ix[x, "fc_exp"]) > fc_min  )
    .select(lambda x: model_stats.ix[x, "pval"] < aim_fit_max_pval )
    .select(lambda x: change_stats.ix[x, "cv_exp"] > cv_min )
    .merge(probeinfo, how='left', left_index=True, right_index=True)
    .sort(["fc_exp", "sign", "pval_vs_AIM"], ascending=True))


def calc_save_datasets():

    xlr = pandas.ExcelWriter("/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/PD_L-DOPA_AIM_corr_Nov2.xls")
    cp101_Low_2xFC_AIMcorr.drop("model", axis=1).to_excel(xlr, "CP101 Low 2xFC AIMcorr")
    cp101_High_2xFC_AIMcorr.drop("model", axis=1).to_excel(xlr, "CP101 High 2xFC AIMcorr")
    cp73_Low_2xFC_AIMcorr.drop("model", axis=1).to_excel(xlr, "CP73 Low 2xFC AIMcorr")
    cp73_High_2xFC_AIMcorr.drop("model", axis=1).to_excel(xlr, "CP73 High 2xFC AIMcorr")
    
    cp101_Low_AIMcorr.drop("model", axis=1).to_excel(xlr, "CP101 Low AIMcorr")
    cp101_High_AIMcorr.drop("model", axis=1).to_excel(xlr, "CP101 High AIMcorr")
    cp73_Low_AIMcorr.drop("model", axis=1).to_excel(xlr, "CP73 Low AIMcorr")
    cp73_High_AIMcorr.drop("model", axis=1).to_excel(xlr, "CP73 High AIMcorr")
    
    xlr.save()
 
def calc_intersections():   
    xlr = pandas.ExcelWriter("/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/PD_L-DOPA_AIM_corr_intersections_Nov3.xls")

    (cp101_High_AIMcorr
        .drop("model", axis=1)
        .merge(cp101_Low_AIMcorr.drop("model", axis=1), suffixes=("cp101hi", "cp101lo"), left_index=True, right_index=True)
        .to_excel(xlr, "cp101 hi + lo"))
    
    (cp73_High_AIMcorr.drop("model", axis=1)
        .merge(cp73_Low_AIMcorr
        .drop("model", axis=1), suffixes=("cp73hi","cp73lo"), left_index=True, right_index=True)
        .to_excel(xlr, "cp73 hi + lo"))
    
    (cp73_High_AIMcorr.drop("model", axis=1)
        .merge(cp101_High_AIMcorr
        .drop("model", axis=1), suffixes=("cp73hi", "cp101hi"), left_index=True, right_index=True)
        .to_excel(xlr, "cp101 hi + cp73 hi"))
    
    (cp73_Low_AIMcorr.drop("model", axis=1)
        .merge(cp101_Low_AIMcorr
        .drop("model", axis=1), suffixes=("cp73hi", "cp101hi"), left_index=True, right_index=True)
        .to_excel(xlr, "cp101 lo + cp73 lo"))
    
    xlr.save()



def makePlots(dataset, covar, resultSet, filename):
    pp = PdfPages(filename)
    for p in resultSet.index:
        f = plotBoth(dataset, covar, p)
        plt.text(0, 0, resultSet.ix[p, :].to_string(), size=8)
        plt.text(0.6, 0.9, time.strftime("%a, %d %b %Y \n %H:%M:%S", time.localtime()))
        plt.savefig(pp, format="pdf")
        plt.clf()
    pp.close()

def plotAIMgraphs():    
    plt.figsize(40, 10)
    plt.rcParams.update({'font.size': 18})
    fig = plt.figure()
    patches_list = []
    label_list = []
    for (mousetype_i, mousetype) in enumerate(["CP101", "CP73"]):
        s = pd_covar.select(lambda x: pd_covar.ix[x, "MouseType"] == mousetype and pd_covar.ix[x, "LesionType"] == "6-OHDA") 
        ax1 = fig.add_subplot(1, 3, mousetype_i + 1)
        ax1.patch.set_facecolor((0.9, 0.9, 0.9))
        ax1.grid(which='both')
        plt.title(mousetype)
        plt.ylabel("AIM score")
        plt.xlabel("Day")
        plt.ylim( (-1, 50) )
        plt.xlim( (-0.5, 4.5) )
        drugtreat_group_color = { "Acute high levodopa" : "green", 
                                  "Chronic saline" : "grey", 
                                  "Acute saline" : "yellow", 
                                  "Chronic high levodopa" : "red", 
                                  "Chronic low levodopa" : "blue"  }
        pd_groups = s.groupby(["MouseType", "LesionType", "DrugTreat"])
        for (mouse, lesion, drug) in pd_groups.groups.keys():
            data = s.ix[pd_groups.groups[(mouse, lesion, drug)],
                        ["Day1AIM", "Day3AIM", "Day4AIM", "Day5AIM", "Day8AIM"]].transpose()
            r = np.repeat(range(len(data.index)), len(data.columns)).reshape( (len(data.index), -1) )
            #xs = [range(len(data))
            r = r + random.random(shape(r)) * 0.15 - 0.075
            ax1.plot(r, data, "o-", color=drugtreat_group_color[drug], lw=1.5, ms=8, alpha=0.5, label=lesion + " " + drug )
            plt.xticks(range(0, 5), ["1", "3", "4", "5", "8"])
            patches, labels = ax1.get_legend_handles_labels()
            patches_list.append(patches[-1])
            label_list.append(labels[-1])
    ax_l = plt.subplot(1,3,3)
    ax_l.set_frame_on(False)
    ax_l.set_xticks([])
    ax_l.set_yticks([])
    plt.legend(patches_list[1:4], label_list[1:4], loc="upper left" )



def compareModels(cell_type_covar_ss, pd_all,  probeset):
    cur_data = pandas.DataFrame( { "expression" : pd_all.ix[probeset, cell_type_covar_ss.filenames],
                               "dose" : [int(i) for i in cell_type_covar_ss["DrugTreat"] == "Chronic high levodopa"], 
                               "AIM" : list(cell_type_covar_ss["totalAIM"].fillna(0)) } )

    aim_models = ["AIM ~ expression + C(dose) + expression:C(dose)", 
                  "AIM ~ expression + C(dose)", "AIM ~ C(dose)",  
                  "AIM ~ expression"]
    expr_models = ["expression ~ AIM*C(dose)", 
                   "expression ~ AIM + C(dose)", 
                   "expression ~ C(dose)"]

    aim_model_fits = {}
    for m in aim_models:
        aim_model_fits[m] = ols(m, cur_data).fit()
        
    expr_model_fits = {}
    for m in expr_models:
        expr_model_fits[m] = ols(m, cur_data).fit()


    return { "aim_models" : aim_models, "aim_model_fits" : aim_model_fits, 
             "expr_models" : expr_models, "expr_model_fits" : expr_model_fits  }


def compareOutlier(cp_type, mouse_id):
    outlier_selection = pd_covar.select( lambda x: pd_covar.MouseType[x] == cp_type 
                                         and pd_covar.MouseID[x] == mouse_id )
    reference_sample = pd_covar.select( lambda x: pd_covar.MouseType[x] == cp_type 
                                        and pd_covar.MouseID[x] != mouse_id 
                                        and pd_covar.DrugTreat[x] == outlier_selection.DrugTreat )
    
    outlier_data = pd_all[outlier_selection.filenames[outlier_selection.index[0]]]
    reference_data = pd_all[reference_sample.filenames]
    
    reference_mean = reference_data.mean(axis = 1)
    reference_std = reference_data.std(axis = 1)
    
    outlier_diff = outlier_data - reference_mean
    outlier_zscore = outlier_diff.div( reference_std )
    
    return pandas.DataFrame( { "difference" : outlier_diff, "zscore" : outlier_zscore } )


def outlierSelect(o, foldchange_cutoff, zscore_cutoff):
    return { "higher" : (o
                .select(lambda x: o.difference[x] > foldchange_cutoff 
                        and abs(o.zscore[x]) > zscore_cutoff).sort_index(by="zscore")
                .merge(mo430info, left_index=True, right_index=True)
                .sort_index(by="zscore", ascending=False)),
            "lower" : (o
                .select(lambda x: o.difference[x] < -foldchange_cutoff 
                        and abs(o.zscore[x]) > zscore_cutoff).sort_index(by="zscore")
                .merge(mo430info, left_index=True, right_index=True)
                .sort_index(by="zscore"))                         
            }
