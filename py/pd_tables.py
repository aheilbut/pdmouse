# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 12:31:02 2013

@author: adrian
"""

import pandas

# load tukeyHSD pvals from file
cp73_tukey = pandas.DataFrame.from_csv("/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/jan30/cp73_tukeyHSD.tab", 
                                        header=None, index_col=None, sep="\t" )
cp73_tukey.columns = ["contrast", "pval", "probeset"]
# pivot so that they're indexed by probeset
cp73_tukey_pivot = cp73_tukey.pivot_table(values="pval", rows="probeset", cols="contrast")
cp73_tukey_pivot.columns = ["CP73 TukeyHSD " + a for a in cp73_tukey_pivot.columns]

# load tukeyHSD pvals from file
cp101_tukey = pandas.DataFrame.from_csv("/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/jan30/cp101_tukeyHSD.tab",
                                        header=None, index_col=None, sep="\t" )
cp101_tukey.columns = ["contrast", "pval", "probeset"]
# pivot so that they're indexed by probeset
cp101_tukey_pivot = cp101_tukey.pivot_table(values="pval", rows="probeset", cols="contrast")
cp101_tukey_pivot.columns = ["CP101 TukeyHSD " + a for a in cp101_tukey_pivot.columns]


cp73_m = cp73_wide_info.merge( cp73_tukey_pivot, left_index=True, right_index=True )
cp101_m = cp101_wide_info.merge( cp101_tukey_pivot, left_index=True, right_index=True )

cp73_m["highlevo_vs_saline_tukey_bh"] = statsmodels.sandbox.stats.multicomp.multipletests(
        cp73_m["CP73 TukeyHSD Chronic saline-Chronic high levodopa"], 
        method="fdr_bh",
        alpha=0.10)[1]
cp73_m["lowlevo_vs_saline_tukey_bh"] = statsmodels.sandbox.stats.multicomp.multipletests(
        cp73_m["CP73 TukeyHSD Chronic saline-Chronic low levodopa"], 
        method="fdr_bh", 
        alpha=0.10)[1]
cp73_m["highlevo_vs_lowlevo_tukey_bh"] = statsmodels.sandbox.stats.multicomp.multipletests(
        cp73_m["CP73 TukeyHSD Chronic saline-Chronic low levodopa"], 
        method="fdr_bh", 
        alpha=0.10)[1]
        
cp73_m["highlevo_vs_saline_tukey_bh"] = statsmodels.sandbox.stats.multicomp.multipletests(
        cp73_m["CP73 TukeyHSD Chronic saline-Chronic high levodopa"], 
        method="fdr_bh",
        alpha=0.10)[1]
cp73_m["lowlevo_vs_saline_tukey_bh"] = statsmodels.sandbox.stats.multicomp.multipletests(
        cp73_m["CP73 TukeyHSD Chronic saline-Chronic low levodopa"], 
        method="fdr_bh", 
        alpha=0.10)[1]
cp73_m["highlevo_vs_lowlevo_tukey_bh"] = statsmodels.sandbox.stats.multicomp.multipletests(
        cp73_m["CP73 TukeyHSD Chronic saline-Chronic low levodopa"], 
        method="fdr_bh", 
        alpha=0.10)[1]        
        
        
t2_up = {
 "CP73" : 
  { "low vs. saline" : cp73_m.select( lambda x: cp73_m.ix[x, "lowlevo_vs_saline_tukey_bh"] <= 0.10 
                            and cp73_m.ix[x, ('cp73_low_vs_saline', 'fc_exp')] >= 1),
  "high vs. saline" : cp73_m.select( lambda x: cp73_m.ix[x, "highlevo_vs_saline_tukey_bh"] <= 0.10 
                            and cp73_m.ix[x, ('cp73_high_vs_saline', 'fc_exp')] >= 1),
  "high vs. low" : cp73_m.select( lambda x: cp73_m.ix[x, "highlevo_vs_lowlevo_tukey_bh"] <= 0.10 
                            and cp73_m.ix[x, ('cp73_high_vs_low', 'fc_exp')] >= 1)
  },
  "CP101" : {  
  "low vs. saline" : cp101_m.select( lambda x: cp101_m.ix[x, "lowlevo_vs_saline_tukey_bh"] <= 0.10 
                            and cp101_m.ix[x, ('cp101_low_vs_saline', 'fc_exp')] >= 1),
  "high vs. saline" :  cp101_m.select( lambda x: cp101_m.ix[x, "highlevo_vs_saline_tukey_bh"] <= 0.10 
                            and cp101_m.ix[x, ('cp101_high_vs_saline', 'fc_exp')] >= 1),
  "high vs. low" : cp101_m.select( lambda x: cp101_m.ix[x, "highlevo_vs_lowlevo_tukey_bh"] <= 0.10 
                            and cp101_m.ix[x, ('cp101_high_vs_low', 'fc_exp')] >= 1)                                           
  }
}        

t2_down = {
 "CP73" : 
  { "low vs. saline" : cp73_m.select( lambda x: cp73_m.ix[x, "lowlevo_vs_saline_tukey_bh"] <= 0.10 
                            and cp73_m.ix[x, ('cp73_low_vs_saline', 'fc_exp')] <= -1),
  "high vs. saline" : cp73_m.select( lambda x: cp73_m.ix[x, "highlevo_vs_saline_tukey_bh"] <= 0.10 
                            and cp73_m.ix[x, ('cp73_high_vs_saline', 'fc_exp')] <= -1),
  "high vs. low" : cp73_m.select( lambda x: cp73_m.ix[x, "highlevo_vs_lowlevo_tukey_bh"] <= 0.10 
                            and cp73_m.ix[x, ('cp73_high_vs_low', 'fc_exp')] <= -1)
  },
  "CP101" : {  
  "low vs. saline" : cp101_m.select( lambda x: cp101_m.ix[x, "lowlevo_vs_saline_tukey_bh"] <= 0.10 
                            and cp101_m.ix[x, ('cp101_low_vs_saline', 'fc_exp')] <= -1),
  "high vs. saline" :  cp101_m.select( lambda x: cp101_m.ix[x, "highlevo_vs_saline_tukey_bh"] <= 0.10 
                            and cp101_m.ix[x, ('cp101_high_vs_saline', 'fc_exp')] <= -1),
  "high vs. low" : cp101_m.select( lambda x: cp101_m.ix[x, "highlevo_vs_lowlevo_tukey_bh"] <= 0.10 
                            and cp101_m.ix[x, ('cp101_high_vs_low', 'fc_exp')] <= -1)                                           
  }
}        

        
count_table_2fold_up = pandas.DataFrame( {
"CP73" :                                  
    { "low vs. saline" : len(t2_up["CP73"]["low vs. saline"].symbol.unique()),
      "high vs. saline" : len(t2_up["CP73"]["high vs. saline"].symbol.unique()),
      "high vs. low" : len(t2_up["CP73"]["high vs. low"].symbol.unique() )
     },

"CP101" :                                  
    { "low vs. saline" : len(t2_up["CP101"]["low vs. saline"].symbol.unique()),
      "high vs. saline" : len(t2_up["CP101"]["high vs. saline"].symbol.unique()),
      "high vs. low" : len(t2_up["CP101"]["high vs. low"].symbol.unique() )
     }
})

count_table_2fold_down = pandas.DataFrame( {
"CP73" :                                  
    { "low vs. saline" : len(t2_down["CP73"]["low vs. saline"].symbol.unique()),
      "high vs. saline" : len( t2_down["CP73"]["high vs. saline"].symbol.unique()),
      "high vs. low" : len( t2_down["CP73"]["high vs. low"].symbol.unique() )
     },

"CP101" :                                  
   { "low vs. saline" : len(t2_down["CP101"]["low vs. saline"].symbol.unique()),
      "high vs. saline" : len( t2_down["CP101"]["high vs. saline"].symbol.unique()),
      "high vs. low" : len( t2_down["CP101"]["high vs. low"].symbol.unique() )
     }
 })


