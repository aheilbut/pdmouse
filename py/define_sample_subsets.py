# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import load_pd_data; reload(load_pd_data); from load_pd_data import *

# <codecell>

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

ss_cp101_allchronic = pd_covar.select(lambda x: (pd_covar.ix[x, "MouseType"] == "CP101" 
                                                and pd_covar.ix[x, "LesionType"] == "6-OHDA" 
                                                and (pd_covar.ix[x, "DrugTreat"] == "Chronic low levodopa" 
                                                     or pd_covar.ix[x, "DrugTreat"] == "Chronic high levodopa")
                                                and pd_covar.ix[x, "MouseID"] != 1343)) 

ss_cp101_allchronic = pd_covar.select(lambda x: (pd_covar.ix[x, "MouseType"] == "CP101" 
                                                 and pd_covar.ix[x, "LesionType"] == "6-OHDA" 
                                                 and (pd_covar.ix[x, "DrugTreat"] == "Chronic low levodopa" 
                                                        or pd_covar.ix[x, "DrugTreat"] == "Chronic high levodopa") 
                                                 and pd_covar.ix[x, "MouseID"] != 1343)) 

ss_cp73_allchronic = pd_covar.select(lambda x: (pd_covar.ix[x, "MouseType"] == "CP73" 
                                                and pd_covar.ix[x, "LesionType"] == "6-OHDA" 
                                                and (pd_covar.ix[x, "DrugTreat"] == "Chronic low levodopa" \
                                                     or pd_covar.ix[x, "DrugTreat"] == "Chronic high levodopa")
                                                and pd_covar.ix[x, "MouseID"] != 1343)) 

ss_cp73_allChronic_LDOPA = pd_covar.select(lambda x: (pd_covar.ix[x, "MouseType"] == "CP73" 
                                                and pd_covar.ix[x, "LesionType"] == "6-OHDA" 
                                                and (pd_covar.ix[x, "DrugTreat"] == "Chronic high levodopa" 
                                                    or pd_covar.ix[x, "DrugTreat"] == "Chronic low levodopa")))
                                          
ss_cp101_allChronic_LDOPA = pd_covar.select(lambda x: (pd_covar.ix[x, "MouseType"] == "CP101" 
                                                and pd_covar.ix[x, "LesionType"] == "6-OHDA" 
                                                and pd_covar.ix[x, "MouseID"] != 1343
                                                and (pd_covar.ix[x, "DrugTreat"] == "Chronic high levodopa" 
                                                     or pd_covar.ix[x, "DrugTreat"] == "Chronic low levodopa")))
                                          

cp73_chronic_saline = pd_covar.select(lambda x: pd_covar.ix[x, "MouseType"] == "CP73" 
                                        and pd_covar.ix[x, "LesionType"] == "6-OHDA" 
                                        and pd_covar.ix[x, "DrugTreat"] == "Chronic saline")
cp101_chronic_saline = pd_covar.select(lambda x: pd_covar.ix[x, "MouseType"] == "CP101" 
                                        and pd_covar.ix[x, "LesionType"] == "6-OHDA" 
                                        and pd_covar.ix[x, "DrugTreat"] == "Chronic saline")

ss_cp73_acuteHigh =  pd_covar.select(lambda x: (pd_covar.ix[x, "MouseType"] == "CP73" 
                                                and pd_covar.ix[x, "LesionType"] == "6-OHDA" 
                                                and pd_covar.ix[x, "DrugTreat"] == "Acute high levodopa"))

ss_cp73_acuteSaline =  pd_covar.select(lambda x: (pd_covar.ix[x, "MouseType"] == "CP73" 
                                                and pd_covar.ix[x, "LesionType"] == "6-OHDA" 
                                                and pd_covar.ix[x, "DrugTreat"] == "Acute saline"))

# <codecell>

ss_cp101_ascorbate_saline = pd_covar.select( lambda x: pd_covar.ix[x, "MouseType"] == "CP101" 
                                                    and pd_covar.ix[x, "DrugTreat"] == "Chronic saline" 
                                                    and pd_covar.ix[x, "LesionType"] == "Ascorbate"  )

ss_cp101_OHDA_saline = pd_covar.select( lambda x: pd_covar.ix[x, "MouseType"] == "CP101" 
                                                    and pd_covar.ix[x, "DrugTreat"] == "Chronic saline" 
                                                    and pd_covar.ix[x, "LesionType"] == "6-OHDA"  )

# <codecell>

ss_cp73_ascorbate_saline = pd_covar.select( lambda x: pd_covar.ix[x, "MouseType"] == "CP73" 
                                                    and pd_covar.ix[x, "DrugTreat"] == "Chronic saline" 
                                                    and pd_covar.ix[x, "LesionType"] == "Ascorbate"  )

ss_cp73_OHDA_saline = pd_covar.select( lambda x: pd_covar.ix[x, "MouseType"] == "CP73" 
                                                    and pd_covar.ix[x, "DrugTreat"] == "Chronic saline" 
                                                    and pd_covar.ix[x, "LesionType"] == "6-OHDA"  )

