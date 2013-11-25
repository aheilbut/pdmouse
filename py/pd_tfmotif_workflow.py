# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 14:33:41 2013

@author: adrian
"""

import sys
sys.path.append("/data/adrian/code/projects/broad/")
import pandas as pd
import pd_analysis as pda
#import alib.plots
import pandasql
import cPickle
import psycopg2
import alib.plots
import numpy as np
import pdgz as gz
from collections import namedtuple
import scipy.stats

conn = psycopg2.connect(database="mousegenome", user="adrian", password="adrian", host="localhost")


class LoadGeneMotifs(gz.GZWorkNode):

    def __init__(self):
        class PortsIn(gz.GZPortGroup):
            pass
        class PortsOut(gz.GZPortGroup):
            mm9_gene_motifs = gz.GZPort(self, "mm9_gene_motifs")
        self.i = PortsIn()
        self.o = PortsOut()
    
    def run(self):
        import pandas.io.sql as sql
        self.o.mm9_gene_motifs.val = sql.read_frame("select * from mm9_gene_motifs_distinct where max_motif_score >= 0.7", conn)
        
class LoadGeneMotifsFlexible(gz.GZWorkNode):
    def __init__(self):
        class PortsIn(gz.GZPortGroup):
            upstreamDistance = gz.GZPort(self, "upstreamDistance")
#            allowInternal = gz.GZPort(self, "allowInternal")
            downstreamDistance = gz.GZPort(self, "downstreamDistance")
            species_genome = gz.GZPort(self, "species_genome")
        self.i = PortsIn()
        
        class PortsOut(gz.GZPortGroup):
            mm9_gene_motifs = gz.GZPort(self, "mm9_gene_motifs")
        self.o = PortsOut()        
        
    def run(self):
        import pandas.io.sql as sql
        upstreamDistance = self.i.upstreamDistance.val
        downstreamDistance = self.i.downstreamDistance.val 

        if self.i.species_genome.val == "mm9":
            srmotifs_tablename = "mm9_gene_srmotifs_detail"
        elif self.i.species_genome.val == "hg19":
            srmotifs_tablename = "hg19_gene_srmotifs_detail"
        
        query = """ select genename, chrom, txstart-motif_start as tx_start_dist, cdsstart-motif_start as cds_start_dist, 
        txend-motif_start as tx_end_dist, txend-motif_end as tx_end_motif_end_dist, motif_end-motif_start,  
        gene_strand, motif_strand, motif_name, motif_score 
        FROM %(srmotifs_tablename)s WHERE motif_score >= 0.7 AND 
        (
          ( gene_strand = '+' AND (motif_start > (txstart - %(upstreamDistance)d) and motif_start < (txend + %(downstreamDistance)d) ) )
        OR
          ( gene_strand = '-' AND (motif_end > (txstart - %(downstreamDistance)d) and motif_end < (txend + %(upstreamDistance)d) )  )
        )
        """ % { "upstreamDistance" : upstreamDistance, 
                "downstreamDistance" : downstreamDistance,
                "srmotifs_tablename" : srmotifs_tablename }
        
        self.o.mm9_gene_motifs.val = sql.read_frame(query, conn)
    


class LoadCheaTable(gz.GZWorkNode):
    def __init__(self):
        class PortsIn(gz.GZPortGroup):
            pass
        class PortsOut(gz.GZPortGroup):
            chea = gz.GZPort(self, "chea")
        self.i = PortsIn()
        self.o = PortsOut()
    
    def run(self):
        chea = pd.read_table("/data/adrian/Dropbox/Data/chea/chea-background.csv", header=None, 
                             names=["chea_id", "tf", "tf_pmid", "target", 
                                    "pmid", "exp_type", "cell_type", "organism", "date"], sep=",")
        self.o.chea.val = chea.astype('unicode')
        
class CalcTFTargetCounts(gz.GZWorkNode):
    def __init__(self):
        class PortsIn(gz.GZPortGroup):
            chea = gz.GZPort(self, "chea")
        class PortsOut(gz.GZPortGroup):
            tf_target_counts = gz.GZPort(self, "tf_target_counts")
        self.i = PortsIn()
        self.o = PortsOut()
    
    def run(self):
        chea = self.i.chea.val 
        tf_target_counts = []
        for (k, g) in chea.groupby("tf"):
            tf_target_counts.append( { "tf" : k, "chea_target_count" : len( pda.uniqueSym( g.target )) } )
        tf_target_counts = pd.DataFrame.from_dict(tf_target_counts)
        tf_target_counts.index = tf_target_counts.tf
        self.o.tf_target_counts.val = tf_target_counts
        
        
class CalcMotifGeneCounts(gz.GZWorkNode):
    
    def __init__(self):
        class PortsIn(gz.GZPortGroup):
            mm9_gene_motifs = gz.GZPort(self, "mm_9_gene_motifs")
        class PortsOut(gz.GZPortGroup):
            mm9_motif_gene_counts = gz.GZPort(self, "mm9_motif_gene_counts")

        self.i = PortsIn()
        self.o = PortsOut()

    def run(self):
        mm9_gene_motifs = self.i.mm9_gene_motifs.val
        
        mm9_motif_gene_counts = []
        for (k, g) in mm9_gene_motifs.groupby("motif_name"):
            mm9_motif_gene_counts.append( { "motif_name" : k, "motif_gene_count" : len(g.genename.unique()) } )
        mm9_motif_gene_counts = pd.DataFrame.from_dict( mm9_motif_gene_counts )
        mm9_motif_gene_counts.index = mm9_motif_gene_counts.motif_name        

        self.o.mm9_motif_gene_counts.val = mm9_motif_gene_counts

        
class LoadMatTFGraph(gz.GZWorkNode):
    def __init__(self):
        class PortsIn(gz.GZPortGroup):
            pass
        
        class PortsOut(gz.GZPortGroup):
            mat_tf_graph = gz.GZPort(self, "mat_tf_graph")
            
        self.i = PortsIn()
        self.o = PortsOut()
    
    def run(self):
        self.o.mat_tf_graph.val = pd.DataFrame.from_csv("/data/adrian/Dropbox/Data/swissregulon/mat_tf_graph.tab", sep="\t")        
        
        
class WideLoad(gz.GZWorkNode):
#    inputs = []
#    outputs = ["cp_both_wide_info"]
 
    def __init__(self):                
        class PortsIn(gz.GZPortGroup):
            pass                
       
        class PortsOut(gz.GZPortGroup):
            cp_both_wide_info = gz.GZPort(self, "cp_both_wide_info")
            
        self.i = PortsIn()
        self.o = PortsOut()

    def run(self):
        result = {}        
        # load amalgamated table of statistics by probeset from pickled pandas file.  This file was generated by the pd_calc_stats notebook or script
        cp_both_wide_info = cPickle.load( open("/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/2013_july_9/cp_both_wide_info.pandas.pickle") )
        self.o.cp_both_wide_info.val = cp_both_wide_info
        return 

class AddEnrichedTFStats(gz.GZWorkNode):
    inputs = ["top_tfs", "tf_target_counts", "total_genes"]
    outputs = ["tf_enrichment_stats"]
    """   
        class PortsIn():
            enriched_tfs = Port("enriched_tfs")
            tf_target_counts = Port("tf_target_counts")
            total_genes = Port("total_genes")
        i = PortsIn()
        
        class PortsOut():
            tf_enrichement_stats = Port()
        o = PortsOut()
    """   

    def __init__(self):
        class PortsIn(gz.GZPortGroup):
            top_tfs = gz.GZPort(self, "top_tfs")
            tf_target_counts = gz.GZPort(self, "tf_target_counts")
            total_genes = gz.GZPort(self, "total_genes")
        self.i = PortsIn()
        
        class PortsOut(gz.GZPortGroup):
            tf_enrichment_stats = gz.GZPort(self, "tf_enrichment_stats")
        self.o = PortsOut()
    
    def run(self):
        top_tfs = self.i.top_tfs.val
        tf_target_counts = self.i.tf_target_counts.val
        total_genes = self.i.total_genes.val        
        
        for k in top_tfs.keys():
            if top_tfs[k] is not None:
                top_tfs[k]["pval"] = top_tfs[k].merge(tf_target_counts, left_on="tf_name", right_index=True).apply( lambda x: (
                 scipy.stats.hypergeom.sf(
                      x.ix["count"]-1, # number of differentially expressed genes in set
                      total_genes,           # total number of genes
                      x.ix["chea_target_count"],   # number of genes in current set
                      x.ix["total_genes_in_group"])),     # total number of genes in test set
                    axis=1 )    
        self.o.tf_enrichment_stats.val = top_tfs
    
    
class AddEnrichedMotifStats(gz.GZWorkNode):
    inputs = ["top_motifs", "mm9_motif_gene_counts", "total_genes"]
    outputs = ["enriched_motifs"]

    def __init__(self):
        class PortsIn(gz.GZPortGroup):
            top_motifs = gz.GZPort(self, "top_motifs")
            mm9_motif_gene_counts = gz.GZPort(self, "mm9_motif_gene_counts")
            total_genes = gz.GZPort(self, "total_genes")
        self.i = PortsIn()
        
        class PortsOut(gz.GZPortGroup):
            enriched_motifs = gz.GZPort(self, "enriched_motifs")
        self.o = PortsOut()

    def run(self):
        top_motifs = self.i.top_motifs.val
        mm9_motif_gene_counts = self.i.mm9_motif_gene_counts.val
        total_genes = self.i.total_genes.val        
        
        for k in top_motifs.keys():
            if top_motifs[k] is not None and len(top_motifs[k].index) > 0:
                top_motifs[k]["pval"] = top_motifs[k].merge(mm9_motif_gene_counts, left_on="motif_name", right_index=True).apply( lambda x: (
                 scipy.stats.hypergeom.sf(
                      x.ix["count"]-1, # number of differentially expressed genes in set
                      total_genes,           # total number of genes
                      x.ix["motif_gene_count"],   # number of genes in current set
                      x.ix["total_genes_in_group"])),     # total number of genes in test set
                    axis=1 )    
        
        self.o.enriched_motifs.val = top_motifs        
    
class GetFactorTypes(gz.GZWorkNode):
    inputs = ["cp_both_wide_info"]
    outputs = ["factortypes"]

    def __init__(self):
        class PortsIn(gz.GZPortGroup):
            cp_both_wide_info = gz.GZPort(self, "cp_both_wide_info")
    
            def __init__(self, parent_node):
                self.parent_node = parent_node
    
        class PortsOut(gz.GZPortGroup):
            factortypes = gz.GZPort(self, "factortypes")
            
            def __init__(self, parent_node):
                self.parent_node = parent_node
    
        self.i = PortsIn(self)
        self.o = PortsOut(self)

    def checkPorts(self):
        pass

    def run(self):
        self.checkPorts()
        cp_both_wide_info = self.i.cp_both_wide_info.val         
        
        factortypes = {}
        for c in cp_both_wide_info.drop(["probe_id", "symbol", "gene_name"], axis=1).columns:
            for (k, v) in dict(pda.tm.d(c)).items():
                if isinstance(k, list):
                    k = tuple(k)
                if isinstance(v, list):
                    v = tuple(v)
                if factortypes.has_key(k):
                    factortypes[k].add(v)
                else:
                    factortypes[k] = set([v])

        self.o.factortypes.val = factortypes
#        return { "factortypes" : factortypes }
        
class GetContrastPatterns(gz.GZWorkNode):    
    def __init__(self):
        class PortsIn(gz.GZPortGroup):
            cp_both_wide_info = gz.GZPort(self, "cp_both_wide_info")
            factortypes = gz.GZPort(self, "factortypes")
        self.i = PortsIn()
        
        class PortsOut(gz.GZPortGroup):
            contrastPatterns = gz.GZPort(self, "contrastPatterns")
            activePatterns = gz.GZPort(self, "activePatterns")
            patternGroups = gz.GZPort(self, "patternGroups")
        self.o = PortsOut()
    
    def run(self):
        cp_both_wide_info = self.i.cp_both_wide_info.val
        factortypes = self.i.factortypes.val
        
        def getContrasts(x):
            test_info = [["mc", "bh"], ["st", "pval"], ["tt", "welch ttest"]]
            fc_info = [["st", "fc_means"]]
            
            res = {}
            for fc_threshold in [1.0]:
                for comparison in factortypes["cmp"]:
                    for celltype in factortypes["ct"]:
                        #print comparison
                        test_query = test_info + [ ["cmp", list(comparison)] ] + [ ["ct", celltype] ]
                        fc_query = fc_info + [ ["cmp", list(comparison)] ] + [ ["ct", celltype] ]
                        #td = dict(tm.d(c))
                        if pda.tm.e(test_query) in set(cp_both_wide_info.columns):
                            res[pda.tm.e(fc_query)] = int((x[ pda.tm.e(test_query) ] < 0.10)) * x[ pda.tm.e(fc_query) ] * int(abs(x[pda.tm.e(fc_query)]) > np.log2(fc_threshold))
            return pd.Series(res)
        
        contrastPatterns = cp_both_wide_info.apply( getContrasts, axis=1 )
        activePatterns = contrastPatterns.select( lambda x: contrastPatterns.ix[x, :].abs().sum() > 0 ) 
        patternGroups = np.sign( activePatterns ).groupby( list(activePatterns.columns), axis=0 )        
        
        self.o.contrastPatterns.val = contrastPatterns
        self.o.activePatterns.val = activePatterns
        self.o.patternGroups.val = patternGroups
                        
class GetEnrichedTFs(gz.GZWorkNode):
    inputs = ["cp_both_wide_info", "chea", "fc_threshold", "factortypes"]
    outputs = ["top_tfs"]

    def __init__(self):
        class PortsIn(gz.GZPortGroup):
            cp_both_wide_info = gz.GZPort(self, "cp_both_wide_info")
            chea = gz.GZPort(self, "chea")
            fc_threshold = gz.GZPort(self, "fc_threshold")
            factortypes = gz.GZPort(self, "factortypes")
        self.i = PortsIn()
        
        class PortsOut(gz.GZPortGroup):
            top_tfs = gz.GZPort(self, "top_tfs")
        self.o = PortsOut()
        
        

    
    def run(self):
        cp_both_wide_info = self.i.cp_both_wide_info.val
        chea = self.i.chea.val
        fc_threshold = self.i.fc_threshold.val
        factortypes = self.i.factortypes.val         
        
        def getTFMatches( genegroup, combination_size=1):
            genegroup["upper_sym"] = [a.upper() if isinstance(a, str) else "-" for a in genegroup.symbol ]
            cz = chea.merge( genegroup, left_on="target", right_on="upper_sym")
            
            total_genes_in_group = len( pda.uniqueSym(genegroup.symbol) )
            
            print(len(cz.index))
            result = []
            for k, z in cz.groupby("tf"):
                result.append( { "tf_name" : k,
                                 "count" : len(pda.uniqueSym(z.symbol)),
                                 "total_genes_in_group" : total_genes_in_group,
                                 "tf_associations" : z
                                 })
            result = pd.DataFrame(result)
            result.index = result.tf_name
            return result    
                

        test_info = [["mc", "bh"], ["st", "pval"], ["tt", "welch ttest"]]
        fc_info = [["st", "fc_means"]]
        
        summary = {}
        top_tfs = {}
        
        for fc_threshold in [1.5]:
            for comparison in factortypes["cmp"]:
                for celltype in factortypes["ct"]:
                    #print comparison
                    test_query = test_info + [ ["cmp", list(comparison)] ] + [ ["ct", celltype] ]
                    fc_query = fc_info + [ ["cmp", list(comparison)] ] + [ ["ct", celltype] ]
                    #td = dict(tm.d(c))
                    if pda.tm.e(test_query) in set(cp_both_wide_info.columns):
                        print pda.tm.e(test_query)
                        current_group = cp_both_wide_info.select(
                                lambda x: cp_both_wide_info.ix[x, pda.tm.e(test_query) ] < 0.10
                                      and abs(cp_both_wide_info.ix[x, pda.tm.e(fc_query) ]) >= np.log2(fc_threshold)
                                        and isinstance(cp_both_wide_info.ix[x, "symbol" ], str) )[["probe_id", "symbol",pda.tm.e(test_query), pda.tm.e(fc_query)]] 
                        k = pda.tm.e([ ["cmp", list(comparison)] ] 
                                         + [ ["ct", celltype] ]
                                         + [ ["dir", "any"]])
                        try:
                            top_tfs[k] = getTFMatches( current_group, 1)
                        except:
                            top_tfs[k] = None
                        
                        
                        current_group = cp_both_wide_info.select(
                                lambda x: cp_both_wide_info.ix[x, pda.tm.e(test_query) ] < 0.10
                                      and cp_both_wide_info.ix[x, pda.tm.e(fc_query) ] >= np.log2(fc_threshold)
                                        and isinstance(cp_both_wide_info.ix[x, "symbol" ], str) )[["probe_id", "symbol",pda.tm.e(test_query), pda.tm.e(fc_query)]]  
                        k = pda.tm.e([ ["cmp", list(comparison)] ] 
                                         + [ ["ct", celltype] ]
                                         + [ ["dir", "up"]])                
                        
                        try:
                            top_tfs[ k ] = getTFMatches( current_group, 1)
                        except:
                            top_tfs[ k ] = None
        
                            
                        current_group = cp_both_wide_info.select(
                                lambda x: cp_both_wide_info.ix[x, pda.tm.e(test_query) ] < 0.10
                                      and cp_both_wide_info.ix[x, pda.tm.e(fc_query) ] <= -np.log2(fc_threshold)
                                        and isinstance(cp_both_wide_info.ix[x, "symbol" ], str) )[["probe_id", "symbol",pda.tm.e(test_query), pda.tm.e(fc_query)]] 
                        k = pda.tm.e([ ["cmp", list(comparison)] ] 
                                         + [ ["ct", celltype] ]
                                         + [ ["dir", "down"]])                
                        try:
                            top_tfs[ k ] = getTFMatches( current_group, 1)
                        except:
                            top_tfs[ k ] = None
                            
                    else:
                        print "--"
            self.o.top_tfs.val = top_tfs

class GetEnrichedSRMotifs(gz.GZWorkNode):
    inputs = ["cp_both_wide_info", "fc_threshold", "mm9_gene_motifs", "factortypes"]
    outputs = ["enriched_motifs"]
           
    def __init__(self):
        class PortsIn(gz.GZPortGroup):
            cp_both_wide_info = gz.GZPort(self, "cp_both_wide_info")
            fc_threshold = gz.GZPort(self, "fc_threshold")
            mm9_gene_motifs = gz.GZPort(self, "mm9_gene_motifs")
            factortypes = gz.GZPort(self, "factortypes")
        self.i = PortsIn()
        
        class PortsOut(gz.GZPortGroup):
            enriched_motifs = gz.GZPort(self, "enriched_motifs")

        self.o = PortsOut()

    def run(self):        
        cp_both_wide_info = self.i.cp_both_wide_info.val
        fc_threshold = self.i.fc_threshold.val
        mm9_gene_motifs = self.i.mm9_gene_motifs.val
        factortypes = self.i.factortypes.val        
        
        def getMotifMatches( genegroup ):
            cz= mm9_gene_motifs.merge(genegroup, left_on="genename", right_on="symbol") 
            
            total_genes_in_group = len( pda.uniqueSym(genegroup.symbol) )
            
            print(len(cz.index))
            
            result = []
            for k, z in cz.groupby("motif_name"):
                result.append( { "motif_name"  : k,
                                 "count" : len(pda.uniqueSym(z.symbol)),
                                 "total_genes_in_group" : total_genes_in_group,
                                 "motif_occurrences": z
                                }  )
                
            result = pd.DataFrame(result)
            result.index = result.motif_name
            return result
        
        test_info = [["mc", "bh"], ["st", "pval"], ["tt", "welch ttest"]]
        fc_info = [["st", "fc_medians"]]
        
        top_motifs = {}
        
        for fc_threshold in [1.5]:
            for comparison in factortypes["cmp"]:
                for celltype in factortypes["ct"]:
                    #print comparison
                    test_query = test_info + [ ["cmp", list(comparison)] ] + [ ["ct", celltype] ]
                    fc_query = fc_info + [ ["cmp", list(comparison)] ] + [ ["ct", celltype] ]
                    #td = dict(tm.d(c))
                    if pda.tm.e(test_query) in set(cp_both_wide_info.columns):
                        print pda.tm.e(test_query)
                        current_group = cp_both_wide_info.select(
                                lambda x: cp_both_wide_info.ix[x, pda.tm.e(test_query) ] < 0.10
                                      and abs(cp_both_wide_info.ix[x, pda.tm.e(fc_query) ]) >= np.log2(fc_threshold)
                                        and isinstance(cp_both_wide_info.ix[x, "symbol" ], str) )[["probe_id", "symbol",pda.tm.e(test_query), pda.tm.e(fc_query)]]
                        k = pda.tm.e([ ["cmp", list(comparison)] ] 
                                         + [ ["ct", celltype] ]
                                         + [ ["dir", "any"]])
                        try:
                            top_motifs[k] = getMotifMatches( current_group )
                        except:
                            print sys.exc_info()
                            top_motifs[k] = None
                        
                        
                        current_group = cp_both_wide_info.select(
                                lambda x: cp_both_wide_info.ix[x, pda.tm.e(test_query) ] < 0.10
                                      and cp_both_wide_info.ix[x, pda.tm.e(fc_query) ] >= np.log2(fc_threshold)
                                        and isinstance(cp_both_wide_info.ix[x, "symbol" ], str) )[["probe_id", "symbol",pda.tm.e(test_query), pda.tm.e(fc_query)]]  
                        k = pda.tm.e([ ["cmp", list(comparison)] ] 
                                         + [ ["ct", celltype] ]
                                         + [ ["dir", "up"]])                
                        
                        try:
                            top_motifs[ k ] = getMotifMatches( current_group )
                        except:
                            print sys.exc_info()
                            top_motifs[ k ] = None
        
                            
                        current_group = cp_both_wide_info.select(
                                lambda x: cp_both_wide_info.ix[x, pda.tm.e(test_query) ] < 0.10
                                      and cp_both_wide_info.ix[x, pda.tm.e(fc_query) ] <= -np.log2(fc_threshold)
                                        and isinstance(cp_both_wide_info.ix[x, "symbol" ], str) )[["probe_id", "symbol",pda.tm.e(test_query), pda.tm.e(fc_query)]]  
                        k = pda.tm.e([ ["cmp", list(comparison)] ] 
                                         + [ ["ct", celltype] ]
                                         + [ ["dir", "down"]])                
                        try:
                            top_motifs[ k ] = getMotifMatches( current_group )
                        except:
                            print sys.exc_info()
                            top_motifs[ k ] = None
        self.o.enriched_motifs.val = top_motifs
