# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

reload(gz)

# <codecell>

from IPython.external.mathjax import install_mathjax

# <codecell>

gz.GZConst(1.5).run()

# <codecell>

import networkx as nx
import statsmodels.sandbox.stats.multicomp
import IPython

# <codecell>

import scipy.cluster
import fastcluster

# <codecell>

import alib.plots

# <codecell>

reload(gz)
reload(tf)
reload(tf.pda)

# <codecell>

c.results.keys()

# <codecell>

pd_tf_enrich = c.results[ ('AddEnrichedTFStats', 1)]["tf_enrichment_stats"]

# <codecell>

c.results[ ('AddEnrichedMotifStats', 1)].keys()

# <codecell>

pd_motif_enrich = c.results[ ('AddEnrichedMotifStats', 1)]["enriched_motifs"]

# <codecell>

for k in pd_motif_enrich.keys():
    if pd_motif_enrich[k] is not None:
        a = pd_motif_enrich[k].motif_occurrences
        pd_motif_enrich[k]["targets"] = [str( list( a[x].symbol.unique() ) ) for x in a.index]
        pd_motif_enrich[k] = pd_motif_enrich[k].merge( mm9_motif_gene_counts, left_on="motif_name", right_on="motif_name")

# <codecell>

pd_motif_enrich['[["cmp", ["6-OHDA, chronicSaline", "Ascorbate, chronicSaline"]], ["ct", "cp73"], ["dir", "up"]]'].sort

# <codecell>

for k in pd_motif_enrich.keys():
    if pd_motif_enrich[k] is not None:
        pd_motif_enrich[k]["bh pval"] = statsmodels.sandbox.stats.multicomp.multipletests(  pd_motif_enrich[k].pval )[1]
        
#statsmodels.sandbox.stats.multicomp.multipletests( 
#    pd_motif_enrich['[["cmp", ["6-OHDA, chronicSaline", "Ascorbate, chronicSaline"]], ["ct", "cp73"], ["dir", "up"]]'].pval )[1]

# <codecell>

motifs = c.results[ ('AddEnrichedMotifStats', 1)]["enriched_motifs"]
chea_binding = c.results[ ("AddEnrichedTFStats", 1) ]["tf_enrichment_stats"]

# <codecell>

contrastPatterns = c.results[ ("GetContrastPatterns",1)]["contrastPatterns"]

# <codecell>

mm9_motif_gene_counts[0:10]

# <codecell>

mm9_motif_gene_counts[0:10]

# <codecell>

hg19_motif_gene_counts = c.results[("CalcMotifGeneCounts",1)]["mm9_motif_gene_counts"]

# <codecell>

mm9_motif_gene_counts = c.results[("CalcMotifGeneCounts",1)]["mm9_motif_gene_counts"]

# <codecell>

tf_target_counts = c.results[("CalcTFTargetCounts",1)]["tf_target_counts"]

# <codecell>

tf_target_counts

# <codecell>

mm9_motif_gene_counts

# <codecell>

motifs[motifs.keys()[0]].merge(mm9_motif_gene_counts, left_index=True, right_index=True)

# <codecell>

k = chea_binding.keys()[0]

# <codecell>

k

# <codecell>

a = motifs[k]

# <codecell>

chea_binding[k]

# <codecell>

corroborated = cur_group.merge( mj, how='outer', left_on="tf_name", right_on="tf_u", suffixes=("_tf", "_motif"))

# <codecell>

f = open("/data/adrian/Dropbox/tf_motif_intersection_sep4.html", "w")

def listformatter(x):
    return "<div style='max-height: 200px; overflow:auto'>" + x + "</div>"

for k in chea_binding.keys():
    if chea_binding[k] is not None:
        cur_group = chea_binding[k].merge(tf_target_counts, left_index=True, right_index=True)
        mj = (motifs[k].merge(mm9_motif_gene_counts, left_index=True, right_index=True, suffixes=('','_y'))
              .merge( motif_tf["mat_tf_graph"], left_on="motif_name", right_on="motif_name" ) )
#              
#              )
#        print mj.columns
        corroborated = cur_group.merge( mj, how='inner', left_on="tf_name", right_on="tf_u", suffixes=("_tf", "_motif"))
    #    corroborated = cur_group.select( lambda x: cur_group.ix[x, "tf_name"] in set( map(lambda x: x.upper(), list(mj.tf_symbol)) ) )
        corroborated["overlap_count"] = corroborated.apply( lambda x: len( set( x.ix["tf_associations"].symbol ).intersection( x.ix["motif_occurrences"].symbol )) if type(x.ix["motif_occurrences"]) != float and type(x.ix["tf_associations"]) != float else None, axis=1 )
        corroborated["overlap_genes"] = corroborated.apply( lambda x: str(sort(list( set( x.ix["tf_associations"].symbol ).intersection( x.ix["motif_occurrences"].symbol )))) if type(x.ix["motif_occurrences"]) != float and type(x.ix["tf_associations"]) != float else None, axis=1 )
        
        
#        print k
        f.write("<hr>")
        f.write("<h2>" + k + "</h2>")
        f.write( (corroborated[["count_tf", "tf_name", "tf_u", "total_genes_in_group_tf", "chea_target_count", "pval_tf", 
                             "motif_name", "count_motif", "motif_gene_count", "pval_motif", "overlap_count", "overlap_genes"]]
        .select(lambda x: corroborated.ix[x, "pval_tf"] < 0.20 or corroborated.ix[x, "pval_motif"] < 0.20)
        .sort_index(by="pval_tf", ascending=True)).to_html(escape=False, formatters={"overlap_genes" : listformatter }) )

# <codecell>

chea_binding[k]

# <codecell>

(cur_group.drop("tf_associations", axis=1)
    .sort_index(by="pval", ascending=True)
    .select(lambda x: cur_group.ix[x, "pval"] < 0.10).to_json(orient='records') )

# <codecell>

chea_binding[k].tf_associations[0]

# <codecell>

import hashlib

# <codecell>

[ hashlib.new('sha1', m).hexdigest() for m in motifs[k].motif_name ]

# <codecell>

motifs[k]["divid"] = [ hashlib.new('sha1', m).hexdigest() for m in motifs[k].motif_name ]

# <codecell>

(motifs[k].merge(mm9_motif_gene_counts, left_index=True, right_index=True, suffixes=('','_y'))
     .select(lambda x: motifs[k].ix[x, "pval"] < 0.20)
     .sort_index(by="pval", ascending=True)
     .drop("motif_occurrences", axis=1).to_json(orient='records') )

# <codecell>

mj.drop("motif_occurrences", axis=1).to_json(orient='records')

# <codecell>

corroborated.index

# <codecell>

corroborated["motif_divid"] = [ hashlib.new('sha1', m).hexdigest() for m in corroborated.motif_name ]

# <codecell>

corroborated[["tf_name", "motif_divid"]].to_json(orient="records")

# <codecell>

f.close()

# <codecell>

pd.DataFrame.to_html(

# <codecell>

import pandas as pd

# <codecell>

pd.set_printoptions(max_colwidth=100000)

# <codecell>

tftarg = set( [a.upper() for a in corroborated.tf_associations[43].symbol] )
motiftarg = set( [a.upper() for a in corroborated.motif_occurrences[43].symbol] )

# <codecell>

#for i in corroborated.index:
corroborated["overlap_size"] = corroborated.apply( lambda x: len( set( x.ix["tf_associations"].symbol ).intersection( x.ix["motif_occurrences"].symbol )), axis=1 )
corroborated["overlap_genes"] = corroborated.apply( lambda x: str(list( set( x.ix["tf_associations"].symbol ).intersection( x.ix["motif_occurrences"].symbol ))), axis=1 )

# <codecell>

list( set( x.ix["tf_associations"].symbol ).intersection( x.ix["motif_occurrences"].symbol ))

# <codecell>

sort(["a", "b", "z", "c"])

# <codecell>

corroborated.tf_associations[43][0:10].target

# <codecell>

corroborated.motif_occurrences[43][0:10]

# <codecell>

 tftarg.intersection(motiftarg) 

# <codecell>

IPython.display.HTML( (corroborated[["count_tf", "tf_name", "total_genes_in_group_tf", "pval_tf", 
                             "motif_name", "count_motif", "pval_motif", "overlap_size"]]
        .select(lambda x: corroborated.ix[x, "pval_tf"] < 0.20 and corroborated.ix[x, "pval_motif"] < 0.20)
        .sort_index(by="pval_tf", ascending=True)).to_html() )

# <codecell>

corroborated.columns

# <codecell>

(corroborated[["count_tf", "tf_name", "pval_tf", "motif_name", "pval_motif"]]
    .select(lambda x: corroborated.ix[x, "pval_tf"] < 0.20 and corroborated.ix[x, "pval_motif"] < 0.20)
    .sort_index(by="pval_tf", ascending=True))

# <codecell>


len( set( chea_binding[motifs.keys()[0]].tf_name).intersection( set( map(lambda x: x.upper(), list(mj.tf_symbol)) ) ))

# <codecell>

w.worknodes["GetEnrichedSRMotifs"][1].inputs

# <codecell>

motif_tf = c.results[ ('LoadMatTFGraph', 1)]

# <codecell>

motif_tf

# <codecell>

cPickle.dump(open("/data/adrian/motif_tf_graph", "w"))

# <codecell>

mj = motifs[motifs.keys()[0]].merge( motif_tf["mat_tf_graph"], left_on="motif_name", right_on="motif_name" )

# <codecell>

mj[["count", "motif_name", "tf_symbol"]][0:10]

# <codecell>

m = c.results[('LoadGeneMotifs',1)]["mm9_gene_motifs"]

# <codecell>

c.results.keys()

# <codecell>

t = c.results[('LoadCheaTable', 1)]["chea"]

# <codecell>

motifs = c.results[ ('GetEnrichedSRMotifs', 1)]["enriched_motifs"]

# <codecell>

t.select(lambda x: t.ix[x, "target"] == "DRD1A")

# <codecell>

m.select( lambda x: m.ix[x, "genename"] == "Drd1a" )
    .merge( motif_tf["mat_tf_graph"], left_on="motif_name", right_on="motif_name" )

# <codecell>

motif_tf["mat_tf_graph"]["tf_u"] = [a.upper() for a in motif_tf["mat_tf_graph"].tf_symbol]

# <codecell>

td = t.select(lambda x: t.ix[x, "target"] == "DRD1A")

# <codecell>

(m.select( lambda x: m.ix[x, "genename"] == "Drd1a" )
 .merge( motif_tf["mat_tf_graph"], left_on="motif_name", right_on="motif_name" ))

# <codecell>

(m.select( lambda x: m.ix[x, "genename"] == "Drd1a" )
 .merge( motif_tf["mat_tf_graph"], left_on="motif_name", right_on="motif_name" )
 .merge( td, left_on="tf_u", right_on="tf")
  )

# <codecell>

motifs[motifs.keys()[0]]

# <codecell>

for (k, v) in c.results.items():
    print k, v.keys()

# <codecell>

figsize(20, 20)
nx.draw_graphviz(g, 'dot', node_shape="s", node_color='lightblue', node_size=500)

# <codecell>

|nx.topological_sort(g)[0]

# <codecell>

dict( g.nodes(data=True) )[nx.topological_sort(g)[0]]["worknode"].run()

# <codecell>

reload(gz)

# <codecell>

w.wires

# <codecell>

a = LoadGeneMotifs()

# <codecell>

import alib.plots

# <codecell>

class EnrichmentFlow(GZWorkFlow):
    """
    functional view:
    GetEnrichedTFs( cp_both_wide_info = WideLoad.cp_both_wide_info,
                    chea = LoadCheaTable.chea,
                    fc_threshold = 1.5,
                    factortypes = GetFactorTypes.factortypes
                    )



    graph view:
    
    WideLoad.cp_both_wide_info -> GetEnrichedTFs.cp_both_wide_info
    LoadCheaTable.chea -> GetEnrichedTFs.chea
    const(1.5) -> GetEnrichedTFs.fc_threshold 
    
    use WideLoad, LoadCheaTable, GetFactorTypes
    
    connect(WideLoad.cp_both_wide_info, GetEnrichedTFs.cp_both_wide_info)
    connect(LoadCheaTable.chea, GetEnrichedTFs.chea)
    connect(WideLoad.cp_both_wide_info, GetFactorTypes.cp_both_wide_info)
    connect(GetFactorTypes.factortypes, GetEnrichedTFs.factortypes)
    connect(Constant(1. 5), GetEnrichedTFs.fc_threshold)

  """

# <codecell>

contrastPatterns = c.results[("GetContrastPatterns",1)]["contrastPatterns"]
activePatterns = c.results[("GetContrastPatterns",1)]["activePatterns"]
patternGroups = c.results[("GetContrastPatterns",1)]["patternGroups"]

# <codecell>

contrastPatterns.drop(contrastPatterns.columns[0:2], axis=1)

# <codecell>

contrastPat_noAcute = contrastPatterns.drop(contrastPatterns.columns[0:2], axis=1)
activePatterns = contrastPat_noAcute.select( lambda x: contrastPat_noAcute.ix[x, :].abs().sum() > 0 ) 

# <codecell>

for i, c in enumerate(activePatterns.columns):
    print i, c

# <codecell>

activeSorted = activePatterns.sort_index( by=[activePatterns.columns[3], 
                                              activePatterns.columns[4], 
                                              activePatterns.columns[0],
                                              activePatterns.columns[5],                                              
                                              activePatterns.columns[6],                                              
                                              activePatterns.columns[7],                                              
#                                              activePatterns.columns[8],                                              
#                                              activePatterns.columns[9],                                              
                                              activePatterns.columns[2],                                              
                                              activePatterns.columns[1]                                          
                                              ]) #.drop([activePatterns.columns[0], activePatterns.columns[1]], axis=1)

# <codecell>

import pd_analysis as pda
import pandas as pd

# <codecell>

gene_reps = []
for symbol, g in pda.mo430symbol.groupby("symbol").groups.items():
    gene_reps.append(  { "symbol" : symbol, "probe_id" : g[0] } )
gene_reps = pd.DataFrame( gene_reps )

# <codecell>

activeSorted = activeSorted.merge(gene_reps, left_index=True, right_on="probe_id")

# <codecell>

activeSorted.index= activeSorted.probe_id

# <codecell>

activeSorted = activeSorted.drop( "probe_id", axis=1)
activeSorted = activeSorted.drop( "symbol", axis=1)

# <codecell>

sign( activeSorted[0:5].as_matrix()[0:5,0:5] )  

# <codecell>

activePatterns["b"] = 3

# <codecell>

activeSorted["b"] = 3

# <codecell>

activeSorted.as_matrix() * [5, 2, 2, 5, 5, 2, 2, 2, 2, 2, 0]

# <codecell>

distances = scipy.cluster.hierarchy.distance.pdist( activeSorted.as_matrix(), # * [3, 5, -4, 6, -3, 4, 1, 1, 1], 
                                                   metric="cosine")

# <codecell>

import sys
sys.setrecursionlimit(10000)

# <codecell>

rowY = fastcluster.linkage(distances)

# <codecell>

rowZ = scipy.cluster.hierarchy.dendrogram(rowY, orientation='right', no_plot=True)

# <codecell>

col_distances = scipy.cluster.hierarchy.distance.pdist( activeSorted.as_matrix().transpose() )
col_rowY = fastcluster.linkage(col_distances)
col_rowZ = scipy.cluster.hierarchy.dendrogram(col_rowY, orientation='right', no_plot=True)

# <codecell>

sum( abs(activeSorted[ activeSorted.columns[0:-1] ].sum(axis=1)) > 0 )

# <codecell>

for c in activePatterns.columns[[9, 7, 5, 0, 3, 1, 4, 6, 8]]:
     print c

# <codecell>

activeSorted

# <codecell>

figsize(40, 8)
plotData =  activeSorted.as_matrix()[rowZ['leaves'], :]
plotData = plotData[:, [7, 5, 3, 1, 2, 4, 6] ] #[9, 7, 5, 0, 3, 1, 4, 6, 8]
imshow( plotData.transpose(), 
           interpolation='nearest', aspect=200, 
           cmap=alib.plots.my_cmap, vmin=-2, vmax=2)
yticks([0, 1, 2, 3, 4, 5, 6], 
       ["Dop. Depletion - dSPN ", 
        "Chronic Low vs. Saline - dSPN ", 
        "Chronic High vs. Saline - dSPN ", # "Acute High vs. Saline - dSPN ",
        "Chronic High vs. Chronic Low - dSPN", # "Acute vs Chronic - dSPN ",
        "Chronic High vs. Saline - iSPN ", 
        "Chronic Low vs. Saline - iSPN ", 
        "Dop Depletion - iSPN "], fontsize=16) 
xticks(arange(0, 9000, 500))
grid(which="both", ls='solid', color='lightgray', alpha=0.5 )
tight_layout()
savefig("/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/2013_10_25/figures/contrast_barcode_geneRep.pdf")

# <codecell>

figsize(30, 10)
plotData =  activeSorted.as_matrix()[rowZ['leaves'], :][0:1000, :]
plotData = plotData[:, [9, 7, 5, 0, 3, 1, 4, 6, 8] ]
imshow( plotData.transpose(), 
           interpolation='nearest', aspect=20, 
           cmap=alib.plots.my_cmap, vmin=-1, vmax=1)
yticks([0, 1, 2, 3, 4, 5, 6, 7, 8], 
       ["Dop. Depletion - dSPN ", "Chronic Low vs. Saline - dSPN ", 
        "Chronic High vs. Saline - dSPN ", "Acute High vs. Saline - dSPN ",
        "Chronic High vs. Chronic Low - dSPN", "Acute vs Chronic - dSPN ",
        "Chronic High vs. Saline - iSPN ", "Chronic Low vs. Saline - iSPN ", 
        "Dop Depletion - iSPN "], fontsize=16)
#xticks(arange(0, 13000, 500))
grid(which="both", ls='solid', color='lightgray', alpha=0.3 )
tight_layout()

# <codecell>

matplotlib.__version__

# <codecell>

activeSorted[0:10]

# <codecell>

alib.plots.clusterHeatmap( activeSorted[0:1000],"", None, None, 
                              cluster_rows=False, cluster_columns=False, width=6, height=15)

# <codecell>

alib.plots.clusterHeatmap( pd.DataFrame( patternGroups.groups.keys(), columns=activePatterns.columns), 
                          "", None, None, width=6, height=15,
                          cluster_rows= True, cluster_columns=True, distmethod='correlation')

# <codecell>

figsize(10, 15)
alib.plots.clusterHeatmap(activePatterns, "", None, None, height=15, width=8,
                          cluster_rows=True, cluster_columns=True,
                          distmethod='euclidean')

# <codecell>

reload(alib.plots)

# <codecell>

import fastcluster

# <codecell>

import scipy

# <codecell>

distances = scipy.cluster.hierarchy.distance.pdist(activePatterns.values[0:5000], "correlation") 

# <codecell>

rowY = fastcluster.linkage( distances )

# <codecell>

rowZ = scipy.cluster.hierarchy.dendrogram(rowY, orientation='right', no_plot=True)

# <codecell>

activePatterns
                    

# <codecell>

ks = top_tfs.keys()
ks.sort()
for k in ks:
    print k
    if top_tfs[k] is None:
        print "n/a"
    else:
        print top_tfs[k][0:10]

# <codecell>

import scipy.stats

# <codecell>

scipy.stats.hypergeom.sf( 33, 22000, 5042, 70) 

# <codecell>

total_genes = 22000

# <codecell>

p = overlaps.apply( lambda x: (
      scipy.stats.hypergeom.sf(
          x.ix["match_count"]-1, # number of differentially expressed genes in set
          total_genes,           # total number of genes
          x.ix["size of set"],   # number of genes in current set
          len( test_set ))),     # total number of genes in test set
                               axis=1)
    p = pd.DataFrame(p, columns=["hypergeom p-val"])

# <codecell>

top_motifs[k]

# <codecell>

for k in top_motifs.keys():
    if top_motifs[k] is not None and len(top_motifs[k].index) > 0:
        print k
        print top_motifs[k].drop("motif_occurrences", axis=1).merge(mm9_motif_gene_counts, 
                                                left_on="motif_name", 
                                                right_index=True).sort_index(by="pval", ascending=True)[0:20]

# <codecell>

top_tfs.items()[0]

# <codecell>


# <codecell>

import scipy.stats
total_genes = 22000

# <codecell>

top_tfs.items()[0]

# <codecell>

"_".join( [str(z[1]) for z in pda.tm.d(k)] )

# <codecell>

e = pd.ExcelWriter("/data/adrian/Dropbox/pd_tf_enrichment_20130820.xls")

# <codecell>

e = open("/data/adrian/Dropbox/ptables/pd_tf_enrichment_20130822.html","w")

# <codecell>

out_path = "/data/adrian/Dropbox/ptables"

# <codecell>

pd.set_printoptions(max_colwidth=10000)

# <codecell>

top_tfs = chea_binding

# <codecell>

for k in top_tfs.keys():
    if top_tfs[k] is not None:
        print k
        t = top_tfs[k].sort_index(by="pval", ascending=True)[0:30].merge(tf_target_counts, left_on="tf_name", right_index=True)
        e.write("<hr>")
        e.write("<h2>" + k + "</h2>")
        
        m = top_motifs[k].merge(mm9_motif_gene_counts, left_on="motif_name", right_index=True).sort_index(by="pval", ascending=True)[0:30].drop(["motif_name_x", "motif_name_y"], axis=1)
        m["lists"] = m.apply( lambda x: "<a href='" + hashlib.new('sha1', k + x["motif_name"]).hexdigest() + ".html" + "'>list</a>",  axis=1)
        t["lists"] = t.apply( lambda x: "<a href='" + hashlib.new('sha1', k + x["tf_name"]).hexdigest() + ".html" + "'>list</a>",  axis=1)

        e.write("<table>")
        e.write("<tr><th> top 30 overrepresented associations from CHEA </th> <th> top 30 over-represented tf binding motifs (SwissRegulon) </th></tr>")
        e.write("<tr><td>" + t.drop("tf_associations", axis=1).to_html(escape=False) + "</td>" + "<td style='padding-left: 25px'>" + m.drop("motif_occurrences", axis=1).to_html(escape=False) + "</td></tr></table>")

        for tf_name in t.tf_name:
            match_table_file = open( out_path + "/" + hashlib.new('sha1', k + tf_name).hexdigest() + ".html", "w")
            match_table_file.write( t.ix[tf_name, "tf_associations"].sort_index(by="symbol").to_html() )
            match_table_file.close()            
        
        for motif_name in m.index:
            match_table_file = open( out_path + "/" + hashlib.new('sha1', k + motif_name).hexdigest() + ".html", "w")
            match_table_file.write( m.ix[motif_name, "motif_occurrences"].sort_index(by="genename").to_html() )
            match_table_file.close()

# <codecell>

(top_tfs[k].sort_index(by="pval", ascending=True)
    .merge(tf_target_counts, left_on="tf_name", right_index=True).drop("tf_associations", axis=1).to_json(orient='records'))

# <codecell>

top_motifs[k]

# <codecell>

m = motifs[k].drop("motif_occurrences",axis=1)

# <codecell>

map(lambda x: dict(zip(motifs[k].columns, x)), motifs[k].drop("motif_occurrences",axis=1)[0:10].to_records(index=False))

# <codecell>

intop_tfs[k]

# <codecell>

top_motifs = motifs

# <codecell>

top_motifs[k].drop("motif_occurrences", axis=1).to_json(orient='records')

# <codecell>

top_tfs['[["cmp", ["6-OHDA, chronicSaline", "Ascorbate, chronicSaline"]], ["ct", "cp73"], ["dir", "up"]]']

# <codecell>

m.apply( lambda x: "<a href='" + hashlib.new('sha1', k + x["motif_name"]).hexdigest() + ".html" + "'>list</a>",  axis=1)

# <codecell>

e.close()

# <codecell>

top_motifs = top_motifs.sort_index(by="pval", ascending=True)

# <codecell>

top_motifs[ '[["cmp", ["6-OHDA, chronicLow", "6-OHDA, chronicSaline"]], ["ct", "cp101"], ["dir", "down"]]' ].sort_index(by="pval", ascending=True).drop("motif_occurrences", axis=1)

# <codecell>

k

# <codecell>

import hashlib

# <codecell>

table_lists = collections.OrderedDict()

# <codecell>

hashlib.new('sha1', k + "MAFB.p2").hexdigest()

# <codecell>

h.hexdigest()

# <codecell>

top_motifs[ k ].ix["MAFB.p2", "motif_occurrences"]

# <codecell>

patterns = cp_both_wide_info[ [c for c in pda.tm.m([("cmp"), 
          ("st", "pval"), 
         ("mc", "bh"), 
          ("tt", "welch ttest")], cp_both_wide_info.columns)] ]  < 0.10

# <codecell>

[c for c in pda.tm.m([("cmp"), 
          ("st", "pval"), 
          ("mc", "bh"), 
          ("tt", "welch ttest")], cp_both_wide_info.columns)]

# <codecell>

table_fc_threshold = log2(1.5)
cn_cp73_bh_ttest_pval = pda.tm.e([("cmp", ["6-OHDA, chronicSaline", "Ascorbate, chronicSaline"]), 
        ("ct", "cp73"), ("mc", "bh"), ("st", "pval"), ("tt", "welch ttest")])
cn_cp73_median_foldchange = pda.tm.e([  ("cmp", ["6-OHDA, chronicSaline", "Ascorbate, chronicSaline"]), 
        ("ct", "cp73"), ("st", "fc_medians")])

cp73_dopamine_depletion_up = (cp_both_wide_info.select(
    lambda x: cp_both_wide_info.ix[x, cn_cp73_bh_ttest_pval] < 0.10
    and cp_both_wide_info.ix[x, cn_cp73_median_foldchange] >= table_fc_threshold
))[["probe_id", "symbol", "gene_name", 
    cn_cp73_bh_ttest_pval,
    cn_cp73_median_foldchange
    ]].sort_index(by=cn_cp73_median_foldchange, ascending=False)

# <codecell>

cn_cp73_chronicHigh_bh_ttest_pval = pda.tm.e([  ("cmp", ["6-OHDA, chronicHigh", "6-OHDA, chronicSaline"]),
            ("ct", "cp73"),
            ("mc", "bh"),
            ("st", "pval"),
            ("tt", "welch ttest")
        ])

cn_cp73_chronicHigh_foldchange = pda.tm.e([  
            ("cmp", ["6-OHDA, chronicHigh", "6-OHDA, chronicSaline"]),
            ("ct", "cp73"),
            ("st", "fc_medians")
        ])

cp73_chronicHigh_up = (cp_both_wide_info.select(
    lambda x: (cp_both_wide_info.ix[x, cn_cp73_chronicHigh_bh_ttest_pval] < 0.10)
    and cp_both_wide_info.ix[x, cn_cp73_chronicHigh_foldchange] >= table_fc_threshold
))[["symbol", cn_cp73_chronicHigh_bh_ttest_pval, cn_cp73_chronicHigh_foldchange]]

cp73_chronicHigh_down = (cp_both_wide_info.select(
    lambda x: (cp_both_wide_info.ix[x, cn_cp73_chronicHigh_bh_ttest_pval] < 0.10)
    and cp_both_wide_info.ix[x, cn_cp73_chronicHigh_foldchange] <= -table_fc_threshold
))[["symbol", cn_cp73_chronicHigh_bh_ttest_pval, cn_cp73_chronicHigh_foldchange]]

# <codecell>

[c for c in pda.tm.m([("cmp"), 
          ("st", "pval"), 
          ("mc", "bh"), 
          ("tt", "welch ttest")], cp_both_wide_info.columns)]

# <codecell>

patterns.columns = [[pda.tm.e(val) for (t, val) in pda.tm.d( c ) if t == "cmp"][0] for c in patterns.columns]

# <codecell>

def patternGroup(x):
    r = patterns.ix[ x, :]
    return int("".join([str(int(a)) for a in r]), 2)

# <codecell>

g = patterns.groupby(patternGroup)

# <codecell>

len(g.groups.items())

# <codecell>

for (k, v) in g.groups.items():
    print k, len(v)

# <codecell>

figsize(10,10)
imshow( g.get_group(656).as_matrix(), interpolation="nearest", cmap=get_cmap("gray"))

# <codecell>

targets = cz.target.unique()
tfs = cz.tf.unique()

# <codecell>

czm = zeros( [len(targets), len(tfs)] )
target_index = dict( (a[1], a[0]) for a in enumerate(targets) )
tf_index = dict( (a[1], a[0]) for a in enumerate(tfs) )

# <codecell>

pandasql.sqldf("select tf, count(*) c from cz group by tf order by c desc", locals())[0:10]

# <codecell>

for i in cz.index:
    czm[target_index[cz.ix[i, "target"]], tf_index[cz.ix[i, "tf"]]] = 1

# <codecell>

czm_df = pandas.DataFrame(data=czm, index=target_index, columns=tf_index)

# <codecell>


# <codecell>


# <codecell>

import IPython

# <codecell>

IPython.display.HTML( cz[0:5].to_html() )

# <codecell>

cp73_dopamine_depletion_up

# <codecell>

import itertools

# <codecell>

177 * 177 * 177

# <codecell>

len( cz.groupby("target").groups.items() )

# <codecell>

len(  list( itertools.combinations( cz.ix[g, "tf"].unique(), 3) ) )

# <codecell>

def countTFcombos( genegroup, combination_size ):
    z = [a.upper() for a in set( uniqueSym( genegroup.symbol ))]
    z_set = set(z)
    cz= chea.select( lambda x: chea.ix[x, "target"] in z_set)
    
    print(len(cz.index))
    
    pair_counts = {}
    for (target_group, g) in cz.groupby("target").groups.items():
        group_tfs = cz.ix[g, "tf"].unique()
        group_tfs.sort()
#        for combo in itertools.combinations( group_tfs , combination_size):
#            pair_counts[combo] = pair_counts.setdefault(combo, 0) + 1
        for tf in group_tfs:
            pair_counts[tf] = pair_counts.setdefault(tf, 0) + 1
    result = pd.DataFrame( pair_counts.items(), columns=["tf", "count"]).sort_index(by="count", ascending=False) 
    result["total_genes_in_group"] = len(z)
    #ret
    return result

# <codecell>

def countMotifCombos( genegroup, combination_size ):
    z = [a for a in set( uniqueSym( genegroup.symbol ))]
    z_set = set(z)
    cz= mm9_gene_motifs.select( lambda x: mm9_gene_motifs.ix[x, "genename"] in z_set)
    
    print(len(cz.index))
    
    pair_counts = {}
    for (target_group, g) in cz.groupby("genename").groups.items():
        group_motifs = cz.ix[g, "motif_name"].unique()
        group_motifs.sort()
#        for combo in itertools.combinations( group_tfs , combination_size):
#            pair_counts[combo] = pair_counts.setdefault(combo, 0) + 1
        for motif_name in group_motifs:
            pair_counts[motif_name] = pair_counts.setdefault(motif_name, 0) + 1
    result = pd.DataFrame( pair_counts.items(), columns=["motif_name", "count"]).sort_index(by="count", ascending=False) 
    result["total_genes_in_group"] = len(z)
    #ret
    return result



# <codecell>

table_text = "<table>"
for row in r.sort_index(by="pval", ascending=True).index[0:30]:
    print row
    table_text += "<tr>"
    table_text += "<td> <a onclick='$(\"#detail_%d\").toggle();'><b>+</b></a> </td>" % row
    table_text += "<td>" + r.ix[row, "motif_name"] + "</td>"
    table_text += "<td>" + str(r.ix[row, "count"]) + "</td>"
    table_text += "<td>" + str(r.ix[row, "pval"]) + "</td>"
    table_text += "<td>" + str(r.ix[row, "total_genes_in_group"]) + "</td>"
    table_text += "</tr>"
    table_text += "<tr><td colspan=4><div id='detail_%d' style='height: 200px; width=400px; overflow:auto; display: none'>" % row 
    table_text += r.ix[row, "motif_occurrences"].to_html() + "</div></td></tr>"
table_text += "</table>"

# <codecell>

f = open("/data/adrian/Dropbox/motif_table_aug22.html", "w")

# <codecell>

f.write(table_text)

# <codecell>

r.ix[1, "motif_occurrences"].sort_index(by=r.ix[1, "motif_occurrences"].columns[8], ascending=False)[0:50]

# <codecell>

test.merge(tf_target_counts, left_on="tf", right_index=True)[0:10]

# <codecell>

pd.set_printoptions(max_columns=100, max_rows=100)

# <codecell>

IPython.display.HTML( mm9_motif_gene_counts.sort_index(by="motif_gene_count", ascending=False).to_html() )

# <codecell>

tf_target_counts.ix["CREM", :]

# <codecell>

down_2 = pd.DataFrame( countTFcombos(cp73_chronicHigh_down).items(), 
             columns=["pair", "c"] ).sort_index(by="c", ascending=False)

# <codecell>

up_2 = pd.DataFrame( countTFcombos(cp73_chronicHigh_up).items(), 
             columns=["pair", "c"] ).sort_index(by="c", ascending=False)

# <codecell>

down_2["frac"] = down_2["c"] / float(len(cp73_chronicHigh_down.index))

# <codecell>

up_2["frac"] = up_2["c"] / float(len(cp73_chronicHigh_up.index))

# <codecell>

both_3 = up_3.merge(down_3, left_on="pair", right_on="pair", how="outer", suffixes=("_up","_down")).fillna(0)

# <codecell>

both_2 = up_2.merge(down_2, left_on="pair", right_on="pair", how="outer", suffixes=("_up","_down")).fillna(0)

# <codecell>

both_3["diff"] = both_3["frac_up"] - both_3["frac_down"] 

# <codecell>

both_2["diff"] = both_2["frac_up"] - both_2["frac_down"] 

# <codecell>

both_2.sort_index(by="diff", ascending=False)[0:20] 

# <codecell>

both_3.sort_index(by="diff", ascending=True)[0:20] 

# <codecell>

pd.DataFrame( pair_counts.items(), columns=["pair", "c"] ).sort_index(by="c", ascending=False)[0:20]

# <codecell>

pair_counts.items()

# <codecell>

len( list( itertools.combinations( cz.tf.unique(), 2) ) )

# <codecell>

cz.groupby("tf").count().sort_index(by="chea_id", ascending=False)[0:20]

# <codecell>

alib.plots.clusterHeatmap( czm_df, "", None, None, 
                            cluster_columns=True, 
                            cluster_rows=True, 
                            distmethod='jaccard')

# <codecell>

import pandas as pd

# <codecell>

tfppi = pd.DataFrame.from_csv("/data/adrian/Dropbox/Data/mouse-2hybrid/mouse2hybid-table3.tab", sep="\t") 

# <codecell>

tfppi

# <codecell>

top_tfs[k].tf_name

# <codecell>

tfppi[ ["geneid1symbol", "geneid2symbol"] ]

# <codecell>

tfppi["tf1"] = tfppi["geneid1symbol"].apply(lambda x: x.upper())

# <codecell>

tfppi["tf2"] = tfppi["geneid2symbol"].apply(lambda x: x.upper())

# <codecell>

tfppi[["tf1", "tf2"]].select( lambda x: tfppi.ix[x, "tf1"] in set( top_tfs[k].tf_name ) and tfppi.ix[x, "tf2"] in set( top_tfs[k].tf_name )  ).to_json(orient='records')

# <codecell>

tfppi

# <codecell>

tf_tf = tfppi[["tf1", "tf2"]]

# <codecell>

set( top_tfs[k].tf_name )

# <codecell>

tf_tf.columns = ["source_id", "target_id"]

# <codecell>

cPickle.dump( tf_tf, open("/data/adrian/c_tf_tf.pickle", "w"))

# <codecell>

len( top_tfs[k].index )

# <codecell>

tfc = zeros( (202, 202) )

# <codecell>

invindex = dict( [ (a[1], a[0]) for a in enumerate( top_tfs[k].index )] ) 

# <codecell>

for i in tfppi.index:
    ppi = tfppi.ix[i, :]
    try:
        tfc[invindex[ppi["tf1"]], invindex[ppi["tf2"]] ] = 1
        tfc[invindex[ppi["tf2"]], invindex[ppi["tf2"]] ] = 1
    except:
        pass

# <codecell>

import scipy.cluster

# <codecell>

distances = scipy.cluster.hierarchy.distance.pdist(tfc + 0.001, "jaccard") 

# <codecell>

rowY = fastcluster.linkage( distances )
rowZ = scipy.cluster.hierarchy.dendrogram(rowY, orientation='right', no_plot=False)

# <codecell>

top_tfs[k].ix[ top_tfs[k].index[ rowZ['leaves'] ], : ].drop("tf_associations", axis=1).to_json(orient='records')

# <codecell>

log2(1.5)

# <codecell>

tf_target_counts.merge( loadCheaTable.o.chea.val, left_on="tf", right_on="target")[0:10]

# <codecell>

loadCheaTable.o.chea.val.tf.unique()

# <codecell>

tf_target_counts

# <codecell>

import cPickle

# <codecell>

cPickle.dump( motifs, open("/data/adrian/c_motifs.pickle","w"), protocol=2)

# <codecell>

tf

# <codecell>

cPickle.dump( top_tfs, open("/data/adrian/c_toptfs.pickle","w"), protocol=2)

# <codecell>

cPickle.dump( mm9_motif_gene_counts, open("/data/adrian/c_mm9_motif_gene_counts.pickle", "w"), protocol=2)

# <codecell>

cPickle.dump( tf_target_counts, open("/data/adrian/c_tf_target_counts", "w"), protocol=2)

# <codecell>

cPickle.dump( motif_tf["mat_tf_graph"], open("/data/adrian/c_motif_tf_graph.pickle", "w"))

# <codecell>

import json

# <codecell>

motif_keys = motifs.keys()

# <codecell>

motif_keys.sort()

# <codecell>

motif_keys

# <codecell>

json.dumps( map(dict, [zip(["name", "key"] , (k, hashlib.new('sha1', k).hexdigest())) for k in motif_keys]) )

# <codecell>

motifs[k].motif_occurrences[0]

# <codecell>


# <codecell>

group_a = cPickle.load( open("/data/adrian/group_a.pickle") )
group_b = cPickle.load( open("/data/adrian/group_b.pickle") )
group_c = cPickle.load( open("/data/adrian/group_c.pickle") )
group_d = cPickle.load( open("/data/adrian/group_d.pickle") )

# <codecell>

chea = loadCheaTable.o.chea.val

# <codecell>

chea.

# <codecell>

import pd_analysis as pda

# <codecell>

specific_groups_enrichment = {}

# <codecell>

specific_groups_motifs = {}

# <codecell>

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

# <codecell>

mm9_gene_motifs = loadGeneMotifs.o.mm9_gene_motifs.val

# <codecell>

group_d.symbol.unique()

# <codecell>

group_d

# <codecell>

getMotifMatches( group_d ).sort_index(by="count", ascending=False).drop("motif_occurrences",axis=1)[0:10]

# <codecell>

for (group_name, g) in [("group A; table 22", group_a),  
                        ("group B; table 23", group_b),  
                        ("group C; table 24", group_c),
                        ("group D; table 25", group_d)]:
    specific_groups_enrichment[group_name] = getTFMatches(g).sort_index( by="count", ascending=False )

# <codecell>

for (group_name, g) in [("group A; table 22", group_a),  
                        ("group B; table 23", group_b),  
                        ("group C; table 24", group_c),
                        ("group D; table 25", group_d)]:
    a = specific_groups_enrichment[group_name].tf_associations
    specific_groups_enrichment[group_name]["targets"] = [str( list( a[x].symbol.unique() ) ) for x in a.index]
    specific_groups_enrichment[group_name] = ( specific_groups_enrichment[group_name]
      .merge( tf_target_counts, left_on="tf_name", right_on="tf")
    )

# <codecell>

cp_73_chronic_high_vs_low = cPickle.load( open("/data/adrian/aim_correlated_table13.pickle") )

# <codecell>

cn_cp73_chronicHigh_vs_Low_foldchange = pda.tm.e([  ("cmp", ["6-OHDA, chronicHigh", "6-OHDA, chronicLow"]),
            ("ct", "cp73"),
            ("st", "fc_medians"),
        ])

# <codecell>

t13_motifs = {}
t13_tfs = {}

t13_up = cp_73_chronic_high_vs_low.select( lambda x: cp_73_chronic_high_vs_low.ix[x, cn_cp73_chronicHigh_vs_Low_foldchange] > 0 )
t13_down = cp_73_chronic_high_vs_low.select( lambda x: cp_73_chronic_high_vs_low.ix[x, cn_cp73_chronicHigh_vs_Low_foldchange] < 0 )

t13_motifs[ "aim_corr_table13 - UP"] = getMotifMatches(t13_up)
t13_motifs[ "aim_corr_table13 - DOWN"] = getMotifMatches(t13_down)
t13_motifs[ "aim_corr_table13 - ALL"] = getMotifMatches(cp_73_chronic_high_vs_low)


t13_tfs[ "aim_corr_table13 - UP"] = getTFMatches(t13_up) 
t13_tfs[ "aim_corr_table13 - DOWN"] = getTFMatches(t13_down) 
t13_tfs[ "aim_corr_table13 - ALL"] = getTFMatches(cp_73_chronic_high_vs_low) 

# <codecell>

t13_tfs

# <codecell>

for k in t13_tfs.keys():
    a = t13_tfs[k].tf_associations
    t13_tfs[k]["targets"] = [str( list( a[x].symbol.unique() ) ) for x in a.index]
    t13_tfs[k] = t13_tfs[k].merge( tf_target_counts, left_on="tf_name", right_on="tf")

# <codecell>

for k in t13_motifs.keys():    
    a = t13_motifs[k].motif_occurrences
    t13_motifs[k]["targets"] = [str( list( a[x].symbol.unique() ) ) for x in a.index]
    t13_motifs[k] = t13_motifs[k].merge( mm9_motif_gene_counts, left_on="motif_name", right_on="motif_name")

# <codecell>

mm9_motif_gene_counts

# <codecell>

for (group_name, g) in [("group A; table 22", group_a),  
                        ("group B; table 23", group_b),  
                        ("group C; table 24", group_c),
                        ("group D; table 25", group_d)]:
    specific_groups_motifs[group_name] = getMotifMatches(g).sort_index( by="count", ascending=False )

# <codecell>

for (group_name, g) in [("group A; table 22", group_a),  
                        ("group B; table 23", group_b),  
                        ("group C; table 24", group_c),
                        ("group D; table 25", group_d)]:
    a = specific_groups_motifs[group_name].motif_occurrences
    specific_groups_motifs[group_name]["targets"] = [str( list( a[x].symbol.unique() ) ) for x in a.index]
    specific_groups_motifs[group_name] = ( specific_groups_motifs[group_name]
      .merge( mm9_motif_gene_counts, left_on="motif_name", right_on="motif_name")
     )

# <codecell>

specific_groups_enrichment["group A; table 22"].merge( tf_target_counts, left_on="tf_name", right_index=True )

# <codecell>

addEnrichedTFStats.i.top_tfs.val = t13_tfs
addEnrichedTFStats.run()

# <codecell>

addEnrichedMotifStats.i.mm9_motif_gene_counts.val

# <codecell>

addEnrichedMotifStats.i.top_motifs.val = t13_motifs
addEnrichedMotifStats.run()

# <codecell>

e = pd.ExcelWriter("/data/adrian/specific_groups_transcripReg_09122013.xls")

# <codecell>

b = a[0]

# <codecell>

pd_tf_enrich.keys()

# <codecell>

pd_tf_enrich['[["cmp", ["6-OHDA, chronicSaline", "Ascorbate, chronicSaline"]], ["ct", "cp73"], ["dir", "down"]]']

# <codecell>

e.save()

# <codecell>

pd_motif_enrich.keys()

# <codecell>

pd_motif_enrich[k]

# <codecell>

e = pd.ExcelWriter("/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/2013_10_11/supptables/pd_heiman2013_supptable_01B.xls")
for (n, k) in [("all", '[["cmp", ["6-OHDA, chronicSaline", "Ascorbate, chronicSaline"]], ["ct", "cp73"], ["dir", "any"]]'),
          ("up", '[["cmp", ["6-OHDA, chronicSaline", "Ascorbate, chronicSaline"]], ["ct", "cp73"], ["dir", "up"]]'),
          ("down", '[["cmp", ["6-OHDA, chronicSaline", "Ascorbate, chronicSaline"]], ["ct", "cp73"], ["dir", "down"]]')]:
#    print k
#    specific_groups_enrichment[k].drop("tf_associations", axis=1).sort_index(by="pval").to_excel( e, sheet_name=k + "; chea")
    pd_motif_enrich[k].drop("motif_occurrences", axis=1).sort_index(by="pval").to_excel( e, sheet_name=n +"; motifs")
e.save()

# <codecell>

e = pd.ExcelWriter("/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/2013_10_11/supptables/pd_heiman2013_supptable_02B.xls")
for (n, k) in [("all", '[["cmp", ["6-OHDA, chronicSaline", "Ascorbate, chronicSaline"]], ["ct", "cp101"], ["dir", "any"]]'),
          ("up", '[["cmp", ["6-OHDA, chronicSaline", "Ascorbate, chronicSaline"]], ["ct", "cp101"], ["dir", "up"]]'),
          ("down", '[["cmp", ["6-OHDA, chronicSaline", "Ascorbate, chronicSaline"]], ["ct", "cp101"], ["dir", "down"]]')]:
#    print k
#    specific_groups_enrichment[k].drop("tf_associations", axis=1).sort_index(by="pval").to_excel( e, sheet_name=k + "; chea")
    pd_motif_enrich[k].drop("motif_occurrences", axis=1).sort_index(by="pval").to_excel( e, sheet_name=n +"; motifs")
e.save()

# <codecell>

for (n, k) in [("all", '[["cmp", ["6-OHDA, chronicSaline", "Ascorbate, chronicSaline"]], ["ct", "cp101"], ["dir", "any"]]'),
          ("up", '[["cmp", ["6-OHDA, chronicSaline", "Ascorbate, chronicSaline"]], ["ct", "cp101"], ["dir", "up"]]'),
          ("down", '[["cmp", ["6-OHDA, chronicSaline", "Ascorbate, chronicSaline"]], ["ct", "cp101"], ["dir", "down"]]')]:
#    print k
#    specific_groups_enrichment[k].drop("tf_associations", axis=1).sort_index(by="pval").to_excel( e, sheet_name=k + "; chea")
    (pd_motif_enrich[k]
     .drop("motif_occurrences", axis=1)
     .sort_index(by="pval")
     .to_csv("/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/2013_10_11/supptables/pd_heiman2013_supptable_02B_" + n + ".tab"))

# <codecell>

pd_motif_enrich[k]

# <codecell>

e = pd.ExcelWriter("/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/2013_10_11/supptables/pd_heiman2013_supptable_07C.xls")
for (n, k) in [("all", '[["cmp", ["6-OHDA, chronicHigh", "6-OHDA, chronicSaline"]], ["ct", "cp73"], ["dir", "any"]]'),
          ("up", '[["cmp", ["6-OHDA, chronicHigh", "6-OHDA, chronicSaline"]], ["ct", "cp73"], ["dir", "up"]]'),
          ("down", '[["cmp", ["6-OHDA, chronicHigh", "6-OHDA, chronicSaline"]], ["ct", "cp73"], ["dir", "down"]]')]:
#    print k
#    specific_groups_enrichment[k].drop("tf_associations", axis=1).sort_index(by="pval").to_excel( e, sheet_name=k + "; chea")
    pd_motif_enrich[k].drop("motif_occurrences", axis=1).drop("targets",axis=1).sort_index(by="pval").to_excel( e, sheet_name=n +"; motifs")
e.save()

# <codecell>

for (n, k) in [("all", '[["cmp", ["6-OHDA, chronicHigh", "6-OHDA, chronicSaline"]], ["ct", "cp73"], ["dir", "any"]]'),
          ("up", '[["cmp", ["6-OHDA, chronicHigh", "6-OHDA, chronicSaline"]], ["ct", "cp73"], ["dir", "up"]]'),
          ("down", '[["cmp", ["6-OHDA, chronicHigh", "6-OHDA, chronicSaline"]], ["ct", "cp73"], ["dir", "down"]]')]:
#    print k
#    specific_groups_enrichment[k].drop("tf_associations", axis=1).sort_index(by="pval").to_excel( e, sheet_name=k + "; chea")
    (pd_motif_enrich[k].drop("motif_occurrences", axis=1)
    .drop("targets",axis=1)
    .sort_index(by="pval")
    .to_csv("/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/2013_10_11/supptables/pd_heiman2013_supptable_07C_" + n + ".tab"))

# <codecell>

e = pd.ExcelWriter("/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/2013_10_11/supptables/pd_heiman2013_supptable_16.xls")
for (n, k) in [("all", '[["cmp", ["6-OHDA, chronicHigh", "6-OHDA, chronicLow"]], ["ct", "cp73"], ["dir", "any"]]'),
          ("up", '[["cmp", ["6-OHDA, chronicHigh", "6-OHDA, chronicLow"]], ["ct", "cp73"], ["dir", "up"]]'),
          ("down", '[["cmp", ["6-OHDA, chronicHigh", "6-OHDA, chronicLow"]], ["ct", "cp73"], ["dir", "down"]]')]:
#    print k
#    specific_groups_enrichment[k].drop("tf_associations", axis=1).sort_index(by="pval").to_excel( e, sheet_name=k + "; chea")
    pd_motif_enrich[k].drop("motif_occurrences", axis=1).drop("targets",axis=1).sort_index(by="pval").to_excel( e, sheet_name=n +"; motifs")
e.save()

# <codecell>

for (n, k) in [("all", '[["cmp", ["6-OHDA, chronicHigh", "6-OHDA, chronicLow"]], ["ct", "cp73"], ["dir", "any"]]'),
          ("up", '[["cmp", ["6-OHDA, chronicHigh", "6-OHDA, chronicLow"]], ["ct", "cp73"], ["dir", "up"]]'),
          ("down", '[["cmp", ["6-OHDA, chronicHigh", "6-OHDA, chronicLow"]], ["ct", "cp73"], ["dir", "down"]]')]:
#    print k
#    specific_groups_enrichment[k].drop("tf_associations", axis=1).sort_index(by="pval").to_excel( e, sheet_name=k + "; chea")
    (pd_motif_enrich[k]
     .drop("motif_occurrences", axis=1)
     .drop("targets",axis=1)
     .sort_index(by="pval")
     .to_csv("/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/2013_10_11/supptables/pd_heiman2013_supptable_16_" + n + ".tab"))

# <codecell>

specific_groups_motifs['group D; table 25'].drop("motif_occurrences", axis=1).sort_index(by="pval")[0:30]

# <codecell>

calcMotifGeneCounts

# <codecell>

mm9_motif_gene_counts

# <codecell>

tf_target_counts

# <codecell>

tf_target_counts

# <codecell>

a = specific_groups_enrichment['group D; table 25'].tf_associations

# <codecell>

specific_groups_enrichment['group D; table 25']["targets"] = [str( list( a[x].symbol.unique() ) ) for x in a.index]

# <codecell>

( specific_groups_enrichment['group D; table 25']
  .merge( tf_target_counts, left_on="tf_name", right_on="tf")
  .drop("tf_associations", axis=1)
  .sort_index(by="pval"))

# <codecell>

specific_groups_motifs['group C; table 24'].drop("motif_occurrences", axis=1).sort_index(by="pval")

# <codecell>

specific_groups_enrichment['group C; table 24'].drop("tf_associations", axis=1).sort_index(by="pval")[0:10]

# <codecell>

def getTFMatches( genegroup, combination_size=1):
    #genegroup["upper_sym"] = [a.upper() if isinstance(a, str) else "-" for a in genegroup.symbol ]
    genegroup["upper_sym"] = [str(a).upper() if isinstance(a, float) is not True else "-" for a in genegroup.symbol ]

    cz = chea.merge( genegroup, left_on="target", right_on="upper_sym")
    
    total_genes_in_group = len( pda.uniqueSym(genegroup.symbol) )
    
    print(len(cz.index))
    if len(cz.index) > 0:
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
    else:
        return None

# <codecell>

e = pd.ExcelWriter("/data/adrian/aim_correlated_transcriptReg.xls")
for k in t13_tfs.keys():
    t13_tfs[k].drop("tf_associations", axis=1).sort_index(by="pval").to_excel( e, sheet_name=k + "; chea")
    t13_motifs[k].drop("motif_occurrences", axis=1).sort_index(by="pval").to_excel( e, sheet_name=k +"; motifs")
e.save()

