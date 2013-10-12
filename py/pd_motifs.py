# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import pandas as pd
import urllib2
import IPython

# <codecell>

mm9_refFlat = pd.read_table("/data/adrian/Dropbox/Data/mouse_genome/refFlat.txt", 
                            names=["geneName", "name", "chrom", "strand", "txStart", 
                                     "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds"],
                            sep="\t") 

# <codecell>

mm9_refFlat[0:10]["txEnd"] - mm9_refFlat[0:10]["txStart"]

# <codecell>

IPython.display.HTML( mm9_refFlat[0:10].to_html() )

# <codecell>

swissreg_sites = pd.read_table("/data/adrian/Dropbox/Data/swissregulon/mm9_sites.gff", sep="\t", skiprows=1,
                               names=["chrom", "alg", "motif_type", 
                                      "start", "end", "score", "strand", "frame", "data"])

swissreg_sites["motif_name"] = swissreg_sites.apply( lambda x: urllib2.unquote( dict([a.split("=") for a in x["data"].split(";")])["Motif"] ), axis=1 )

# <codecell>

swissreg_sites_hg19 = pd.read_table("/data/adrian/Dropbox/Data/swissregulon/hg19_sites.gff", 
                                    sep="\t",
                                    skiprows=1,
                                    names=["chrom", "alg", "motif_type", 
                                      "start", "end", "score", "strand", "frame", "data"])

# <codecell>

swissreg_sites_hg19["motif_name"] = swissreg_sites_hg19.apply( lambda x: urllib2.unquote( dict([a.split("=") for a in x["data"].split(";")])["Motif"] ), axis=1 )

# <codecell>

r = dict([a.split("=") for a in swissreg_sites.data[0].split(";")])["Motif"]

# <codecell>

r["Motif"].decode("utf-8")

# <codecell>

urllib2.unquote(r["Motif"])

# <codecell>

import IPython

# <codecell>

IPython.display.HTML( swissreg_sites[0:10].to_html() )

# <codecell>

swissreg_sites.to_csv("/data/adrian/Dropbox/Data/swissregulon/mm9_sites.tab", sep="\t", index=False)

# <codecell>

swissreg_sites_hg19.to_csv("/data/adrian/Dropbox/Data/swissregulon/hg19_sites.tab", sep="\t", index=False)

# <codecell>

swissreg_sites[0:100].apply( lambda x: x["data"], axis=1) 

# <codecell>

mat_tf = []
for l in open("/data/adrian/Dropbox/Data/swissregulon/mat_TF_associations.mm").readlines():
    fields = l.split("\t")
    motif = fields[0]
    tfs = fields[2:]
    for tf in tfs:
        tfdata = tf.split(":")
        mat_tf.append( {"motif_name" : motif,
                        "tf_id_to_knownbind" : tfdata[0], 
                        "tf_entrez_gene_id" : tfdata[1], 
                        "tf_symbol" : tfdata[2],
                        "tf_description" : tfdata[3] } )
mat_tf = pd.DataFrame(mat_tf)

# <codecell>

IPython.display.HTML( mat_tf[0:10].to_html() )

# <codecell>

mat_tf.to_csv("/data/adrian/Dropbox/Data/swissregulon/mat_tf_graph.tab", sep="\t")

# <codecell>


