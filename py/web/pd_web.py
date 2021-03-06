# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 13:13:24 2013

@author: adrian
"""

import cherrypy
from cherrypy.process import wspbus, plugins
import mako
import sqlalchemy
from sqlalchemy.sql.expression import *

from lib.sa.saplugin import SAEnginePlugin
from lib.sa.satool import SATool
from sqlalchemy import and_

import glob
import re
import sys

import urllib
import hashlib

from mako.template import Template
from mako.lookup import TemplateLookup

import pd_analysis as pda
import pd_locals
import pandas

import netfigdb

import numpy as np

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from cherrypy.lib import file_generator
import StringIO

import cPickle

import define_sample_subsets as ss

#import alib.wikipath


mo430symbol = pandas.read_table(pd_locals.datadir +  "/2012_10_29/mo4302symbols.tab")
mo430symbol.index = mo430symbol.probe_id

mo430names = pandas.read_table(pd_locals.datadir + "/2012_10_29/mo4302genenames.tab")
mo430names.index =  mo430names.probe_id

mo430info = mo430symbol.merge(mo430names)
mo430info.index = mo430info.probe_id

# load pd data
print "load pd data..."
pd_all = pandas.read_table(pd_locals.datadir + "/2012_10_29/PD_arraydata.tab")
pd_covar = pandas.read_table(pd_locals.datadir + "/2012_10_29/pd.covar.tab")

dim_descriptions = dict([map(str.rstrip, a.split("\t")) for a in open(pd_locals.datadir + "/2013_02_12/dim_description.txt").readlines()])


cp73_wide_info = cPickle.load(open(pd_locals.datadir + "/jan30/cp73_m.pandas.pickle"))
cp101_wide_info = cPickle.load(open(pd_locals.datadir + "/jan30/cp101_m.pandas.pickle"))

t2_up = cPickle.load(open(pd_locals.datadir + "/2013_02_12/t2_up.pickle"))
t2_down = cPickle.load(open(pd_locals.datadir + "/2013_02_12/t2_down.pickle"))

cp73_modelcomp = cPickle.load(open(pd_locals.datadir + "/2013_04_03/cp73_model_comp_stats.pickle"))
cp101_modelcomp = cPickle.load(open(pd_locals.datadir + "/2013_04_03/cp101_model_comp_stats.pickle"))

print "load motifs..."
motifs = cPickle.load(open(pd_locals.datadir + "/2013_sep_8/c_motifs.pickle"))
toptfs = cPickle.load(open(pd_locals.datadir + "/2013_sep_8/c_toptfs.pickle"))
motif_tf_graph = cPickle.load(open(pd_locals.datadir + "/2013_sep_8/c_motif_tf_graph.pickle"))
tf_tf_graph = cPickle.load(open(pd_locals.datadir + "/2013_sep_8/c_tf_tf.pickle"))

tl = TemplateLookup(directories=["templates"],
                    module_directory="tmp/mako_mod",
                    default_filters=["decode.utf8"],
                    output_encoding="utf-8")


wp = 0 #alib.wikipath.WikiPathSets()

class PDC():
    def __init__(self):
        pass


    def index(self):
        t = tl.get_template("index.html")
        return t.render(t2_up = t2_up, t2_down = t2_down)

    def probeset_report(self, probeset):

        cp73_modelComp = pda.compareModels(ss.ss_cp73_allchronic, pd_all, probeset)
        cp101_modelComp = pda.compareModels(ss.ss_cp101_allchronic, pd_all, probeset)

        try:
            (gene_symbol, gene_name) = mo430info.ix[ probeset, ["symbol", "gene_name"]]
        except:
            (gene_symbol, gene_name) = ("symbol_not_found", "name_not_found")
        
        # get figure
        aim_models_fig = 1
        expr_models_fig = 1


        t = tl.get_template("probeset_report.html")
        return t.render_unicode(probeset=probeset, gene_symbol=gene_symbol, gene_name=gene_name, cp73_modelComp = cp73_modelComp, cp101_modelComp = cp101_modelComp)

    def render_figure(self, probeset, fig_type, format, r=None):

        if fig_type == "scatter":
            probeset = urllib.unquote_plus(probeset)
            # generate figure 
            fig = Figure()
            fig.set_size_inches(16, 6)
            
            ax_cp73 = fig.add_subplot(121)
            ax_cp101 = fig.add_subplot(122)

            if r == "a": # auto range 
                xlim = None
            elif r == "f":  # full range
                xlim = (0, 15)
            
            pda.plotProbe(pd_all, pd_covar, probeset, "CP73", ax_cp73, legend=False, xlim=xlim)
            pda.plotProbe(pd_all, pd_covar, probeset, "CP101", ax_cp101, legend=True, xlim=xlim)

            fig.tight_layout(rect=[0, 0, 0.7, 1], h_pad=0.5)
            fig.set_facecolor('w')
            canvas = FigureCanvas(fig)
        # render

        if fig_type == "boxplot":
            probeset = urllib.unquote_plus(probeset)

            fig = Figure()
            fig.set_size_inches(15, 8)
            
            ax_cp73 = fig.add_subplot(121)
            ax_cp101 = fig.add_subplot(122)
            
            cp73_groups = [("ASC / chronic saline", ss.ss_cp73_ascorbate_saline.filenames),
                           ("OHDA / chronic saline", ss.cp73_chronic_saline.filenames),
                           ("OHDA / acute saline", ss.ss_cp73_acuteSaline.filenames),
                           ("OHDA / acute High L-DOPA", ss.ss_cp73_acuteHigh.filenames),
                           ("OHDA / chronic Low L-DOPA", ss.ss_cp73_chronicLow.filenames),
                           ("OHDA / chronic High L-DOPA", ss.ss_cp73_chronicHigh.filenames)
                           ]
            
            cp101_groups = [("ASC / chronic saline", ss.ss_cp101_ascorbate_saline.filenames),
                           ("OHDA / chronic saline", ss.cp101_chronic_saline.filenames),
                           ("OHDA / chronic Low L-DOPA", ss.ss_cp101_chronicLow.filenames),
                           ("OHDA / chronic High L-DOPA", ss.ss_cp101_chronicHigh.filenames)
                           ]
        
            pda.probe_boxPlots(probeset, pd_all, cp73_groups, "CP73", ax_cp73)
            pda.probe_boxPlots(probeset, pd_all, cp101_groups, "CP101", ax_cp101)
            
            fig.tight_layout(rect=[0, 0, 0.7, 1], h_pad=0.5)
            fig.set_facecolor('w')
            canvas = FigureCanvas(fig)
        
        if format == "png":
            cherrypy.response.headers['Content-Type'] = "image/png"
            buffer = StringIO.StringIO()
            canvas.print_png(buffer)
            buffer.seek(0)

            return file_generator(buffer)

        
    def overlap_table(self, cell_type, contrast, direction):
        t = tl.get_template("enrichment_table.html")

        if direction == "UP":
            test_set = t2_up[cell_type][contrast].symbol.unique()
        elif direction == "DOWN":
            test_set = t2_down[cell_type][contrast].symbol.unique()

        s = alib.wikipath.calc_overlap_stats(test_set, wp.mouse_wp_symbols, 20826)

        return t.render_unicode(overlaps = s)


    def result_table(self, result_set, format="HTML", sort_field=0, sort_dir="DESC", rows=100, dimlist=None, page=0, cell_type=None, contrast=None):
        if result_set == "cp73":
            s = cp73_wide_info
        elif result_set == "cp101":
            s = cp101_wide_info
        elif result_set == "t2_up":
            s = t2_up[cell_type][contrast]		
        elif result_set == "t2_down":
            s = t2_down[cell_type][contrast]      


        if sort_dir == "ASC":
            sort_asc = True
        else:
            sort_asc = False

        if sort_field is not None:
            s = s.sort([s.columns.values[int(sort_field)]], ascending=sort_asc)
        else:
            s = s


        if format == "txt":
            cherrypy.response.headers['Content-Type'] = "application/octet-stream"
            filename = "PDAIMS_2013_%s_%s_%s.txt" % (result_set, cell_type, contrast) 
            cherrypy.response.headers["Content-Disposition"] = "attachment; filename=%s" % filename
            buffer = StringIO.StringIO()
            s.to_csv(buffer, sep="\t")
            buffer.seek(0)

            return file_generator(buffer)

        s["probe_img"] = ["<img width=350 src='/probeset/figure/%s/scatter/png'/>" % urllib.quote_plus(probe_id) for probe_id in s.probe_id]

        s["probe_id"] = ["<a target='_detail'  href='/probeset/%s'>%s</a>" % (urllib.quote_plus(probe_id), probe_id) for probe_id in s.probe_id]
        s["symbol"] = ["<a target='_ncbi' href='http://www.ncbi.nlm.nih.gov/gene/?term=%s'>%s</a>" % (sym, sym) for sym in s.symbol]


#        c = [s.columns[-1]]
#        c.extend(s.columns[0:-1])
#        s = s[c]
        row_start = int(page) * 100
        #print 

        t = tl.get_template("result_table.html")

        if dimlist is None:
            dimarray = range(0, len(s.columns)-1)
            dimlist = ",".join([str(q) for q in dimarray])
        else:
            dimarray = [int(i) for i in dimlist.split(",")]
            dimarray.sort()


        if dimarray is not None:
            column_list = s.columns[dimarray]
        else:
            column_list = s.columns

        return t.render_unicode( fulltable = s.ix[s.index[row_start:row_start + 100],:],
                                 displaytable = s.ix[s.index[row_start:row_start + 100], column_list],
                                 page=int(page), 
                                 col_widths={"gene_name" : "300px"},
                                 sort_field=int(sort_field),
                                 sort_dir=sort_dir,
                                 dimlist=dimlist,
                                 dimarray=dimarray,
                                 dim_descriptions=dim_descriptions,
                                 result_set=result_set,
                                 cell_type=cell_type,
                                 contrast=contrast)


    
    def gsea_list(self):
        t = tl.get_template("gsea_list.html")
        gsea_directories = glob.glob(pd_locals.gsea_dir + "*")
        gsea_directories.sort()
        
        return t.render_unicode(gsea_directories = gsea_directories)

    def dimlist(self):
        t = tl.get_template("dim_list.html")
        
        return t.render_unicode(dim_descriptions=dim_descriptions)


    
    def outliers(self, cptype, mouseid, foldmin=2, zmin=3):
        t = tl.get_template("outliers.html")
        
        outliers = pda.outlierSelect(pda.compareOutlier(cptype, int(mouseid)), foldmin, zmin)
        
        return t.render_unicode(cptype=cptype, mouseid=int(mouseid), outliers=outliers)
        
    
    def search(self, q):
        t= tl.get_template("search_results.html")
        
        def matchProbe(x, query):
            if ( re.match(".*%s.*" % query.lower(), mo430info.gene_name[x].lower()) 
                or re.match(".*%s.*" % query.lower(), mo430info.symbol[x].lower())
                or query == mo430info.probe_id[x]):
                return True
            else:
                return False        
            
        matches = mo430info.select(lambda x: matchProbe(x, q))
        
        return t.render_unicode(matches = matches.sort("gene_name"))
    
    
    def methods(self):
        t = tl.get_template("methods.html")
        return t.render_unicode()
    
    def aims(self):
        t = tl.get_template("aims.html")
        
        return t.render_unicode(pd_covar=pd_covar)
    
    def trm(self):
        t = tl.get_template("linked_tables.html")
        return t.render_unicode()
    

    def table_view(self, queryspec=None):
        """ return the html page with the tableviewer, hooked up to the requested data """
        t = tl.get_template("table_view.html")

        return t.render_unicode()
    
#    @cherrypy.tools.json_in()
    @cherrypy.tools.json_out()
    def cube_info(self, dataid):
#        try:
#            query = cherrypy.request.json
#        except:
#            print "no query specified"
            
        
        if dataid == "CP101":
            dataset = cp101_modelcomp
            title = "CP101 linear model comparison statistics"
            description = "nominal and bh adjusted p-values for F-tests comparing AIM ~ expression + C(dose) vs AIM ~ C(dose)"

        elif dataid == "CP73":
            dataset = cp101_modelcomp
            title = "CP73 linear model comparison statistics",
            description = "nominal and bh adjusted p-values for F-tests comparing AIM ~ expression + C(dose) vs AIM ~ C(dose)"
            
        elif dataid == "cp73_wide":
            dataset = cp73_wide_info
            dataset  = dataset.reindex(pandas.Series(dataset.index.values, name=dataset.index.name + "_index"))
            title = "cp73 wide"
            description = ""
                        
        recs = dataset.to_records()
        fieldset = [str(n) for n in recs.dtype.names]
        #columndef_dict = dict( [(f, {"id" : f, "name" : f, "field" : f, "sortable" : True  }) for f in fieldset] )

        columns = [ {"id" : f, "name" : f, "field" : f, "sortable" : True  } for f in fieldset] 
        
        return {                
            "columns" : columns,
            'success' : True,
            "title" : title,
            "description" : description,
            "total" : len(dataset)
        }




    @cherrypy.tools.json_in()
    @cherrypy.tools.json_out()
    def gzdata(self):
        try:
            query = cherrypy.request.json
            dataset = query["dataset"]
            print query
        except:
            print sys.exc_info()
            print "no query specified"
                    
        if dataset == "motifs":
            selection = query["selection"]        
     
            selectormap = dict([ (hashlib.new('sha1', k).hexdigest(), k) for k in motifs.keys()])
            key = selectormap[selection]
            motifdata = motifs[key].drop("motif_occurrences", axis=1).sort_index(by="pval")
            motifdata = motifdata.select(lambda x: motifdata.ix[x, "pval"] < 0.20)
            motifdata["divid"] = [ hashlib.new("sha1", m).hexdigest() for m in motifdata.motif_name ]
            result = map(lambda x: dict(zip(motifdata.columns, x)), motifdata.to_records(index=False))            
            columns = list( map(str, motifdata.columns))
            columns.remove("divid") 
            
            return { "success" : True, "data" : result, "columns" : columns, "idcol" : "divid"  }
            
        if dataset == "tfs":
            dataset = query["dataset"]
            selection = query["selection"]        
     
            selectormap = dict([ (hashlib.new('sha1', k).hexdigest(), k) for k in motifs.keys()])
            key = selectormap[selection]
            tfdata = toptfs[key].drop("tf_associations", axis=1).sort_index(by="pval")
            tfdata = tfdata.select(lambda x: tfdata.ix[x, "pval"] < 0.20)
            result = map(lambda x: dict(zip(tfdata.columns, x)), tfdata.to_records(index=False))            

            columns = map(str, tfdata.columns)            
            
            return { "success" : True, "data" : result, "columns" : columns, "idcol" : "tf_name"  }
            
            
        if dataset == "tf_motif_associations":
            dataset = query["dataset"]
            selection = query["selection"]             

            selectormap = dict([ (hashlib.new('sha1', k).hexdigest(), k) for k in motifs.keys()])
            key = selectormap[selection]
            motif_tf_graph["tf_u"] = [a.upper() for a in motif_tf_graph.tf_symbol];
    
            tfdata = toptfs[key]
            tfdata = tfdata.select(lambda x: tfdata.ix[x, "pval"] < 0.20)
            
            mj = motifs[key].merge( motif_tf_graph, left_on="motif_name", right_on="motif_name")
            mj["divid"] = [ hashlib.new("sha1", m).hexdigest() for m in mj.motif_name ]
            mj = mj.select(lambda x: mj.ix[x, "pval"] < 0.20)            
            
            
            corroborated = tfdata.merge( mj, left_on="tf_name", right_on="tf_u", suffixes=("_tf", "_motif"))
            
            
            result = corroborated[["tf_name", "divid"]]
            result.columns = ["source_id", "target_id"]
            result = map(lambda x: dict(zip(result.columns, x)), result.to_records(index=False))
            return { "success" : True, 
                    "connections" : result, 
                    "connector_class" : "connected_across",
                    "connector_color" : "rgba(243,140,18,0.5)", 
                    "connector_locations" : ["Right", "Left"] } 


        if dataset == "tf_tf_physical":
            dataset = query["dataset"]
            selection = query["selection"]        
         
            result = map(lambda x: dict(zip(tf_tf_graph.columns, x)), tf_tf_graph.to_records(index=False))
            return { "success" : True, 
                     "connections" : result,
                     "connector_class" : "connected_same",
                     "connector_color" : "rgba(243, 40, 18, 0.5)",
                     "connector_locations" : ["Left", "Left"]                    
                    }            
            
            
        if dataset == "tf_targets":
            selectormap = dict([ (hashlib.new('sha1', k).hexdigest(), k) for k in motifs.keys()])
            key = selectormap[query["group"]]

            tfdata = toptfs[key]
            targets = tfdata.ix[query["tf"], "tf_associations"]

            
            result = map(lambda x: dict(zip(targets.columns, x)), targets.to_records(index=False))
            
            return { "success" : True,
                    "data" : result, 
                    "columns" : map(str, targets.columns),
                    "idcol" : "chea_id" 
                    }
            
        if dataset == "motif_targets":
            selectormap = dict([ (hashlib.new('sha1', k).hexdigest(), k) for k in motifs.keys()])
            key = selectormap[query["group"]]
            motif_key = query["motif_key"]

            motifdata = motifs[key].sort_index(by="pval")
            motifdata["divid"] = [ hashlib.new("sha1", m).hexdigest() for m in motifdata.motif_name ]
            motif_occurrences = motifdata.select(lambda x: motifdata.ix[x, "divid"] == motif_key)["motif_occurrences"][0]
 #           motif_occurrences["id"] = [ hashlib.new("sha1", m).hexdigest() for m in motif_occurrences.index ]
            
            result = map(lambda x: dict(zip(motif_occurrences.columns, x)), motif_occurrences.to_records(index=False))            
            
            return {
                "success" : True,
                "data" : result,
                "columns" : map(str, motif_occurrences.columns),
                "idcol" : "id"
            }
    
    
    @cherrypy.tools.json_in()
    @cherrypy.tools.json_out()
    def cube_select(self):
        """ return json formatted slickgrid table data for the specified query """  
        try: 
            query = cherrypy.request.json

        except:
            print "no query specified" 
        
        def cleanFormat(z):
            if type(z) is not str and np.isnan(z):
                return ""
            elif type(z) is float:
                return round(z,3)
            else:
                return z

 
        start = int(query.get("offset", 0))
        finish = start + int(query.get("count", 10))

        print query["cptype"]
            
        if query["cptype"] == "CP101":
            dataset = cp101_modelcomp
            recs = dataset[start:finish].to_records()
            datadesc = "CP101 linear model comparison statistics"
        elif query["cptype"] == "CP73":
            dataset = cp73_modelcomp
            recs = dataset[start:finish].to_records()
            datadesc = "CP73 linear model comparison statistics"
        elif query["cptype"] == "cp73_wide":
            dataset = cp73_wide_info
            dataset  = dataset.reindex(pandas.Series(dataset.index.values, name=dataset.index.name + "_index"))
            datadesc = "CP73 wide"
            recs = dataset[start:finish].to_records()
            
            
        fieldset = [str(n) for n in recs.dtype.names]    
        reclist = [dict( zip( fieldset, 
                              [cleanFormat(z) for z in row]) )
                   for row in recs.tolist()]
        
        #        columndef_dict = dict( [(f, {"id" : f, "name" : f, "field" : f, "sortable" : True }) for f in fieldset] )
        
        columns = [ {"id" : f, "name" : f, "field" : f, "sortable" : True } for f in fieldset]

#        columndef_dict["probe_id"]["formatter"] = probeset_formatter
        
#        columndef = columndef_dict.values()
       
        print "start: ", start
        print "finish: ", finish
        
        result = { "success" : True,  
                   "data" : reclist, 
                   "data_description" : datadesc, 
                   "start": start, 
                   "total": len(dataset), 
                   "count": finish-start, 
                   "columns" : columns 
                   }  
            
        return result
        
    
class NetView():
    def __init__(self):
        pass
        
    def netView_app(self):
        t = tl.get_template("netfig.template.html")
        
        return t.render_unicode()

    @cherrypy.tools.json_out()
    def netView_netfig(self, netfig_id):
        s = cherrypy.request.db

        network_fig = s.query( netfigdb.Netfig ).filter( netfigdb.Netfig.netfig_id == int(netfig_id)).one()

        return network_fig.toDict()

    @cherrypy.tools.json_out()
    def netView_expressionchanges(self, netfig_id):
        s = cherrypy.request.db
        fig = s.query(netfigdb.Netfig).filter(netfigdb.Netfig.netfig_id == int(netfig_id)).one()
        gene_symbols = [n.node_obj_id for n in fig.nodes if n.node_obj_idtype == "entrez_gene_symbol"]
        expression_changes = dict([ (c.gene_symbol, c.to_dict())
                                    for c in
                                    s.query(netfigdb.NetfigMeanFoldChanges)
                                    .filter(netfigdb.NetfigMeanFoldChanges.gene_symbol.in_(gene_symbols) )
                                    .all()])
        return expression_changes

    @cherrypy.tools.json_out()
    def netView_netfigList(self):
        s = cherrypy.request.db
        figlist = [{ "netfig_id" :f.netfig_id, "creator" : f.creator_id, "title" : f.title, "description" : f.description} for f in s.query(netfigdb.Netfig).all()]
        return figlist

    @cherrypy.tools.json_in()
    @cherrypy.tools.json_out()
    def netView_update_node(self, netfig_node_id):
        try:
            s = cherrypy.request.db
            d = cherrypy.request.json
            print d
            netfig_node_id = int(netfig_node_id)
            node = s.query(netfigdb.NetfigNode).filter(netfigdb.NetfigNode.netfig_node_id == netfig_node_id).one()

            for k in d.keys():
                if k is not None:
                    setattr(node, k, d[k])
            print node.node_pos_top
            print node.node_pos_left
            s.commit()
            return node.toDict()

        except:
            s.rollback()
            print sys.exc_info()
            return "error"

pdc = PDC()
netview = NetView()

def setup_routes():
    dispatch = cherrypy.dispatch.RoutesDispatcher()

    dispatch.connect(name = "probeset_report",
                     route = "/probeset/{probeset}",
                     controller = pdc,
                     action = "probeset_report")

    dispatch.connect(name = "render_figure",
                     route = "/probeset/figure/{probeset}/{fig_type}/{r}/{format}",
                     controller = pdc,
                     action = "render_figure" )

    dispatch.connect(name = "result_table",
                     route = "/resulttable/{result_set}",
                     controller = pdc, 
                     action = "result_table")

    dispatch.connect(name = "overlaps",
                     route = "/overlaps",
                     controller = pdc, 
                     action = "overlap_table")

    dispatch.connect(name = "gsea_list",
                     route = "/gsea",
                     controller = pdc, 
                     action = "gsea_list")
    
    dispatch.connect(name = "dim_list",
                     route= "/dimlist",
                     controller = pdc,
                     action = "dimlist")
    
    dispatch.connect(name = "search",
                     route = "/search",
                     controller = pdc,
                     action = "search")
    
    dispatch.connect(name = 'methods',
                     route = '/methods',
                     controller = pdc,
                     action = 'methods')
    
    dispatch.connect(name = "outliers",
                     route = "/outliers",
                     controller = pdc, 
                     action = "outliers")
    
    dispatch.connect(name = "index",
                     route = "/",
                     controller = pdc,
                     action = "index")

    dispatch.connect(name = "aims",
                     route = "/aims",
                     controller = pdc,
                     action = "aims")
    dispatch.connect(name = "tableview",
                     route = "/table",
                     controller = pdc,
                     action = "table_view")

    dispatch.connect(name = "cube_select",
                     route = "/cube_select", 
                     controller = pdc,
                     action = 'cube_select')

    dispatch.connect(name = "cube_info",
                     route = "/cube_info/{dataid}",
                     controller = pdc,
                     action = "cube_info")
    
    dispatch.connect(name = "gzdata",
                     route = "/gzdata",
                     controller = pdc,
                     action = "gzdata")
                     
    dispatch.connect(name = "trm",
                     route = "/trm",
                     controller = pdc, 
                     action = "trm")
    

    dispatch.connect(name = "network_figure_list",
                     route = "/network/figures",
                     controller = netview,
                     action = "netView_netfigList")

    dispatch.connect(name = "network", 
                     route = "/network",
                     controller = netview,
                     action = "netView_app")

    dispatch.connect(name = "network_expression",
                     route = "/network/expression/{netfig_id}",
                     controller = netview,
                     action = "netView_expressionchanges")

    dispatch.connect(name="network_figure",
                     route="/network/figure/{netfig_id}",
                     controller=netview,
                     action="netView_netfig")

    dispatch.connect(name="network_node_update",
                     route="/network/nodes/{netfig_node_id}",
                     controller=netview,
                     action="netView_update_node")

    return dispatch





def start(config=None):
    conf = {
        '/' : {
            'request.dispatch' : setup_routes(),
            'tools.staticdir.root' : pd_locals.staticdir,
            'tools.gzip.on' : True,
            'tools.db.on' : True
            },

        '/static' : { 'tools.staticdir.on' : True, 
                      'tools.staticdir.dir' : 'static' }
    }

    if config:
        cherrypy.config.update(config)

#    SATool.SAEnginePlugin(cherrypy.engine).subscribe()
#    cherrypy.tools.db = SATool.SATool()

    SAEnginePlugin(cherrypy.engine, pd_locals.dbstring).subscribe()
    cherrypy.tools.db = SATool()

    cherrypy.server.socket_host = "0.0.0.0"
    cherrypy.server.socket_port = pd_locals.web_port 

    app = cherrypy.tree.mount(None, config=conf)
    cherrypy.quickstart(app)

if __name__ == '__main__':
    start()
