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

import glob
import re


import urllib

from mako.template import Template
from mako.lookup import TemplateLookup

import pd_analysis as pda
import pd_locals
import pandas

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from cherrypy.lib import file_generator
import StringIO

import cPickle

import alib.wikipath


mo430symbol = pandas.read_table(pd_locals.datadir +  "/Oct29/mo4302symbols.tab")
mo430symbol.index = mo430symbol.probe_id

mo430names = pandas.read_table(pd_locals.datadir + "/Oct29/mo4302genenames.tab")
mo430names.index =  mo430names.probe_id

mo430info = mo430symbol.merge(mo430names)
mo430info.index = mo430info.probe_id

# load pd data
pd_all = pandas.read_table(pd_locals.datadir + "/Oct29/PD_arraydata.tab")
pd_covar = pandas.read_table(pd_locals.datadir + "/Oct29/pd.covar.tab")

dim_descriptions = dict([map(str.rstrip, a.split("\t")) for a in open(pd_locals.datadir + "/feb25/dim_description.txt").readlines()])


cp73_wide_info = cPickle.load(open(pd_locals.datadir + "/jan30/cp73_m.pandas.pickle"))
cp101_wide_info = cPickle.load(open(pd_locals.datadir + "/jan30/cp101_m.pandas.pickle"))

t2_up = cPickle.load(open(pd_locals.datadir + "/feb25/t2_up.pickle"))
t2_down = cPickle.load(open(pd_locals.datadir + "/feb25/t2_down.pickle"))


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


cp73_chronic_saline = pd_covar.select(lambda x: pd_covar.ix[x, "MouseType"] == "CP73" and pd_covar.ix[x, "LesionType"] == "6-OHDA" and pd_covar.ix[x, "DrugTreat"] == "Chronic saline")
cp101_chronic_saline = pd_covar.select(lambda x: pd_covar.ix[x, "MouseType"] == "CP101" and pd_covar.ix[x, "LesionType"] == "6-OHDA" and pd_covar.ix[x, "DrugTreat"] == "Chronic saline")


tl = TemplateLookup(directories=["templates"],
                    module_directory="tmp/mako_mod",
                    default_filters=["decode.utf8"],
                    output_encoding="utf-8")


wp = alib.wikipath.WikiPathSets()

class PDC():
    def __init__(self):
        pass


    def index(self):
        t = tl.get_template("index.html")
        return t.render(t2_up = t2_up, t2_down = t2_down)

    def probeset_report(self, probeset):

        cp73_modelComp = pda.compareModels(ss_cp73_allchronic, pd_all, probeset)
        cp101_modelComp = pda.compareModels(ss_cp101_allchronic, pd_all, probeset)

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
            fig.set_size_inches(15, 6)
            
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
    
    
    def aims(self):
        t = tl.get_template("aims.html")
        
        return t.render_unicode()
    
    
pdc = PDC()

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
    
    

    return dispatch





def start(config=None):
    conf = {
        '/' : {
            'request.dispatch' : setup_routes(),
            'tools.staticdir.root' : pd_locals.staticdir
            },

        '/static' : { 'tools.staticdir.on' : True, 
                      'tools.staticdir.dir' : 'static' }
    }

    if config:
        cherrypy.config.update(config)

#    SATool.SAEnginePlugin(cherrypy.engine).subscribe()
#    cherrypy.tools.db = SATool.SATool()

    cherrypy.server.socket_host = "0.0.0.0"
    cherrypy.server.socket_port = 9080

    app = cherrypy.tree.mount(None, config=conf)
    cherrypy.quickstart(app)

if __name__ == '__main__':
    start()