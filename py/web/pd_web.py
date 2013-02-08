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

tl = TemplateLookup(directories=["templates"],
                    module_directory="tmp/mako_mod",
                    default_filters=["decode.utf8"],
                    output_encoding="utf-8")

# load pd data
pd_all = pandas.read_table(pd_locals.datadir + "/Oct29/PD_arraydata.tab")
pd_covar = pandas.read_table(pd_locals.datadir + "/Oct29/pd.covar.tab")


cp73_wide_info = cPickle.load(open(pd_locals.datadir + "/jan30/cp73_m.pandas.pickle"))
cp101_wide_info = cPickle.load(open(pd_locals.datadir + "/jan30/cp101_m.pandas.pickle"))

t2_up = cPickle.load(open(pd_locals.datadir + "/jan30/t2_up.pickle"))
t2_down = cPickle.load(open(pd_locals.datadir + "/jan30/t2_down.pickle"))

wp = alib.wikipath.WikiPathSets()

class PDC():
    def __init__(self):
        pass
    

    def index(self):
        t = tl.get_template("index.html")
        return t.render(t2_up = t2_up, t2_down = t2_down)

    def probeset_report(self, probeset):
        t = tl.get_template("probeset_report.html")
        return t.render()

    def render_figure(self, probeset, fig_type, format):

        if fig_type == "scatter":
            probeset = urllib.unquote_plus(probeset)
            # generate figure 
            fig = Figure()
            fig.set_size_inches(15, 6)
            ax_cp73 = fig.add_subplot(121)
            ax_cp101 = fig.add_subplot(122)

            pda.plotProbe(pd_all, pd_covar, probeset, "CP73", ax_cp73)
            pda.plotProbe(pd_all, pd_covar, probeset, "CP101", ax_cp101)

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


    def result_table(self, result_set, sort_field=0, sort_dir="DESC", rows=100, dimlist=None, page=0, cell_type=None, contrast=None):
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

        s["probe_img"] = ["<img width=350 src='/probeset/figure/%s/scatter/png'/>" % urllib.quote_plus(probe_id) for probe_id in s.probe_id]

        s["probe_id"] = ["<a href='/probeset/figure/%s/scatter/png'>%s</a>" % (urllib.quote_plus(probe_id), probe_id) for probe_id in s.probe_id]
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
                                 result_set=result_set,
                                 cell_type=cell_type,
                                 contrast=contrast)



pdc = PDC()

def setup_routes():
    dispatch = cherrypy.dispatch.RoutesDispatcher()

    dispatch.connect(name = "probeset_report",
                     route = "/probeset/{probeset}",
                     controller = pdc,
                     action = "probeset_report")

    dispatch.connect(name = "render_figure",
                     route = "/probeset/figure/{probeset}/{fig_type}/{format}",
                     controller = pdc,
                     action = "render_figure")

    dispatch.connect(name = "result_table",
                     route = "/resulttable/{result_set}",
                     controller = pdc, 
                     action = "result_table")

    dispatch.connect(name = "overlaps",
                     route = "/overlaps",
                     controller = pdc, 
                     action = "overlap_table")



    dispatch.connect(name = "index",
                     route = "/",
                     controller = pdc,
                     action = "index")


    return dispatch





def start(config=None):
    conf = {
        '/' : {
        'request.dispatch' : setup_routes(),
        'tools.staticdir.root' : "/data/adrian/code/projects/broad/pdmouse/py/web"
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
