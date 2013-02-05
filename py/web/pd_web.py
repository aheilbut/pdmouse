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

from mako.template import Template
from mako.lookup import TemplateLookup

import pd_analysis as pda
import pandas

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from cherrypy.lib import file_generator
import StringIO

import cPickle

tl = TemplateLookup(directories=["templates"],
                    module_directory="tmp/mako_mod",
                    default_filters=["decode.utf8"],
                    output_encoding="utf-8")

# load pd data
pd_all = pandas.read_table("/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/Oct29/PD_arraydata.tab")
pd_covar = pandas.read_table("/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/Oct29/pd.covar.tab")


#cp73_wide_info = cPickle.load(open("/data/adrian/data/temp/cp73_wide.pandas.pickle"))
#cp101_wide_info = cPickle.load(open("/data/adrian/data/temp/cp101_wide.pandas.pickle"))

t2_up = cPickle.load(open("/home/adrian/Dropbox/Projects/Broad/PD_mouse/results/jan30/t2_up.pickle"))

class PDC():
    def __init__(self):
        pass
    
    def probeset_report(self, probeset):
        t = tl.get_template("probeset_report.html")
        return t.render()

    def render_figure(self, probeset, fig_type, format):

        if fig_type == "scatter":
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

    def geneset_enrichment_table(self, cell_type, contrast):
        pass


    def result_table(self, result_set, sort_field=0, sort_dir="DESC", rows=100, page=0):
        if result_set == "cp73":
            s = cp73_wide_info
        elif result_set == "cp101":
            s = cp101_wide_info
    	elif result_set == "t2_up_cp73_low_vs_saline":
            s = t2_up["CP73"]["low vs. saline"]		
        elif result_set == "t2_up_cp73_high_vs_saline":
            s = t2_up["CP73"]["high vs. saline"]     
        elif result_set == "t2_up_cp73_low_vs_high":
            s = t2_up["CP73"]["low vs. high"]     


        if sort_dir == "ASC":
            sort_asc = True
        else:
            sort_asc = False

        if sort_field is not None:
            s = s.sort([s.columns.values[int(sort_field)]], ascending=sort_asc)
        else:
            s = s

        s["probe_id"] = ["<a href='/probeset/figure/%s/scatter/png'>%s</a>" % (probe_id, probe_id) for probe_id in s.probe_id]
        s["symbol"] = ["<a target='_ncbi' href='http://www.ncbi.nlm.nih.gov/gene/?term=%s'>%s</a>" % (sym, sym) for sym in s.symbol]


#        c = [s.columns[-1]]
#        c.extend(s.columns[0:-1])
#        s = s[c]
        row_start = int(page) * 100
        #print 

        t = tl.get_template("result_table.html")

        return t.render_unicode(table = s.ix[s.index[row_start:row_start + 100],:],
                                 page=int(page), 
                                 col_widths={"gene_name" : "300px"},
                                 sort_field=int(sort_field),
                                 sort_dir=sort_dir,
                                 result_set=result_set)



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


    return dispatch





def start(config=None):
    conf = {
        '/' : {
        'request.dispatch' : setup_routes(),
        'tools.staticdir.root' : "/data/adrian/code/projects/broad/pdmouse/web"
        },
        
        '/static' : { 'tools.staticdir.on' : True, 'tools.staticdir.dir' : 'static' }
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
