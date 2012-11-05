import cherrypy
from cherrypy.process import wspbus, plugins
import mako
import sqlalchemy
from sqlalchemy.sql.expression import *

from mako.template import Template
from mako.lookup import TemplateLookup

def setup_routes():
    dispatcher = cherrypy.dispatch.RoutesDispatcher()
    return dispatcher



def start(config=None):
    conf = {
        '/' : {
        'request.dispatch' : setup_routes(),
        'tools.staticdir.root' : "/data/adrian/code/projects/broad/pdmouse/cf"
        },
        
        '/static' : { 'tools.staticdir.on' : True, 'tools.staticdir.dir' : 'static' }
        }
        
    if config:
        cherrypy.config.update(config)

#    SATool.SAEnginePlugin(cherrypy.engine).subscribe()
#    cherrypy.tools.db = SATool.SATool()
        
    cherrypy.server.socket_host = "0.0.0.0"
    cherrypy.server.socket_port = 9999
        
    app = cherrypy.tree.mount(None, config=conf)
    cherrypy.quickstart(app)

if __name__ == '__main__':
    start()
