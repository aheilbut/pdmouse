# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 14:38:52 2013

@author: adrian
"""
import sys

import networkx as nx


class GZPort():
    def __init__(self, parent_node, name="", port_type=None, description=""):
        self.value = None
        self.parent_node=parent_node
        self.port_name=name
        self.port_type = port_type
        self.port_description = description

    @property    
    def val(self):
        return self.value
        
    @val.setter
    def val(self, val):
        self.value = val

class GZPortGroup():
    def __init__(self):
        pass
    
    def listports(self):
        return [getattr(self, a) for a in dir(self) if isinstance(getattr(self, a), GZPort)]

    def resultdict(self):
        return dict( [(a, getattr(self, a).val) for a in dir(self) if isinstance(getattr(self, a), GZPort)]  ) 


class GZWorkNode():
    pass

class GZConst(GZWorkNode):    
    def __init__(self, value):
        class PortsIn(GZPortGroup):
            pass
        self.i = PortsIn()
        
        class PortsOut(GZPortGroup):
            value = GZPort(self, "value")
        self.o = PortsOut()

        self.value = value
        
    def run(self):
        self.o.value.val = self.value

class GZWorkFlow():            
    def __init__(self):
        self.worknodes = {} # keyed by worknode class name to a list for each type
        self.worknodes_ids = {}
        self.wires = []
        self.wires_in = {}
        self.wires_out = {}
        self.worknode_types = {} # dict of types with current_id
        
    def addWorkNodes(self, worknodeList):
        for w in worknodeList:
            self.addWorkNode(w)
        
    def addWorkNode(self, worknode):
        # assign this node its id (just a sequence based on its name)
        worknode_id_type = worknode.__class__.__name__
        worknode_id_num = self.worknode_types.setdefault(worknode_id_type, 0) + 1
        self.worknode_types[worknode_id_type] = worknode_id_num

        self.worknodes.setdefault(worknode_id_type, {})[worknode_id_num] = worknode
        self.worknodes_ids[worknode] = (worknode_id_type, worknode_id_num)
        
    def connectPorts(self, source_node, source_port, sink_node, sink_port):
        source_port_id = self.worknodes_ids[source_node] + (source_port,)
        sink_port_id = self.worknodes_ids[sink_node] + (sink_port,)
        
        self.wires.append( (source_port_id, sink_port_id ) )
        self.wires_in.setdefault(sink_port_id, []).append( source_port_id )
        self.wires_out.setdefault(source_port_id, []).append( sink_port_id) 

    def plug(self, connections):
        for (source_port, sink_port) in connections:
            source_node = source_port.parent_node
            source_port_name = source_port.port_name
            
            sink_node = sink_port.parent_node
            sink_port_name = sink_port.port_name
            
            self.connectPorts(source_node, source_port_name, sink_node, sink_port_name)
        
    def verify():
        pass
    
    

class GZWorkRunner():
    def __init__(self, workflow):
        self.workflow = workflow
    
    def worknode_graph(self):
        g = nx.DiGraph()
        # construct graph of port connections
        # each worknode / port gets a node in the graph, with edges from wires
#        for worknode_type in self.workflow.worknodes.keys():
#            for (worknode_idnum, worknode) in self.workflow.worknodes[worknode_type].items():
#                g.add_vertices( [str((worknode_type, worknode_idnum, in_port)) for in_port in worknode.inputs] )
#                g.add_vertices( [str((worknode_type, worknode_idnum, out_port)) for out_port in worknode.outputs] )
#        g.add_edges( [(str(e[0]), str(e[1])) for e in self.workflow.wires] )
 
        # construct dependency graph between WorkNodes based on port connections
        # first, just add all the worknodes to the graph
        for worknode_type in self.workflow.worknodes.keys():
            for (worknode_idnum, worknode) in self.workflow.worknodes[worknode_type].items():
                g.add_node( (worknode_type, worknode_idnum), { "worknode" : worknode } )

        # add dependency edges
        for e in self.workflow.wires:
            if not g.has_edge( e[0][0:2], e[1][0:2] ):
                g.add_edge( e[0][0:2], e[1][0:2] )
       
        return g
        

    def runFlow(self, clean=False):
        # create new context to store the data
        context = GZRunContext(self.workflow)
        
        # get the ordering on worknodes to run from dependency graph
        g = self.worknode_graph()

        worknode_dict = dict(g.nodes(data=True))                
                
        # for each node in ordering, run it (using inputs from context) and add its outputs to the context
        for worknode_id in nx.topological_sort(g):
            print context.results.keys()
            print "running ", worknode_id
            current_worknode = worknode_dict[worknode_id]["worknode"]
            #args = []
            
            try:
                for in_port in current_worknode.i.listports():
                    # bind the values to this port
                    source_port_id = self.workflow.wires_in[worknode_id + (in_port.port_name,)][0]
                    print source_port_id
                    in_port.val = context.results[source_port_id[0:2]][source_port_id[2]]
#                    print self.workflow.wires_in[worknode_id + (port,)]
#                    args.extend( [ (port, context.results[source_port_id[0:2]][source_port_id[2]])
#                                for source_port_id in self.workflow.wires_in[worknode_id + (port,)] ])

            except:
                print sys.exc_info()
                return context
 #           print dict(args).keys()
#            print current_worknode.run(**dict(args))
#            context.results[worknode_id] =  current_worknode.run(**dict(args))
            current_worknode.run()
            context.results[worknode_id] = current_worknode.o.resultdict()
        return context             
        
        
class GZRunContext():
    def __init__(self, workflow):
        self.results = {}
        self.workflow = workflow

    
class GZResult():
    def __init__(self, name, timestamp, data, description=""):
        self.name = name
        self.timestamp = timestamp
        self.data = data
        self.description = description
        
    