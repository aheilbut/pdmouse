
from sqlalchemy import create_engine
from sqlalchemy import Table, Column, Integer, Float, String, Date, BigInteger, DateTime, Boolean, MetaData, ForeignKey, UniqueConstraint, ForeignKeyConstraint
from sqlalchemy.orm import mapper, sessionmaker, relationship, backref

import datetime

import pd_locals
import sys
metadata = MetaData()

engine = create_engine(pd_locals.dbstring)
#engine = create_engine('sqlite:///:memory:', echo=True)
Session = sessionmaker(engine)



from sqlalchemy.ext.declarative import declarative_base
Base = declarative_base()

class PDUser(Base):
    __tablename__ = 'pd_user'

    pd_user_id = Column('pd_user_id', Integer, primary_key=True)

    email = Column('email', String)
    password = Column("password", String)
    firstname = Column('firstname', String)
    lastname = Column("lastname", String)
    institution = Column("institution", String)
    create_time = Column("create_time", DateTime)
    last_login = Column("last_login", DateTime)


    def toDict(self):
        d = {
            "email" : self.email,
            "firstname" : self.firstname,
            "lastname" : self.lastname,
            "institution" : self.institution,
            "last_login" : self.last_login
        }
        return d

class Netfig(Base):
    __tablename__ = "netfig"

    netfig_id = Column("netfig_id", BigInteger, primary_key=True)
    creator_id = Column("creator", Integer, ForeignKey("pd_user.pd_user_id"))
    create_time = Column("create_time", DateTime)
    last_edit = Column("last_edit", DateTime)
    title = Column("title", String)
    description = Column("description", String)

    creator = relationship("PDUser", backref="figures")
#    def

    def __init__(self, creator_id, title, description):
        self.creator_id = creator_id,
        self.create_time = datetime.datetime.now()
        self.last_edit = datetime.datetime.now()
        self.title = title
        self.description = description

    def toDict(self):
        d = {
            "netfig_id" : self.netfig_id,
            "creator_id" : self.creator_id,
            "create_time" : str(self.create_time),
            "last_edit" : str(self.last_edit),
            "title" : self.title,
            "description" : self.description,
            "nodes" : [n.toDict() for n in self.nodes],
            "edges" : [e.toDict() for e in self.edges]
        }

        return d

class NetfigNodeType(Base):
    __tablename__ = "netfig_node_types"

    netfig_node_type = Column("netfig_node_type", String, primary_key=True)

class NetfigObjIDType(Base):
    __tablename__ = "netfig_obj_idtypes"

    netfig_obj_idtype = Column("netfig_obj_idtype", String, primary_key=True)
    netfig_obj_idtype_desc = Column("netfig_obj_idtype_description", String)


class NetfigNode(Base):
    __tablename__ = "netfig_nodes"

    netfig_node_id = Column("netfig_node_id", BigInteger, primary_key=True)
    netfig_id = Column("netfig_id", BigInteger, ForeignKey("netfig.netfig_id"))

    node_type = Column("node_type", String, ForeignKey("netfig_node_types.netfig_node_type"))
    node_obj_idtype = Column("node_obj_idtype", String, ForeignKey("netfig_obj_idtypes.netfig_obj_idtype"))
    node_obj_id = Column("node_obj_id", String)

    node_title = Column("node_title", String)
    node_description = Column("node_description", String)
    node_data = Column("node_data", String)

    node_compartment = Column("node_compartment", String)

    node_creator = Column("node_creator", Integer, ForeignKey("pd_user.pd_user_id"))
    node_create_time = Column("node_create_time", DateTime)

    node_pos_top = Column("node_pos_top", Integer)
    node_pos_left = Column("node_pos_left", Integer)

    netfig = relationship("Netfig", backref="nodes")

    def toDict(self):
        d = {
            "netfig_node_id": self.netfig_node_id,
            "netfig_id": self.netfig_id,

            "node_type": self.node_type,
            "node_obj_idtype": self.node_obj_idtype,
            "node_obj_id": self.node_obj_id,

            "node_title": self.node_title,
            "node_description": self.node_description,
            "node_data": self.node_data,

            "node_compartment": self.node_compartment,
            "node_creator": self.node_creator,
            "node_create_time": str(self.node_create_time) if self.node_create_time is not None else None,

            "node_pos_top" : self.node_pos_top,
            "node_pos_left" : self.node_pos_left
        }

        return d

class NetfigEdgeType(Base):
    __tablename__ = "netfig_edge_types"

    netfig_edge_type = Column("netfig_edge_type", String, primary_key=True)


class NetfigEdgeOrigins(Base):
    __tablename__ = "netfig_edge_origins"

    netfig_edge_source = Column("netfig_edge_origin", String, primary_key=True)
#    netfig_edge_source_desc = Column("netfig_edge_source_decc", String)


class NetfigEdge(Base):

    __tablename__ = "netfig_edges"

    netfig_edge_id = Column("netfig_edge_id", BigInteger, primary_key=True)
    netfig_id = Column("netfig_id", BigInteger, ForeignKey("netfig.netfig_id"))

    edge_type = Column("edge_type", String, ForeignKey("netfig_edge_types.netfig_edge_type"))
    directed = Column("directed", Boolean)

    edge_origin = Column("edge_origin", String, ForeignKey("netfig_edge_origins.netfig_edge_origin"))
    edge_origin_id = Column("edge_origin_id", String)

    edge_title = Column("edge_title", String)
    edge_description = Column("edge_description", String)
    edge_data = Column("edge_data", String)

    source_node_id = Column("source_node_id", BigInteger, ForeignKey("netfig_nodes.netfig_node_id"))
    target_node_id = Column("target_node_id", BigInteger, ForeignKey("netfig_nodes.netfig_node_id"))

    netfig = relationship("Netfig", backref="edges")
    source = relationship("NetfigNode", foreign_keys=[source_node_id], backref=backref("out_edges", cascade="all,delete" ) )
    target = relationship("NetfigNode", foreign_keys=[target_node_id], backref=backref("in_edges", cascade="all,delete" ))

    def toDict(self):
        d = {
            "netfig_edge_id": self.netfig_edge_id,
            "netfig_id": self.netfig_id,

            "edge_type": self.edge_type,
            "directed" : self.directed,

            "edge_origin": self.edge_origin,
            "edge_origin_id": self.edge_origin_id,

            "edge_title": self.edge_title,
            "edge_description": self.edge_description,
            "edge_data": self.edge_data,

            "source_node_id": self.source_node_id,
            "target_node_id": self.target_node_id,

        }
        return d


class NetfigMeanFoldChanges(Base):
    __tablename__ = "mean_fold_change"

    gene_symbol = Column("gene_symbol", String, primary_key=True)
    drd2_dopdepletion = Column("drd2_dopdepletion", Float)
    drd2_chronic_low = Column("drd2_chronic_low", Float)
    drd2_chronic_high = Column("drd2_chronic_high", Float)
    drd2_chronic_high_vs_low = Column("drd2_chronic_high_vs_low", Float)
    drd1a_dopdepletion = Column("drd1a_dopdepletion", Float)
    drd1a_chronic_low = Column("drd1a_chronic_low", Float)
    drd1a_chronic_high = Column("drd1a_chronic_high", Float)
    drd1a_chronic_high_vs_low = Column("drd1a_chronic_high_vs_low", Float)

    def to_dict(self):
        return {
            "gene_symbol": self.gene_symbol,
            "drd2_dopdepletion": self.drd2_dopdepletion,
            "drd2_chronic_low": self.drd2_chronic_low,
            "drd2_chronic_high": self.drd2_chronic_high,
            "drd2_chronic_high_vs_low" : self.drd2_chronic_high_vs_low,
            "drd1a_dopdepletion": self.drd1a_dopdepletion,
            "drd1a_chronic_low": self.drd1a_chronic_low,
            "drd1a_chronic_high": self.drd1a_chronic_high,
            "drd1a_chronic_high_vs_low" : self.drd1a_chronic_high_vs_low
        }


def addGeneProductNode(s, netfig_id, gene_symbol, node_compartment="cytoplasm"):
    nn = NetfigNode()
    nn.netfig_id = netfig_id
    nn.node_type = "gene product"
    nn.node_compartment = node_compartment
    nn.node_creator = 1
    nn.node_obj_idtype = "entrez_gene_symbol"
    nn.node_obj_id = gene_symbol
    nn.node_pos_left = 50
    nn.node_pos_top = 50

    s.add(nn)
    s.commit()


def addEdge(s, netfig_id, source_node_id, target_node_id, edge_type, directed, edge_origin):
    try:
        source_node = s.query(NetfigNode).filter(NetfigNode.netfig_node_id == source_node_id).one()
        target_node = s.query(NetfigNode).filter(NetfigNode.netfig_node_id == target_node_id).one()

        e = NetfigEdge()
        e.netfig_id = netfig_id
        e.source = source_node
        e.target = target_node
        e.edge_type = edge_type
        e.directed = directed
        e.edge_origin = edge_origin

        s.add(e)
        s.commit()

    except:
        print "error adding edge"
        print sys.exc_info()


def create_demo_data():
    pass