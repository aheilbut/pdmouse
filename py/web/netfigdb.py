
from sqlalchemy import create_engine
from sqlalchemy import Table, Column, Integer, String, Date, BigInteger, DateTime, Boolean, MetaData, ForeignKey, UniqueConstraint, ForeignKeyConstraint
from sqlalchemy.orm import mapper, sessionmaker, relationship, backref

import datetime

import pd_locals

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

    netfig = relationship("Netfig", backref="nodes")

    def toDict(self):
        d = {
            "node_id" : self.netfig_node_id,
            "netfig_id" : self.netfig_id,

            "node_type" : self.node_type,
            "node_obj_idtype" : self.node_obj_idtype,
            "node_obj_id" : self.node_obj_id,

            "node_title" : self.node_title,
            "node_description" : self.node_description,
            "node_data" : self.node_data,

            "node_compartment" : self.node_compartment,
            "node_creator" : self.node_creator,
            "node_create_time" : str(self.node_create_time),
        }

        return d

class NetfigEdgeType(Base):
    __tablename__ = "netfig_edge_types"

    netfig_edge_type = Column("netfig_edge_type", String, primary_key=True)

class NetfigEdgeSource(Base):
    __tablename__ = "netfig_edge_sources"

    netfig_edge_source = Column("netfig_edge_source", String, primary_key=True)
    netfig_edge_source_desc = Column("netfig_edge_source_decc", String)


class NetfigEdge(Base):

    __tablename__ = "netfig_edges"

    netfig_edge_id = Column("netfig_edge_id", BigInteger, primary_key=True)
    netfig_id = Column("netfig_id", BigInteger, ForeignKey("netfig.netfig_id"))

    edge_type = Column("edge_type", String, ForeignKey("netfig_edge_types.netfig_edge_type"))
    directed = Column("directed", Boolean)

    edge_source = Column("edge_source", String, ForeignKey("netfig_edge_sources.netfig_edge_source"))
    edge_source_id = Column("edge_source_id", String)

    edge_title = Column("edge_title", String)
    edge_description = Column("edge_description", String)
    edge_data = Column("edge_data", String)

    source_node_id = Column("source_node_id", BigInteger, ForeignKey("netfig_nodes.netfig_node_id"))
    target_node_id = Column("target_node_id", BigInteger, ForeignKey("netfig_nodes.netfig_node_id"))

    source_node_obj_idtype = Column("source_node_obj_idtype", String,
                                    ForeignKey("netfig_obj_idtypes.netfig_obj_idtype"))
    source_node_obj_id = Column("source_node_obj_id", String)

    target_node_obj_idtype = Column("target_node_obj_idtype", String,
                                ForeignKey("netfig_obj_idtypes.netfig_obj_idtype"))
    target_node_obj_id = Column("target_node_obj_id", String)

    netfig = relationship("Netfig", backref="edges")
    source = relationship("NetfigNode", foreign_keys=[source_node_id], backref=backref("out_edges", cascade="all,delete" ) )
    target = relationship("NetfigNode", foreign_keys=[target_node_id], backref=backref("in_edges", cascade="all,delete" ))

    def toDict(self):
        d = {
            "netfig_edge_id" : self.netfig_edge_id,
            "netfig_id" : self.netfig_id,

            "edge_type" : self.edge_type,

            "edge_source" : self.edge_source,
            "edge_source_id" : self.edge_source_id,

            "edge_title" : self.edge_title,
            "edge_description" : self.edge_description,

            "edge_data" : self.edge_data,
            "source_node_id" : self.source_node_id,
            "target_node_id" : self.target_node_id,

            "source_node_obj_idtype" : self.source_node_obj_idtype,
            "source_node_obj_id" : self.source_node_obj_id,
            "target_node_obj_idtype" : self.target_node_obj_idtype,
            "target_node_obj_id" : self.target_node_obj_id

        }
        return d