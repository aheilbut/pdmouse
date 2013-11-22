CREATE TABLE pd_user (
    pd_user_id INTEGER NOT NULL PRIMARY KEY,
    email VARCHAR,
    password VARCHAR,
    firstname VARCHAR,
    lastname VARCHAR,
    institution VARCHAR,
    create_time TIMESTAMP,
    last_login TIMESTAMP
);

INSERT INTO pd_user
(pd_user_id, email, password, firstname, lastname, institution)
  VALUES (1, 'heilbut@broadinstitute.org', '', 'Adrian', 'Heilbut', 'Broad Institute');
INSERT INTO pd_user
  (pd_user_id, email, password, firstname, lastname, institution)
  VALUES
(2, 'heiman@broadinstitute.org', '', 'Myriam', 'Heiman', 'Broad Institute');


CREATE TABLE netfig (
    netfig_id SERIAL NOT NULL PRIMARY KEY,
    creator INTEGER NOT NULL REFERENCES pd_user(pd_user_id),
    create_time TIMESTAMP,
    last_edit TIMESTAMP,
    title VARCHAR,
    description VARCHAR
);

CREATE TABLE netfig_node_types (
    netfig_node_type VARCHAR NOT NULL PRIMARY KEY
);

INSERT INTO netfig_node_types VALUES ('gene');
INSERT INTO netfig_node_types VALUES ('gene product');
INSERT INTO netfig_node_types VALUES ('small molecule');

CREATE TABLE netfig_obj_idtypes (
    netfig_obj_idtype VARCHAR NOT NULL PRIMARY KEY,
    netfig_obj_idtype_description VARCHAR NOT NULL
);

INSERT INTO netfig_obj_idtypes VALUES ('entrez_gene_id', 'Entrez Gene');
INSERT INTO netfig_obj_idtypes VALUES ('entrez_gene_symbol', 'Gene Symbol');


CREATE TABLE netfig_compartments (
  compartment VARCHAR NOT NULL PRIMARY KEYm
  compartment_name VARCHAR
);

INSERT INTO netfig_compartments VALUES ('extracellular', 'Extracellular');
INSERT INTO netfig_compartments VALUES ('surfmem', 'Cell Surface Membrane');
INSERT INTO netfig_compartments VALUES ('cytoplasm', 'Cytoplasm');
INSERT INTO netfig_compartments VALUES ('nucmem', 'Nuclear Membrane');
INSERT INTO netfig_compartments VALUES ('nucleus', 'Nucleus');
INSERT INTO netfig_compartments VALUES ('genome', 'Genome');

CREATE TABLE netfig_nodes (
    netfig_node_id SERIAL NOT NULL PRIMARY KEY,
    netfig_id INTEGER REFERENCES netfig(netfig_id),
    
    node_type VARCHAR REFERENCES netfig_node_types(netfig_node_type),

    node_obj_idtype VARCHAR REFERENCES netfig_obj_idtypes(netfig_obj_idtype),
    node_obj_id VARCHAR,

    node_title VARCHAR,
    node_description VARCHAR,
    node_data TEXT,

    node_pos_top FLOAT,
    node_pos_left FLOAT,
    
    node_compartment VARCHAR,
    
    node_creator INTEGER REFERENCES pd_user(pd_user_id),
    node_create_time TIMESTAMP
);

CREATE TABLE netfig_edge_origins (
  netfig_edge_origin VARCHAR NOT NULL PRIMARY KEY
);

CREATE TABLE netfig_edges (
    netfig_edge_id SERIAL NOT NULL PRIMARY KEY,
    netfig_id INTEGER REFERENCES netfig(netfig_id),

    directed BOOLEAN NOT NULL,

    edge_origin VARCHAR REFERENCES netfig_edge_origins(netfig_edge_origin),
    edge_origin_id VARCHAR,

    edge_type VARCHAR REFERENCES netfig_edge_types(netfig_edge_type),

    edge_title VARCHAR,
    edge_description VARCHAR,
    edge_data TEXT,

    source_node_id BIGINT NOT NULL REFERENCES netfig_nodes(netfig_node_id),
    target_node_id BIGINT NOT NULL REFERENCES netfig_nodes(netfig_node_id),

    edge_creator INTEGER REFERENCES pd_user(pd_user_id),
    edge_create_time TIMESTAMP   
);


    source_node_id BIGINT REFERENCES netfig_nodes(netfig_node_id),
    target_node_id BIGINT REFERENCES netfig_nodes(netfig_node_id),

CREATE TABLE mean_fold_change (
  gene_symbol VARCHAR PRIMARY KEY,
  drd2_dopdepletion FLOAT,
  drd2_chronic_low FLOAT,
  drd2_chronic_high FLOAT,
  drd1a_dopdepletion FLOAT,
  drd1a_chronic_low FLOAT,
  drd1a_chronic_high FLOAT
);

\copy mean_fold_change FROM '/data/adrian/Dropbox/Projects/Broad/PD_mouse/results/2013_11_18/mean_fold_change.tab' WITH DELIMITER '       ' CSV HEADER


CREATE TABLE netfig_edge_types (
  netfig_edge_type VARCHAR NOT NULL PRIMARY KEY
);

INSERT INTO netfig_edge_types VALUES ('binds');
INSERT INTO netfig_edge_types VALUES ('activates');
INSERT INTO netfig_edge_types VALUES ('inhibits');
INSERT INTO netfig_edge_types VALUES ('regulates');
INSERT INTO netfig_edge_types VALUES ('in complex');
INSERT INTO netfig_edge_types VALUES ('co-regulated');
INSERT INTO netfig_edge_types VALUES ('co-occurs in literature');


CREATE TABlE pd_log (
    
);


    source_node_type VARCHAR REFERENCES netfig_node_types(netfig_node_type),
    source_node_id VARCHAR NOT NULL,

    source_node_obj_idtype VARCHAR REFERENCES netfig_obj_idtypes(netfig_obj_idtype),
    source_node_obj_id VARCHAR NOT NULL,

    target_node_type VARCHAR REFERENCES netfig_node_types(netfig_node_type),
    target_node_id VARCHAR NOT NULL,

    target_node_obj_idtype VARCHAR REFERENCES netfig_obj_idtypes(netfig_obj_idtype),
    target_node_obj_id VARCHAR,
