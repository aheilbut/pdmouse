CREATE TABLE motifmap_loc (
    zscore FLOAT,
    bbls FLOAT,
    fdr_upper FLOAT,
    stop BIGINT,
    fdr FLOAT,
    strand CHAR(1),
    bls FLOAT,
    mm_motif_id VARCHAR(50),
    fdr_lower FLOAT,
    cid BIGINT,
    medianhits FLOAT,
    start BIGINT,
    mm_tf_name VARCHAR(50),
    orientation INTEGER,
    chromosome VARCHAR(20),
    stdevhits FLOAT,
    lod FLOAT,
    nlod FLOAT,
    realhits INTEGER
);

CREATE TABLE motifmap_motifs (
    mm_motif_id VARCHAR(50),
    mm_motif_name VARCHAR(50),
    consensus VARCHAR(200),
    length INTEGER
);

\copy motifmap_loc from '/data/adrian/data/MotifMap/MotifMap_MOUSE_mm9.multiz30way.tsv' with delimiter as '       ' csv header;

\copy motifmap_motifs from '/data/adrian/Dropbox/Data/motifmap/motifmap_mm9_tf_motifs.txt' with delimiter as '  ' csv header;


create index mm_loc_chrom ON motifmap_loc(chromosome);
create index mm_loc_start ON motifmap_loc(start);
create index mm_loc_stop ON motifmap_loc(stop);


genename   | character varying(255) | 
 name       | character varying(255) | 
 chrom      | character varying(255) | 
 strand     | character(1)           | 
 txstart    | bigint                 | 
 txend      | bigint                 | 
 cdsstart   | bigint                 | 
 cdsend     | bigint                 | 
 exoncount  | integer                | 
 exonstarts | text                   | 
 exonends   | text                   

select genename, mm9_refflat.chrom, cdsstart 
 
 FROM mm9_refflat 
 INNER JOIN motifmap_motifs ON 
  motifmap_motifs.chromosome = mm9_refflat.chrom 
  AND
  (
  (
    mm9_refFlat.strand = '+' 
    AND 
    motifmap_motifs.start < mm9_refflat.txstart
    AND 
    motifmap_motifs.start > ( mm9_refflat.txstart - 2000 ) 
  )
  OR
  (
    mm9_refFlat.strand = '-' 
    AND 
    motifmap_motifs.end > mm9_refflat.txend
    AND 
    motifmap_motifs.end < ( mm9_refflat.txend + 2000 ) 
  
  )
  )