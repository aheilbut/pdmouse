CREATE TABLE mm9_refFlat (
  geneName VARCHAR(255),
  name VARCHAR(255),
  chrom VARCHAR(255),
  strand CHAR(1),
  txStart BIGINT,
  txEnd BIGINT,
  cdsStart BIGINT,
  cdsEnd BIGINT,
  exonCount INTEGER,
  exonStarts TEXT,
  exonEnds TEXT
);

create index i_cdsStart on mm9_refFlat(cdsStart);
create index i_cdsEnd on mm9_refFlat(cdsEnd);
create index i_geneName on mm9_refFlat(geneName);
create index i_chrom on mm9_refFlat(chrom);

\copy mm9_refFlat from '/data/adrian/Dropbox/Data/mouse_genome/refFlat.txt' with delimiter as '   ' csv header;

CREATE TABLE swissreg_sites (
	chrom VARCHAR(255),
 	alg VARCHAR(255),
  	motif_type VARCHAR(255),
	motif_start BIGINT,
	motif_end BIGINT, 
 	score FLOAT,
	strand CHAR(1),
	frame VARCHAR(5),
 	data VARCHAR(2000),
	motif_name VARCHAR(50)
)

create index i_ss_chrom on swissreg_sites(chrom);
create index i_ss_motifstart on swissreg_sites(motif_start);
create index i_ss_motifname on swissreg_sites(motif_name);


CREATE TABLE mm9_gene_motifs (
	geneName VARCHAR(255),
 	chrom VARCHAR(255),
	cdsStart BIGINT,
	cdsEnd BIGINT,
	gene_strand CHAR(1),
	motif_type VARCHAR(255),
	motif_score FLOAT,
	motif_start BIGINT,
	motif_end BIGINT,
	motif_name VARCHAR(50)
)

create index i_mg_geneName on mm9_gene_motifs(geneName);
create index i_mg_motif_name on mm9_gene_motifs(motif_name);


select distinct geneName, motif_name from mm9_gene_motifs


insert into mm9_gene_motifs

select geneName, mm9_refFlat.chrom, cdsStart, cdsEnd, mm9_refFlat.strand, motif_type, score, 
	motif_start, motif_end, motif_name
 from mm9_refFlat 
 INNER JOIN swissreg_sites
  ON 
  swissreg_sites.chrom = mm9_refFlat.chrom
  AND
  (
  (mm9_refFlat.strand = '+' 
  	AND
   swissreg_sites.motif_start < mm9_refFlat.cdsStart
  	AND
   swissreg_sites.motif_start > (mm9_refFlat.cdsStart - 10000))
  OR
  (mm9_refFlat.strand = '-' 
  	AND
   swissreg_sites.motif_start > mm9_refFlat.cdsEnd
  	AND
   swissreg_sites.motif_start < (mm9_refFlat.cdsEnd + 10000))
  )


CREATE TABLE mm9_gene_srmotifs_detail (
	geneName VARCHAR(255),
 	chrom VARCHAR(255),
	cdsStart BIGINT,
	cdsEnd BIGINT,
	txStart BIGINT,
	txEnd BIGINT,
	exonStarts TEXT,
	exonEnds TEXT,
	gene_strand CHAR(1),
	motif_strand CHAR(1),
	motif_type VARCHAR(255),
	motif_score FLOAT,
	motif_start BIGINT,
	motif_end BIGINT,
	motif_name VARCHAR(50)
)



insert into mm9_gene_srmotifs_detail

select geneName, mm9_refFlat.chrom, cdsStart, cdsEnd, txStart, txEnd, exonStarts, exonEnds,
    mm9_refFlat.strand, swissreg_sites.strand, motif_type, score, motif_start, motif_end, motif_name
 from mm9_refFlat 
 INNER JOIN swissreg_sites
  ON 
  swissreg_sites.chrom = mm9_refFlat.chrom
  AND
  (  
  ( 
   mm9_refFlat.strand = '+' 
  	AND
   swissreg_sites.motif_start > (mm9_refFlat.txStart - 10000)
  	AND
   swissreg_sites.motif_start < (mm9_refFlat.txEnd + 10000)
  )
  OR
  (
   mm9_refFlat.strand = '-' 
  	AND
   swissreg_sites.motif_start < (mm9_refFlat.txEnd + 10000)
  	AND
   swissreg_sites.motif_start > (mm9_refFlat.txStart - 10000)
  )
  )




