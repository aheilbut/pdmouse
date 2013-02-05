require(affy)
require(samr)
require(mouse4302.db)

mo430genenames <- toTable(mouse4302GENENAME)
mo430symbol <- toTable(mouse4302SYMBOL)

loadCELs <- function() {
  cp101files <- list.celfiles("~/Dropbox/Projects/Broad/PD_mouse/CP101", full.names=TRUE)
  cp73files <- list.celfiles("~/Dropbox/Projects/Broad/PD_mouse/CP73", full.names=TRUE)
  
  allfiles <- c(cp73files, cp101files)
  
  cel_data <- affy::ReadAffy(filenames=allfiles)
  pd.eset <- affy::rma(cel_data)
  return(pd.eset)
}

pd.eset <- loadCELs()

# make a simple linear model
aimcorr_cp73_chron_high_f <- lapply(1:length(rownames(exprs(pd.eset))), 
  function(i) {
    return( lm(aim ~ expression, 
              data.frame(expression = exprs(pd.eset)[i,d_cp73_chronic_high_levo$filenames], 
                         aim = d_cp73_chronic_high_levo$totalAIM)) )
  })

aimcorr_cp73_chron_low_f <- lapply(1:length(rownames(exprs(pd.eset))), 
  function(i) {
    return( lm(aim ~ expression, 
               data.frame(expression=exprs(pd.eset)[i,d_cp73_chronic_low_levo$filenames], 
                          aim=d_cp73_chronic_low_levo$totalAIM)) )    
  })

aimcorr_cp73_chron_high_f_pval <- sapply( aimcorr_cp73_chron_high_f, lmp)
aimcorr_cp73_chron_low_f_pval <- sapply( aimcorr_cp73_chron_low_f, lmp)

mo430genenames[ rownames(exprs(pd.eset))[ which( aimcorr_cp73_chron_high_f_pval < 0.001 ) ], ]

plotcorr <- function(x) {
  plot( aimcorr_cp73_chron_high_f[[x]]$model$expression, aimcorr_cp73_chron_high_f[[x]]$model$aim)
  lines( aimcorr_cp73_chron_high_f[[x]]$model$expression, aimcorr_cp73_chron_high_f[[x]]$fitted.values, lty=8, col="red")
}


loadData <- function() { 
  
  pd.covar <- read.table("~/Dropbox/Projects/Broad/PD_mouse/behav_data.tab", sep="\t", comment.char="^", header=TRUE)
  files <- data.frame(celfile=sampleNames(pd.eset))
  filesplit <- data.frame(matrix(unlist(strsplit(as.character(files$celfile), "_|\\.")), nrow=nrow(files), byrow=T), row.names=files$celfile)
  colnames(filesplit) <- c("Lesion", "P", "N", "MouseType", "MouseID","Filetype")
  filesplit$filenames <- rownames(filesplit)
  pd.covar <- merge(pd.covar, filesplit, by.x=c("MouseType","MouseID"),by.y=c("MouseType", "MouseID"))
  
  pd.edge.covar <- t(pd.covar)
  colnames(pd.edge.covar) <- pd.covar$filenames
  
  write.table(pd.edge.covar, 
              "~/Dropbox/Projects/Broad/PD_mouse/pd_covar_edge.tab", 
              quote=FALSE,
              sep="\t")
  return(pd.covar)
                           
}

pd.covar <- loadData()
pd.covar$totalAIM <- apply(pd.covar[c("Day1AIM","Day3AIM","Day4AIM","Day5AIM","Day8AIM")], 1, sum, na.rm=TRUE)


ascorbate <- subset(pd.covar, MouseType=="CP101" & DrugTreat=="Chronic saline" & LesionType=="Ascorbate")
ohda <- subset(pd.covar, MouseType=="CP101" & DrugTreat=="Chronic saline" & LesionType=="6-OHDA")

ascorbate.data <- exprs(pd.eset)[,ascorbate$filenames]
ohda.data <- exprs(pd.eset)[,ohda$filenames]

dat.medians <- data.frame(ascorbate=apply(ascorbate.data, 1, median), 
                          ascorbate.sd=apply(ascorbate.data, 1, sd),
                          ohda=apply(ohda.data, 1, median),
                          ohda.sd=apply(ohda.data, 1, sd)
                          )


dat.medians$ttest <- apply(exprs(pd.eset), 1, 
  function(x) {
    t.test(x[ohda$filenames], x[ascorbate$filenames])$p.value
  })


data <- exprs(pd.eset[,c(ascorbate$filenames, ohda$filenames)])
y <- c(rep(1,10), rep(2,7))
samfit <- SAM(data, y, resp.type="Two class unpaired", geneid=row.names(data), nperms=50)

# 6-OHDA/saline vs 6-OHDA/chronic low levodopa vs. chronic high levodopa
d_cp73_chronic_saline <- subset(pd.covar, MouseType=="CP73" & DrugTreat=="Chronic saline" & LesionType=="6-OHDA")
d_cp73_chronic_low_levo <- subset(pd.covar, MouseType=="CP73" & DrugTreat=="Chronic low levodopa" & LesionType=="6-OHDA")
d_cp73_chronic_high_levo <- subset(pd.covar, MouseType=="CP73" & DrugTreat=="Chronic high levodopa" & LesionType=="6-OHDA")
d_cp73_acute_high_levo <- subset(pd.covar, MouseType=="CP73" & DrugTreat=="Acute high levodopa" & LesionType=="6-OHDA")
d_cp73_acute_saline <- subset(pd.covar, MouseType=="CP73" & DrugTreat=="Acute saline" & LesionType=="6-OHDA")

d_cp101_chronic_saline <- subset(pd.covar, MouseType=="CP101" & DrugTreat=="Chronic saline" & LesionType=="6-OHDA")
d_cp101_chronic_low_levo <- subset(pd.covar, MouseType=="CP101" & DrugTreat=="Chronic low levodopa" & LesionType=="6-OHDA")
d_cp101_chronic_high_levo <- subset(pd.covar, MouseType=="CP101" & DrugTreat=="Chronic high levodopa" & LesionType=="6-OHDA")
d_cp101_acute_high_levo <- subset(pd.covar, MouseType=="CP101" & DrugTreat=="Acute high levodopa" & LesionType=="6-OHDA")
d_cp101_acute_saline <- subset(pd.covar, MouseType=="CP101" & DrugTreat=="Acute saline" & LesionType=="6-OHDA")

d_cp101_all_levo <- subset(pd.covar, MouseType=="CP101" & (DrugTreat=="Chronic high levodopa" | DrugTreat=="Chronic low levodopa" ) & LesionType=="6-OHDA")
d_cp73_all_levo <- subset(pd.covar, MouseType=="CP73" & (DrugTreat=="Chronic high levodopa" | DrugTreat=="Chronic low levodopa" ) & LesionType=="6-OHDA")

d_cp101_all_saline <- subset(pd.covar, MouseType=="CP101" & (DrugTreat=="Acute saline" | DrugTreat=="Chronic saline" ) & LesionType=="6-OHDA")
d_cp73_all_saline <- subset(pd.covar, MouseType=="CP73" & (DrugTreat=="Acute saline" | DrugTreat=="Chronic saline" ) & LesionType=="6-OHDA")

d_allMice_allChronicLevo <- subset(pd.covar, (MouseType=="CP73" | MouseType=="CP101" ) & (DrugTreat=="Chronic high levodopa" | DrugTreat=="Chronic low levodopa" ) & LesionType=="6-OHDA")


for (i in 1:8) { 
  par(new=T); 
  plot( x = c(1,3,4,5,8), 
        y = d_c[i, c("Day1AIM","Day3AIM", "Day4AIM", "Day5AIM", "Day8AIM")], 
        type="o", 
        col="red", 
        ylim=c(0, 40)) 
}


function doSAM(data, response, comparison) { 

data <- exprs(pd.eset[,c(d_b$filenames, d_c$filenames)])
y <- c(d_b$totalAIM, d_c$totalAIM)

data <- exprs(pd.eset[,d_chronic_high_levo$filenames])
y <- c(d_chronic_high_levo$totalAIM)

fdr_cutoff = 0.2


SAM.bothMice.chronic_levo_vs_AIM <- SAM(exprs(pd.eset[,d_allMice_allChronicLevo$filenames]),
                d_allMice_allChronicLevo$totalAIM,
                resp.type="Quantitative", 
                geneid=row.names(exprs(pd.eset)),  
                fdr.output=0.10, nperm=1000)

SAM.bothMice.chronic_levo_vs_AIMDay8 <- SAM(exprs(pd.eset[,d_allMice_allChronicLevo$filenames]),
                d_allMice_allChronicLevo$Day8AIM,
                resp.type="Quantitative", 
                geneid=row.names(exprs(pd.eset)),  
                fdr.output=0.15, nperm=200)

SAM.cp73.chronic_high_vs_AIM <- SAM(exprs(pd.eset[,d_cp73_chronic_high_levo$filenames]),
                d_cp73_chronic_high_levo$totalAIM,
                resp.type="Quantitative", 
                geneid=row.names(exprs(pd.eset)),  
                fdr.output=fdr_cutoff)

SAM.cp73.all_levo_vs_AIM <- SAM(exprs(pd.eset[,d_cp73_all_levo$filenames]),
                                 d_cp73_all_levo$totalAIM,
                                 resp.type="Quantitative", 
                                 geneid=row.names(exprs(pd.eset)),
                                fdr.output=fdr_cutoff, nperm=1000)

SAM.cp73.all_levo_vs_AIM.df.UP <- as.data.frame( SAM.cp73.all_levo_vs_AIM$siggenes.table$genes.up, stringsAsFactors=FALSE)
SAM.cp73.all_levo_vs_AIM.df.DN <- as.data.frame( SAM.cp73.all_levo_vs_AIM$siggenes.table$genes.lo, stringsAsFactors=FALSE)


SAM.cp101.chronic_high_vs_AIM <- SAM(exprs(pd.eset[,d_cp101_chronic_high_levo$filenames]),
                                     d_cp101_chronic_high_levo$totalAIM,
                                     resp.type="Quantitative", 
                                     geneid=row.names(exprs(pd.eset)),
                                     fdr.output=fdr_cutoff)

SAM.cp101.chronic_low_vs_AIM <- SAM(exprs(pd.eset[,d_cp101_chronic_low_levo$filenames]),
                                     d_cp101_chronic_low_levo$totalAIM,
                                     resp.type="Quantitative", 
                                     geneid=row.names(exprs(pd.eset)),
                                     fdr.output=fdr_cutoff)

SAM.cp73.chronic_low_vs_AIM <- SAM(exprs(pd.eset[,d_cp73_chronic_low_levo$filenames]),
                                     d_cp73_chronic_low_levo$totalAIM,
                                     resp.type="Quantitative", 
                                     geneid=row.names(exprs(pd.eset)),
                                     fdr.output=fdr_cutoff)


SAM.cp101.all_levo_vs_AIM <- SAM(exprs(pd.eset[,d_cp101_all_levo$filenames]),
                                 d_cp101_all_levo$totalAIM,
                                 resp.type="Quantitative", 
                                 geneid=row.names(exprs(pd.eset)),
                                 fdr.output=fdr_cutoff, nperm=1000)
SAM.cp101.all_levo_vs_AIM.df.UP <- as.data.frame( SAM.cp101.all_levo_vs_AIM$siggenes.table$genes.up, stringsAsFactors=FALSE)
SAM.cp101.all_levo_vs_AIM.df.DN <- as.data.frame( SAM.cp101.all_levo_vs_AIM$siggenes.table$genes.lo, stringsAsFactors=FALSE)



SAM.cp101.all_levo_vs_saline <- SAM(exprs(pd.eset[,c(d_cp101_all_levo$filenames,d_cp101_all_saline$filenames)]),
                                 c( rep(2, dim(d_cp101_all_levo)[1]), rep(1, dim(d_cp101_all_saline)[1])),
                                 resp.type="Two class unpaired", 
                                 geneid=row.names(exprs(pd.eset)),
                                 fdr.output=fdr_cutoff)

SAM.cp73.all_levo_vs_saline <- SAM(exprs(pd.eset[,c(d_cp73_all_levo$filenames,d_cp73_all_saline$filenames)]),
                                    c( rep(2, dim(d_cp73_all_levo)[1]), rep(1, dim(d_cp73_all_saline)[1])),
                                    resp.type="Two class unpaired", 
                                    geneid=row.names(exprs(pd.eset)),
                                    fdr.output=fdr_cutoff)

levo_vs_AIM_bothUP <- annotate_samtable( merge(SAM.cp101.all_levo_vs_AIM.df.UP, SAM.cp73.all_levo_vs_AIM.df.UP, by="Gene Name") )
levo_vs_AIM_bothDN <- annotate_samtable( merge(SAM.cp101.all_levo_vs_AIM.df.DN, SAM.cp73.all_levo_vs_AIM.df.DN, by="Gene Name") )


SAM.cp101.chronic_high_vs_saline <- SAM(exprs(pd.eset[,c(d_cp101_chronic_high_levo$filenames, d_cp101_chronic_saline$filenames)]),
                                        c(rep(2, dim(d_cp101_chronic_high_levo)[1]), rep(1, dim(d_cp101_chronic_saline)[1])),
                                        resp.type="Two class unpaired", 
                                        geneid=row.names(exprs(pd.eset)),
                                        fdr.output=0.20, nperm=100)  

SAM.cp73.chronic_high_vs_saline <- SAM(exprs(pd.eset[,c(d_cp73_chronic_high_levo$filenames, d_cp73_chronic_saline$filenames)]),
                                        c(rep(2, dim(d_cp73_chronic_high_levo)[1]), rep(1, dim(d_cp73_chronic_saline)[1])),
                                        resp.type="Two class unpaired", 
                                        geneid=row.names(exprs(pd.eset)),
                                        fdr.output=0.20, nperm=100)


SAM.cp73.chronicHigh_vs_chronicLow <- SAM(exprs(pd.eset[,c(d_cp73_chronic_high_levo$filenames, d_cp73_chronic_low_levo$filenames)]),
                                        c(rep(2, dim(d_cp73_chronic_high_levo)[1]), rep(1, dim(d_cp73_chronic_low_levo)[1])),
                                        resp.type="Two class unpaired", 
                                        geneid=row.names(exprs(pd.eset)),
                                        fdr.output=0.20, nperm=1500)
SAM.cp73.chronicHigh_vs_chronicLow.df.DN <- as.data.frame( SAM.cp73.chronicHigh_vs_chronicLow$siggenes.table$genes.lo, stringsAsFactors=FALSE)
SAM.cp73.chronicHigh_vs_chronicLow.df.UP <- as.data.frame( SAM.cp73.chronicHigh_vs_chronicLow$siggenes.table$genes.up, stringsAsFactors=FALSE)


SAM.cp101.chronicHigh_vs_chronicLow <- SAM(exprs(pd.eset[,c(d_cp101_chronic_high_levo$filenames, d_cp101_chronic_low_levo$filenames)]),
                                        c(rep(2, dim(d_cp101_chronic_high_levo)[1]), rep(1, dim(d_cp101_chronic_low_levo)[1])),
                                        resp.type="Two class unpaired", 
                                        geneid=row.names(exprs(pd.eset)),
                                        fdr.output=0.20, nperm=1500)
SAM.cp101.chronicHigh_vs_chronicLow.df.UP <- as.data.frame( SAM.cp101.chronicHigh_vs_chronicLow$siggenes.table$genes.up, stringsAsFactors=FALSE)
SAM.cp101.chronicHigh_vs_chronicLow.df.DN <- as.data.frame( SAM.cp101.chronicHigh_vs_chronicLow$siggenes.table$genes.lo, stringsAsFactors=FALSE)

SAM.cp101.chronhighup_chronhighlevo_vs_AIM <- SAM(exprs(pd.eset[SAM.cp101.chronic_high_vs_saline$siggenes.table$genes.up[,"Gene Name"],
                                                                d_cp101_chronic_high_levo$filenames]),
                                                d_cp101_chronic_high_levo$totalAIM,
                                 resp.type="Quantitative", 
                                 geneid=row.names(exprs(pd.eset[SAM.cp101.chronic_high_vs_saline$siggenes.table$genes.up[,"Gene Name"],
                                                                d_cp101_chronic_high_levo$filenames])),
                                 fdr.output=0.20)


SAM.cp101.chronhighup_chronhighlevo_vs_AIMDay8 <- SAM(exprs(pd.eset[SAM.cp101.chronic_high_vs_saline$siggenes.table$genes.up[,"Gene Name"],
                                                                d_cp101_chronic_high_levo$filenames]),
                                                  d_cp101_chronic_high_levo$Day8AIM,
                                                  resp.type="Quantitative", 
                                                  geneid=row.names(exprs(pd.eset[SAM.cp101.chronic_high_vs_saline$siggenes.table$genes.up[,"Gene Name"],
                                                                                 d_cp101_chronic_high_levo$filenames])),
                                                  fdr.output=0.20)


SAM.cp73.chronhighup_chronhighlevo_vs_AIM <- SAM(exprs(pd.eset[SAM.cp73.chronic_high_vs_saline$siggenes.table$genes.up[,"Gene Name"],d_cp73_chronic_high_levo$filenames]),
                                                  d_cp73_chronic_high_levo$totalAIM,
                                                  resp.type="Quantitative", 
                                                  geneid=row.names(exprs(pd.eset[SAM.cp73.chronic_high_vs_saline$siggenes.table$genes.up[,"Gene Name"],
                                                                                 d_cp73_chronic_high_levo$filenames])),
                                                  fdr.output=0.20)


SAM.cp73.chronhighup_chronhighlevo_vs_AIMDay8 <- SAM(exprs(pd.eset[SAM.cp73.chronic_high_vs_saline$siggenes.table$genes.up[,"Gene Name"],d_cp73_chronic_high_levo$filenames]),
                                                 d_cp73_chronic_high_levo$Day8AIM,
                                                 resp.type="Quantitative", 
                                                 geneid=row.names(exprs(pd.eset[SAM.cp73.chronic_high_vs_saline$siggenes.table$genes.up[,"Gene Name"],
                                                                                d_cp73_chronic_high_levo$filenames])),
                                                 fdr.output=0.20)



SAM.cp73.chronhigh_down_chronhighlevo_vs_AIM <- SAM(exprs(pd.eset[SAM.cp73.chronic_high_vs_saline$siggenes.table$genes.up[,"Gene Name"],
                                                                  d_cp73_chronic_high_levo$filenames]),
                                                 d_cp73_chronic_high_levo$totalAIM,
                                                 resp.type="Quantitative", 
                                                 geneid=row.names(exprs(pd.eset[SAM.cp73.chronic_high_vs_saline$siggenes.table$genes.lo[,"Gene Name"],
                                                                                d_cp73_chronic_high_levo$filenames])),
                                                 fdr.output=0.20)


SAM.cp73.chronhigh_down_chronhighlevo_vs_AIMDay8 <- SAM(exprs(pd.eset[SAM.cp73.chronic_high_vs_saline$siggenes.table$genes.up[,"Gene Name"],d_cp73_chronic_high_levo$filenames]),
                                                     d_cp73_chronic_high_levo$Day8AIM,
                                                     resp.type="Quantitative", 
                                                     geneid=row.names(exprs(pd.eset[SAM.cp73.chronic_high_vs_saline$siggenes.table$genes.lo[,"Gene Name"],
                                                                                    d_cp73_chronic_high_levo$filenames])),
                                                     fdr.output=0.20)



responsive_set <- SAM.cp73.all_levo_vs_saline$siggenes.table$genes.up[as.numeric(as.data.frame(SAM.cp73.all_levo_vs_saline$siggenes.table$genes.up, stringsAsFactors=FALSE )$'q-value(%)') < 10, "Gene Name"  ]
SAM.cp73.all_levo_up_and_all_levo_AIM <- SAM(exprs(pd.eset[responsive_set,d_cp73_all_levo$filenames]),
                                            d_cp73_all_levo$totalAIM,
                                             resp.type="Quantitative", 
                                             geneid=responsive_set,
                                            fdr.output=0.05, nperm=1000)


responsive_set <- SAM.cp73.all_levo_vs_saline$siggenes.table$genes.lo[as.numeric(as.data.frame(SAM.cp73.all_levo_vs_saline$siggenes.table$genes.lo, stringsAsFactors=FALSE )$'q-value(%)') < 10, "Gene Name"  ]
SAM.cp73.all_levo_down_and_all_levo_AIM <- SAM(exprs(pd.eset[responsive_set,d_cp73_all_levo$filenames]),
                                            d_cp73_all_levo$totalAIM,
                                             resp.type="Quantitative", 
                                             geneid=responsive_set,
                                            fdr.output=0.05, nperm=1000)


responsive_set <- SAM.cp101.all_levo_vs_saline$siggenes.table$genes.up[as.numeric(as.data.frame(SAM.cp101.all_levo_vs_saline$siggenes.table$genes.up, stringsAsFactors=FALSE )$'q-value(%)') < 10, "Gene Name"  ]
SAM.cp101.all_levo_up_and_all_levo_AIM <- SAM(exprs(pd.eset[responsive_set,d_cp101_all_levo$filenames]),
                                             d_cp101_all_levo$totalAIM,
                                             resp.type="Quantitative", 
                                             geneid=responsive_set,
                                             fdr.output=0.05, nperm=1000)

responsive_set <- SAM.cp101.all_levo_vs_saline$siggenes.table$genes.lo[as.numeric(as.data.frame(SAM.cp101.all_levo_vs_saline$siggenes.table$genes.lo, stringsAsFactors=FALSE )$'q-value(%)') < 10, "Gene Name"  ]
SAM.cp101.all_levo_down_and_all_levo_AIM <- SAM(exprs(pd.eset[responsive_set,d_cp101_all_levo$filenames]),
                                            d_cp101_all_levo$totalAIM,
                                             resp.type="Quantitative", 
                                             geneid=responsive_set,
                                            fdr.output=0.05, nperm=1000)



# plot gene VS. AIM scores for 

bothUP_doseDep <- merge(SAM.cp73.chronicHigh_vs_chronicLow.df.UP, SAM.cp101.chronicHigh_vs_chronicLow.df.UP, by="Gene Name")
bothUP_doseDep <- bothUP_doseDep[ order(as.numeric(bothUP_doseDep[,"q-value(%).x"]) + as.numeric(bothUP_doseDep[,"q-value(%).y"])), ]

bothDN_doseDep <- merge(SAM.cp73.chronicHigh_vs_chronicLow.df.DN, SAM.cp101.chronicHigh_vs_chronicLow.df.DN, by="Gene Name")
bothDN_doseDep <- bothDN_doseDep[ order(as.numeric(bothDN_doseDep[,"q-value(%).x"]) + as.numeric(bothDN_doseDep[,"q-value(%).y"])), ]

UPDN_doseDep <- merge(SAM.cp73.chronicHigh_vs_chronicLow.df.UP, SAM.cp101.chronicHigh_vs_chronicLow.df.DN, by="Gene Name")
DNUP_doseDep <- merge(SAM.cp73.chronicHigh_vs_chronicLow.df.DN, SAM.cp101.chronicHigh_vs_chronicLow.df.UP, by="Gene Name")

write.table(annotate_samtable(bothDN_doseDep),
            file="~/Dropbox/bothDN_doseDep_oct22.xls", quote=FALSE, sep="\t")

write.table(annotate_samtable(bothUP_doseDep),
            file="~/Dropbox/bothUP_doseDep_oct22.xls", quote=FALSE, sep="\t")

write.table(annotate_samtable(DNUP_doseDep),
            file="~/Dropbox/DNUP_doseDep_oct22.xls", quote=FALSE, sep="\t")



all_levo_files <- c(d_chronic_low_levo$filenames, d_chronic_high_levo$filenames)
data.all_levo <- exprs(pd.eset[,all_levo_files])
responses.all_levo <- c(d_chronic_low_levo$totalAIM, d_chronic_high_levo$totalAIM)
SAM.allLevo_vs_AIM <- SAM( data.all_levo,
                           responses.all_levo,
                           resp.type="Quantitative",
                           geneid=row.names(data.all_levo)
                           )

up <- as.data.frame(samfit.2$siggenes.table$genes.up)
up$probe_id <- rownames(exprs(pd.eset))[as.integer(as.matrix(up$"Gene Name"))]

down <- as.data.frame(samfit$siggenes.table$genes.lo)
down$probe_id <- rownames(exprs(pd.eset))[down$"Gene Name"]

merge(up, mo430genenames, by.x="Gene Name", by.y="probe_id", sort=FALSE)[1:50,]
merge(dat.medians[order(dat.medians$ttest, decreasing=FALSE),][1:30,], mo430genenames, by.x=0, by.y="probe_id", sort=FALSE

annotate_samtable <- function(sam_result) {
  return( merge( merge(sam_result, 
               mo430genenames,
               by.x="Gene Name", 
               by.y="probe_id",
               all.x=TRUE,
               sort=FALSE,), mo430symbol, by.x="Gene Name", by.y="probe_id", all.x=TRUE, sort=FALSE)[1:min(10000, dim(sam_result)[1]),] )
  
}
      
print_samtable <- function(sam_result) {
  print("Down:")
  if (!is.null(sam_result$siggenes.table$genes.lo)) {
    print(  annotate_samtable(sam_result$siggenes.table$genes.lo) )
  }
      
  print("Up:")
  if (!is.null(sam_result$siggenes.table$genes.up)) {
    print(  annotate_samtable(sam_result$siggenes.table$genes.up) )
    
    }      
}
      
      
plot_ExprAIM <- function(probe_id, mouse_type) { 
 genename <- mo430genenames[mo430genenames$probe_id == probe_id,"gene_name"]
 symbol <-   mo430symbol[mo430symbol$probe_id == probe_id, "symbol"]

  plist = list()
 
  title <- paste(paste(mouse_type, probe_id, symbol, sep=" "), genename, sep="\n")
  pd_set <- subset(pd.covar, MouseType==mouse_type)
#  print(exprs(pd.eset)[probe_id,pd_set$filenames])
#  print(probe_id)
#  print(pd_set)
 
  p <- qplot( x = exprs(pd.eset)[probe_id,pd_set$filenames], 
       y=pd_set$totalAIM, 
       color=pd_set$DrugTreat, 
       shape=pd_set$LesionType, 
       size=I(3), 
       xlab="log2 expression", 
        ylab="integral of AIM score",
       main=title, 
       geom="jitter"
#       xlim=c(0,15)
        ) 
  
  return(p)
 }      
      
boxplot_Expression <- function(probe_id) {
  genename <- mo430genenames[mo430genenames$probe_id == probe_id,"gene_name"]
  symbol <-   mo430symbol[mo430symbol$probe_id == probe_id, "symbol"]
  title <- paste( paste(probe_id, symbol, sep=" "), genename, sep="\n")
  data <- merge( melt(exprs(pd.eset)[probe_id,]), pd.covar, by.x="row.names", by.y="filenames")
  par(mar = c(10,4,4,2) + 0.1)
  boxplot(value ~ DrugTreat + MouseType, data, las=2, xlab="treatment group", ylab="log2 expression", cex.axis=1, main=title)  
  beeswarm(value ~ DrugTreat + MouseType, data, las=2, add=TRUE)
}
      
      
#for 
   
printProbeGraphs <- function(probe) {
    do.call( arrange_ggplot2, list( plot_ExprAIM(probe, "CP73"),  plot_ExprAIM(probe, "CP101"), ncol=1))
    boxplot_Expression(probe)   
}      
      
printGraphForList <- function(probeList, filename) {
      pdf(file=paste("~/Dropbox/", filename, sep=""))
      for (p in probeList) {
        print(p)
        printProbeGraphs(p)        
      }
      dev.off()
}      
      
      
printGraphs <- function(samset, filename, description, maxgraphs) {       
  pdf(file=paste("~/Dropbox/", filename, sep=""))
  current <- samset
  #grid.table(annotate_samtable(current$siggenes.table$genes.up), gp=gpar(fontsize=6))      
  current.up <- as.data.frame(samset$siggenes.table$genes.up, stringsAsFactors=FALSE)
  current.lo <- as.data.frame(samset$siggenes.table$genes.lo, stringsAsFactors=FALSE )
  
  if (!is.null(current$siggenes.table$genes.up)) {
  for (g in rownames(current.up[1:min(maxgraphs, dim(current.up)[1]),])) { 
    probe <- current.up[g, "Gene Name"]
    print(probe)
    printProbeGraphs(probe)
    #do.call(arrange_ggplot2, list(do.call(plot_ExprAIM, list(probe, "CP73")), do.call(plot_ExprAIM, list(probe, "CP101")), ncol=1))
  }
  }
      
  if (!is.null(current$siggenes.table$genes.lo)) {
  for (g in rownames(current.lo[1:min(maxgraphs, dim(current.lo)[1]),])) { 
    probe <- current.lo[g, "Gene Name"]
    print(probe)
#    do.call( arrange_ggplot2, list( do.call(plot_ExprAIM, list(probe, "CP73")), do.call(plot_ExprAIM, list(probe, "CP101")), ncol=1))
    arrange_ggplot2( plot_ExprAIM(probe, "CP73"),  plot_ExprAIM(probe, "CP101"), ncol=1)
    boxplot_Expression(probe)
#    plot_ExprAIM(probe, "CP73")
#    plot_ExprAIM(probe, "CP101")
    
  }
  }
  dev.off()      
}
   

printGraphs(SAM.cp101.chronhighup_chronhighlevo_vs_AIM, "sam_cp101_chron_high_levo_vs_saline_VS_AIM.pdf", "Filter genes up vs saline, then scored by regression to AIM", 10)
      
printGraphs(SAM.cp73.chronhighup_chronhighlevo_vs_AIM, "sam_cp73_chron_high_levo_vs_saline_VS_AIM.pdf", "Filter genes up vs saline, then scored by regression to AIM", 10)

printGraphs(SAM.cp73.chronhighup_chronhighlevo_vs_AIM, "sam_cp73_chron_high_levo_vs_saline_VS_AIM.pdf", "Filter genes up vs saline, then scored by regression to AIM", 10)
      
printGraphs(SAM.cp73.all_levo_up_and_all_levo_AIM, "SAM.cp73.all_levo_up_and_all_levo_AIM.pdf", "Filter genes up vs saline, then scored by regression to AIM", 100)
printGraphs(SAM.cp101.all_levo_up_and_all_levo_AIM, "SAM.cp101.all_levo_up_and_all_levo_AIM.pdf", "Filter genes up vs saline, then scored by regression to AIM", 100)
 
printGraphs(SAM.cp73.all_levo_down_and_all_levo_AIM, "SAM.cp73.all_levo_down_and_all_levo_AIM.pdf", "Filter genes up vs saline, then scored by regression to AIM", 100)
printGraphs(SAM.cp101.all_levo_down_and_all_levo_AIM, "SAM.cp101.all_levo_down_and_all_levo_AIM.pdf", "Filter genes up vs saline, then scored by regression to AIM", 100)
      
annotate_samtable(SAM.cp73.all_levo_up_and_all_levo_AIM$siggenes.table$genes.up)      

write.table(annotate_samtable(SAM.cp73.all_levo_up_and_all_levo_AIM$siggenes.table$genes.up),
            file="~/Dropbox/CP73_up_allLevo_AIM.xls", quote=FALSE, sep="\t")
write.table(annotate_samtable(SAM.cp101.all_levo_up_and_all_levo_AIM$siggenes.table$genes.up),
            file="~/Dropbox/CP101_up_allLevo_AIM.xls", quote=FALSE, sep="\t")
write.table(annotate_samtable(SAM.cp73.all_levo_down_and_all_levo_AIM$siggenes.table$genes.lo),
            file="~/Dropbox/CP73_down_allLevo_AIM.xls", quote=FALSE, sep="\t")
write.table(annotate_samtable(SAM.cp101.all_levo_down_and_all_levo_AIM$siggenes.table$genes.lo),
            file="~/Dropbox/CP101_down_allLevo_AIM.xls", quote=FALSE, sep="\t")

write.table(annotate_samtable(SAM.cp73.all_levo_vs_saline$siggenes.table$genes.up),
            file="~/Dropbox/CP73_all_levo_up.xls", quote=FALSE, sep="\t")      

write.table(annotate_samtable(SAM.cp101.all_levo_vs_saline$siggenes.table$genes.up),
            file="~/Dropbox/CP101_all_levo_up.xls", quote=FALSE, sep="\t")      

      
# get expression dataset by gene symbols

gene_reps = read.csv(file="~/Dropbox/Projects/Broad/PD_mouse/results/gene_reps.tab", header=FALSE, sep="\t")      
expression.symbols <- merge( exprs(pd.eset), gene_reps, by.x="row.names", by.y="probeset")
expression.symbols <- expression.symbols[,!(names(expression.symbols) %in% c("symbol"))]
rownames(expression.symbols) <- toupper(rownames(expression.symbols))    
      
      
# run GSA

gs.mf <- GSA.read.gmt("~/Dropbox/Projects/Broad/msigdb/c5.mf.v3.0.symbols.gmt")
gs.cgp <- GSA.read.gmt("~/Dropbox/Projects/Broad/msigdb/c2.cgp.v3.1.symbols.gmt")
      
cp101_gsa_mf <- GSA( x = expression.symbols[, c( d_cp101_chronic_high_levo$filenames, d_cp101_chronic_saline$filenames)], 
     y = c( rep(2, dim(d_cp101_chronic_high_levo)[1]), rep(1, dim(d_cp101_chronic_saline)[1])), 
     genesets = gs.mf$genesets,
     genenames = rownames(expression.symbols),
     resp.type = "Two class unpaired",
     nperms = 100)

cp101_gsa_cgp <- GSA( x = expression.symbols[, c( d_cp101_chronic_high_levo$filenames, d_cp101_chronic_saline$filenames)], 
     y = c( rep(2, dim(d_cp101_chronic_high_levo)[1]), rep(1, dim(d_cp101_chronic_saline)[1])), 
     genesets = gs.cgp$genesets,
     genenames = rownames(expression.symbols),
     resp.type = "Two class unpaired",
     nperms = 100)      
      
      
pdf("~/Dropbox/test.pdf")
  for (i in c("1422507_at", "1452731_x_at")) {
       plot_ExprAIM(i, "CP73")      
  }
dev.off()
#mo430genenames[mo430genenames$probe_id %in% rownames(exprs(pd.eset))[as.integer(samfit$siggenes.table$genes.up[,"Gene Name"])]

calc_anova_stats <- function(covar_subset, probeset) {
  cur_expression <- as.data.frame( alldata[probeset,] )
  colnames(cur_expression) <- "expression"
  cur_merged <- merge(covar_subset, cur_expression, by.x="filenames", by.y=0)
  a <- aov(expression ~ DrugTreat, data=cur_merged)
  return(a)
}
      
write_tukeyHSDpvals = function(covar_subset, subset_name) {
  for (p in rownames(alldata)) {
    # get the Tukey stats
    a <- calc_anova_stats(covar_subset, p)
    pv <- as.data.frame(TukeyHSD(a)$DrugTreat[,"p adj"])
    colnames(pv) <- "tukeyHSD"
    pv$probeset <- p
    # append to file
    write.table(pv, file=paste("~/Dropbox/Projects/Broad/PD_mouse/results/jan30/", subset_name, "_tukeyHSD.tab", sep=""), append=TRUE, sep="\t", col.names=FALSE, quote=FALSE)
  }      
}
      
cp73 <- subset(pd.covar, MouseType=="CP73" & LesionType=="6-OHDA")      
cp101 <- subset(pd.covar, MouseType=="CP101" & LesionType=="6-OHDA")      
      
      