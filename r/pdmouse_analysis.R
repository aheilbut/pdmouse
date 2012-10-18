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

SAM.cp73.chronic_high_vs_AIM <- SAM(exprs(pd.eset[,d_cp73_chronic_high_levo$filenames]),
                d_cp73_chronic_high_levo$totalAIM,
                resp.type="Quantitative", 
                geneid=row.names(data),  
                fdr.output=fdr_cutoff)

SAM.cp73.all_levo_vs_AIM <- SAM(exprs(pd.eset[,d_cp73_all_levo$filenames]),
                                 d_cp73_all_levo$totalAIM,
                                 resp.type="Quantitative", 
                                 geneid=row.names(data),
                                fdr.output=fdr_cutoff)

SAM.cp101.chronic_high_vs_AIM <- SAM(exprs(pd.eset[,d_cp101_chronic_high_levo$filenames]),
                                     d_cp101_chronic_high_levo$totalAIM,
                                     resp.type="Quantitative", 
                                     geneid=row.names(data),
                                     fdr.output=fdr_cutoff)


SAM.cp101.chronic_high_vs_AIMDay8 <- SAM(exprs(pd.eset[,d_cp101_chronic_high_levo$filenames]),
                                     d_cp101_chronic_high_levo$Day8AIM,
                                     resp.type="Quantitative", 
                                     geneid=row.names(data),
                                     fdr.output=fdr_cutoff)



SAM.cp101.chronic_low_vs_AIM <- SAM(exprs(pd.eset[,d_cp101_chronic_low_levo$filenames]),
                                     d_cp101_chronic_low_levo$totalAIM,
                                     resp.type="Quantitative", 
                                     geneid=row.names(data),
                                     fdr.output=fdr_cutoff)

SAM.cp73.chronic_low_vs_AIM <- SAM(exprs(pd.eset[,d_cp73_chronic_low_levo$filenames]),
                                     d_cp73_chronic_low_levo$totalAIM,
                                     resp.type="Quantitative", 
                                     geneid=row.names(data),
                                     fdr.output=fdr_cutoff)


SAM.cp101.all_levo_vs_AIM <- SAM(exprs(pd.eset[,d_cp101_all_levo$filenames]),
                                 d_cp101_all_levo$totalAIM,
                                 resp.type="Quantitative", 
                                 geneid=row.names(data),
                                 fdr.output=fdr_cutoff)


SAM.cp101.all_levo_vs_saline <- SAM(exprs(pd.eset[,c(d_cp101_all_levo$filenames,d_cp101_all_saline$filenames)]),
                                 c( rep(2, dim(d_cp101_all_levo)[1]), rep(1, dim(d_cp101_all_saline)[1])),
                                 resp.type="Two class unpaired", 
                                 geneid=row.names(data),
                                 fdr.output=fdr_cutoff)

SAM.cp73.all_levo_vs_saline <- SAM(exprs(pd.eset[,c(d_cp73_all_levo$filenames,d_cp73_all_saline$filenames)]),
                                    c( rep(2, dim(d_cp73_all_levo)[1]), rep(1, dim(d_cp73_all_saline)[1])),
                                    resp.type="Two class unpaired", 
                                    geneid=row.names(data),
                                    fdr.output=fdr_cutoff)




SAM.cp101.chronic_high_vs_saline <- SAM(exprs(pd.eset[,c(d_cp101_chronic_high_levo$filenames, d_cp101_chronic_saline$filenames)]),
                                        c(rep(2, dim(d_cp101_chronic_high_levo)[1]), rep(1, dim(d_cp101_chronic_saline)[1])),
                                        resp.type="Two class unpaired", 
                                        geneid=row.names(data),
                                        fdr.output=0.20)  

SAM.cp73.chronic_high_vs_saline <- SAM(exprs(pd.eset[,c(d_cp73_chronic_high_levo$filenames, d_cp73_chronic_saline$filenames)]),
                                        c(rep(2, dim(d_cp73_chronic_high_levo)[1]), rep(1, dim(d_cp73_chronic_saline)[1])),
                                        resp.type="Two class unpaired", 
                                        geneid=row.names(data),
                                        fdr.output=0.20)


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
                                            fdr.output=0.05)


responsive_set <- SAM.cp101.all_levo_vs_saline$siggenes.table$genes.up[as.numeric(as.data.frame(SAM.cp101.all_levo_vs_saline$siggenes.table$genes.up, stringsAsFactors=FALSE )$'q-value(%)') < 10, "Gene Name"  ]
SAM.cp101.all_levo_up_and_all_levo_AIM <- SAM(exprs(pd.eset[responsive_set,d_cp101_all_levo$filenames]),
                                             d_cp101_all_levo$totalAIM,
                                             resp.type="Quantitative", 
                                             geneid=responsive_set,
                                             fdr.output=0.05)



# plot gene VS. AIM scores for 


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
               sort=FALSE), mo430symbol, by.x="Gene Name", by.y="probe_id", sort=FALSE)[1:min(100, dim(sam_result)[1]),] )
  
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
       size=I(2), 
       xlab="log2 expression", 
        ylab="integral of AIM score",
       main=title, 
       geom="jitter"
        ) 
  
  return(p)
 }      
      

#for 
      
printGraphs <- function(samset, filename, description, maxgraphs) {       
  pdf(file=paste("~/Dropbox/", filename))
  current <- samset
  #grid.table(annotate_samtable(current$siggenes.table$genes.up), gp=gpar(fontsize=6))      
  current.up <- as.data.frame(samset$siggenes.table$genes.up, stringsAsFactors=FALSE)
  current.lo <- as.data.frame(samset$siggenes.table$genes.lo, stringsAsFactors=FALSE )
  
  if (!is.null(current$siggenes.table$genes.up)) {
  for (g in rownames(current.up[1:min(maxgraphs, dim(current.up)[1]),])) { 
    probe <- current.up[g, "Gene Name"]
    #do.call(arrange_ggplot2, list(do.call(plot_ExprAIM, list(probe, "CP73")), do.call(plot_ExprAIM, list(probe, "CP101")), ncol=1))
    arrange_ggplot2( plot_ExprAIM(probe, "CP73"),  plot_ExprAIM(probe, "CP101"), ncol=1)
  }
  }
      
  if (!is.null(current$siggenes.table$genes.lo)) {
  for (g in rownames(current.lo[1:min(maxgraphs, dim(current.lo)[1]),])) { 
    probe <- current.lo[g, "Gene Name"]
#    do.call( arrange_ggplot2, list( do.call(plot_ExprAIM, list(probe, "CP73")), do.call(plot_ExprAIM, list(probe, "CP101")), ncol=1))
    arrange_ggplot2( plot_ExprAIM(probe, "CP73"),  plot_ExprAIM(probe, "CP101"), ncol=1)
    
#    plot_ExprAIM(probe, "CP73")
#    plot_ExprAIM(probe, "CP101")
    
  }
  }
  dev.off()      
}
   

printGraphs(SAM.cp101.chronhighup_chronhighlevo_vs_AIM, "sam_cp101_chron_high_levo_vs_saline_VS_AIM.pdf", "Filter genes up vs saline, then scored by regression to AIM", 10)
      
printGraphs(SAM.cp73.chronhighup_chronhighlevo_vs_AIM, "sam_cp73_chron_high_levo_vs_saline_VS_AIM.pdf", "Filter genes up vs saline, then scored by regression to AIM", 10)

printGraphs(SAM.cp73.chronhighup_chronhighlevo_vs_AIM, "sam_cp73_chron_high_levo_vs_saline_VS_AIM.pdf", "Filter genes up vs saline, then scored by regression to AIM", 10)
      
printGraphs(SAM.cp73.all_levo_up_and_all_levo_AIM, "SAM.cp73.all_levo_up_and_all_levo_AIM.pdf", "Filter genes up vs saline, then scored by regression to AIM", 10)
printGraphs(SAM.cp101.all_levo_up_and_all_levo_AIM, "SAM.cp101.all_levo_up_and_all_levo_AIM.pdf", "Filter genes up vs saline, then scored by regression to AIM", 10)
      
pdf("~/Dropbox/test.pdf")
  for (i in c("1422507_at", "1452731_x_at")) {
       plot_ExprAIM(i, "CP73")      
  }
dev.off()
#mo430genenames[mo430genenames$probe_id %in% rownames(exprs(pd.eset))[as.integer(samfit$siggenes.table$genes.up[,"Gene Name"])]
