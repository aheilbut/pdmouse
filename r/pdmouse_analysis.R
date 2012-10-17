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
pd.covar$totalAIM <- apply(pd.covar[c("Day1AIM","Day3AIM","Day4AIM","Day5AIM","Day8AIM")], 1, sum)


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
d_chronic_saline <- subset(pd.covar, MouseType=="CP73" & DrugTreat=="Chronic saline" & LesionType=="6-OHDA")
d_chronic_low_levo <- subset(pd.covar, MouseType=="CP73" & DrugTreat=="Chronic low levodopa" & LesionType=="6-OHDA")
d_chronic_high_levo <- subset(pd.covar, MouseType=="CP73" & DrugTreat=="Chronic high levodopa" & LesionType=="6-OHDA")
d_acute_high_levo <- subset(pd.covar, MouseType=="CP73" & DrugTreat=="Acute high levodopa" & LesionType=="6-OHDA")
d_acute_saline <- subset(pd.covar, MouseType=="CP73" & DrugTreat=="Acute saline" & LesionType=="6-OHDA")


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
# y <- c(rep(0, dim(a)[1]), d_c$totalAIM)
y <- c(d_b$totalAIM, d_c$totalAIM)

data <- exprs(pd.eset[,d_chronic_high_levo$filenames])
y <- c(d_chronic_high_levo$totalAIM)

SAM.chronic_high_vs_AIM <- SAM(exprs(pd.eset[,d_chronic_high_levo$filenames]),
                d_chronic_high_levo$totalAIM,
                resp.type="Quantitative", 
                geneid=row.names(data))


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
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      

#mo430genenames[mo430genenames$probe_id %in% rownames(exprs(pd.eset))[as.integer(samfit$siggenes.table$genes.up[,"Gene Name"])]
