require(affy)

loadCELs <- function() {
  cp101files <- list.celfiles("~/Dropbox/Projects/Broad/PD_mouse/CP101", full.names=TRUE)
  cp73files <- list.celfiles("~/Dropbox/Projects/Broad/PD_mouse/CP73", full.names=TRUE)
  
  allfiles <- c(cp73files, cp101files)
  
  cel_data <- affy::ReadAffy(filenames=allfiles)
  pd.eset <- affy::rma(cel_data)
  return(pd.eset)
}

pd.eset <- loadCELs()