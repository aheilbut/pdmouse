# Drd1 differentially expressed genes by treatment group
venn.plot <- draw.quad.venn(
  area1=226,
  area2=3100,
  area3=4603,
  area4=4773,
  n12=121,
  n13=146,
  n14=95,
  n23=2738,
  n24=1723,
  n34=2508,
  n123=114,
  n124=68,
  n134=81,
  n234=1611,
  n1234=68
  title="blah",
  category = c("dopamine\n depleted (220)", "acute high\nL-DOPA\n(4727)",  "chronic low\nL-DOPA (2980)", "chronic high\nL-DOPA\n(4511)"),
  col = c("yellow", "red", "blue", "purple"),
  fontfamily = "sansserif",
  cat.fontfamily = "sansserif",
  lwd = 8, 
  lty = "solid",
  cex = 5.0,
  alpha = 0.05, 
  cat.cex = 2.8,
  fill = c("yellow", "purple", "red", "blue" ),
#  cat.dist= 0.18
#  cat.pos= 0 
)

# without acute data
venn.plot <- draw.triple.venn(
  area1=226,
  area2=3100,
  area3=4603,
  n12=121,
  n23=2738,
  n13=146,
  n123=114,
  title="blah",
  category = c("dopamine\n depleted (226)", "chronic low\nL-DOPA (3100)", "chronic high\nL-DOPA\n(4603)"),
  col = c("darkorange", "red", "blue"),
  fontfamily = "sansserif",
  cat.fontfamily = "sansserif",
  lwd = 8, 
  lty = "solid",
  cex = 5.0,
  alpha = 0.05, 
  cat.cex = 3.25,
  fill = c("orange", "red", "blue" ),
  cat.dist= 0.12
  #  cat.pos= 0 
)

dev.copy(device=cairo_pdf, file="~/Dropbox/Projects/Broad/PD_mouse/results/2013_10_25/Drd1a_byTreatment.pdf", width=20, height=20 )
dev.off()
dev.off()


# Drd2 differentially expressed genes by treatment group
venn.plot <- draw.triple.venn(
       area1=156,
       area2=72,
       area3=415,
       n12=8,
       n23=53,
       n13=37,
       n123=7,
       category = c("dopamine\ndepleted (156)", "chronic low\n l-DOPA (72)", "chronic high\nl-DOPA (415)"),
       cat.fontfamily = rep("sansserif", 3),
       fontfamily = "sanserif",
       col = c("darkorange", "darkred", "darkblue"),
       lwd = 8,
       lty = "solid",
       cex = 5.0,
       alpha = 0.05,
       cat.cex = 3.25,
       fill = c("darkorange", "darkred", "darkblue"),
       cat.dist = 0.12,
#       margin = 0.05
);


dev.copy(device=cairo_pdf, file="~/Dropbox/Projects/Broad/PD_mouse/results/2013_10_25/Drd2_byTreatment.pdf", width=20, height=20 )
dev.off()
dev.off()
  
# dop depleted, up-regulated
venn.plot <- draw.pairwise.venn(
  area1=79,
  area2=103,
  cross.area=2,
  category = c("Drd2 (iSPNs)\n(103) ", "Drd1a (dSPNs)\n(79)"),
  col = c("darkorange", "orange"),
  fill = c("darkorange", "orange"),
  lty = "solid",
  lwd = 8,
  cex = 4,
  alpha = 0.05,
  fontfamily = "sanserif",
  cat.fontfamily = "sanserif",
  margin = 0.05,
  ext.text = FALSE,
  cat.cex = 4,
  cat.dist=c(0.05,0.05),
  cat.pos=c(0,0),
  inverted=TRUE
);

dev.copy(device=cairo_pdf, file="~/Dropbox/Projects/Broad/PD_mouse/results/2013_10_25/venn_dopDep_UP.pdf", width=16, height=16 )
dev.off()
dev.off()


# dop depleted downregulated
venn.plot <- draw.pairwise.venn(
  area1=147,
  area2=53,
  cross.area=2,
  category=c( "Drd1a (dSPNs)\n(147)", "Drd2 (iSPNs)\n(53)"  ),
  col = c("orange", "darkorange"),
  fill = c("orange", "darkorange"),
  lty = "solid",
  lwd = 8,
  cex = 4,
  alpha = 0.05,
  fontfamily = "sanserif",
  cat.fontfamily = "sanserif",
  margin = 0.05,
  ext.text = FALSE,
  cat.cex = 4,
  cat.dist=c(0.05,0.05),
  cat.pos=c(0,0),  
);

dev.copy(device=cairo_pdf, file="~/Dropbox/Projects/Broad/PD_mouse/results/2013_10_25/venn_dopDep_DOWN.pdf", width=16, height=16 )
dev.off();
dev.off();



# chronic low, up-regulated
venn.plot <- draw.pairwise.venn(
  area1=1352,
  area2=48,
  cross.area=23,
  category=c( "Drd1a (dSPNs)\n(1352)", "Drd2\n(iSPNs)\n  (48)"  ),
  col = c("red", "darkred"),
  fill = c("red", "darkred"),
  lty = "solid",
  lwd = 8,
  cex = 4,
  alpha = 0.05,
  fontfamily = "sanserif",
  cat.fontfamily = "sanserif",
  margin = 0.10,
  ext.text = FALSE,
  cat.cex = 4,
  cat.dist=c(0.06,0.04),
  cat.pos=c(45,45),  
#  cat.pos=c(25, 25),
#  inverted=TRUE
);
dev.copy(device=cairo_pdf, file="~/Dropbox/Projects/Broad/PD_mouse/results/2013_10_25/venn_chronicLow_UP.pdf", width=20, height=16 )
dev.off();
dev.off();



# chronic low, down-regulated
venn.plot <- draw.pairwise.venn(
  area1=1758,
  area2=24,
  cross.area=7,
  category=c( "Drd1a (dSPNs)\n(1758)", "Drd2\n(iSPNs)\n(24)"  ),
  col = c("red", "darkred"),
  fill = c("red", "darkred"),
  lty = "solid",
  lwd = 8,
  cex = 4,
  alpha = 0.05,
  fontfamily = "sanserif",
  cat.fontfamily = "sanserif",
  margin = 0.10,
  ext.text = FALSE,
  cat.cex = 4,
  cat.dist=c(0.05,0.05),
  cat.pos=c(45,45),  
  #  cat.pos=c(25, 25),
  #  inverted=TRUE
);

dev.copy(device=cairo_pdf, file="~/Dropbox/Projects/Broad/PD_mouse/results/2013_10_25/venn_chronicLow_DOWN.pdf", width=20, height=16 )
dev.off();
dev.off();



# chronic High, up-regulated
venn.plot <- draw.pairwise.venn(
  area1=1898,
  area2=244,
  cross.area=102,
  category=c( "Drd1a (dSPNs)\n(1898)", "Drd2 (iSPNs)\n(244)"  ),  
  col = c("blue", "darkblue"),
  fill = c("blue", "darkblue"),
  lty = "solid",
  lwd = 8,
  cex = 4,
  alpha = 0.05,
  fontfamily = "sanserif",
  cat.fontfamily = "sanserif",
  margin = 0.05,
  ext.text = FALSE,
  cat.cex = 4,
  cat.dist=c(0.04,0.04),
  cat.pos=c(20,20),  
  #  cat.pos=c(25, 25),
  #  inverted=TRUE
);

dev.copy(device=cairo_pdf, file="~/Dropbox/Projects/Broad/PD_mouse/results/2013_10_25/figures/venn_chronicHigh_UP.pdf", width=16, height=16 )
dev.off();
dev.off();



# chronic High, up-regulated
venn.plot <- draw.pairwise.venn(
  area1=1871,
  area2=243,
  cross.area=105,
  category=c( "Drd1a (dSPNs)", "Drd2 (iSPNs)"  ),  
  col = c("blue", "darkblue"),
  fill = c("blue", "darkblue"),
  lty = "solid",
  lwd = 8,
  cex = 5,
  alpha = 0.05,
  fontfamily = "sanserif",
  cat.fontfamily = "sanserif",
  margin = 0.05,
  ext.text = FALSE,
  cat.cex = 4,
  cat.dist=c(0.04,0.04),
  cat.pos=c(25,25),  
  #  cat.pos=c(25, 25),
  #  inverted=TRUE
);


# chronic High, down-regulated
venn.plot <- draw.pairwise.venn(
  area1=2733,
  area2=172,
  cross.area=72,
  category=c( "Drd1a (dSPNs)\n(2733)", "Drd2\n(iSPNs)\n(172)"  ),
  
  col = c("blue", "darkblue"),
  fill = c("blue", "darkblue"),
  lty = "solid",
  lwd = 8,
  cex = 5,
  alpha = 0.05,
  fontfamily = "sanserif",
  cat.fontfamily = "sanserif",
  margin = 0.05,
  ext.text = FALSE,
  cat.cex = 4,
  cat.dist=c(0.04,0.05),
  cat.pos=c(25,20),  
  #  cat.pos=c(25, 25),
  #  inverted=TRUE
);


dev.copy(device=cairo_pdf, file="~/Dropbox/Projects/Broad/PD_mouse/results/2013_10_25/figures/venn_chronicHigh_DOWN.pdf", width=16, height=16 )
dev.off();
dev.off();


