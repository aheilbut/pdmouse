# Drd1 differentially expressed genes by treatment group
venn.plot <- draw.quad.venn(
  area1=220,
  area2=2980,
  area3=4511,
  area4=4727,
  n12=123,
  n13=142,
  n14=100,
  n23=2610,
  n24=1654,
  n34=2433,
  n123=113,
  n124=72,
  n134=86,
  n234=1540,
  n1234=71,
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



# Drd2 differentially expressed genes by treatment group
venn.plot <- draw.triple.venn(
       area1=154,
       area2=76,
       area3=408,
       n12=9,
       n23=53,
       n13=38,
       n123=7,
       category = c("dopamine\ndepleted (154)", "chronic low\n l-dopa (76)", "chronic high\nl-dopa (408)"),
       cat.fontfamily = rep("sansserif", 3),
       fontfamily = "sanserif",
       col = c("yellow", "red", "blue"),
       lwd = 8,
       lty = "solid",
       cex = 5.0,
       alpha = 0.05,
       cat.cex = 2.8,
       fill = c("yellow", "red", "blue"),
       cat.dist = 0.15,
       margin = 0.05
);
  
# dop depleted, up-regulated
venn.plot <- draw.pairwise.venn(
  area1=70,
  area2=103,
  cross.area=2,
  category = c("Drd2 (iSPNs)\n(103) ", "Drd1a (dSPNs)\n(70)"),
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

dev.copy(device=cairo_pdf, file="~/Dropbox/Projects/Broad/PD_mouse/Manuscript/figures_0929/venn_dopDep_UP.pdf", width=16, height=16 )
dev.off()
dev.off()


# dop depleted downregulated
venn.plot <- draw.pairwise.venn(
  area1=150,
  area2=51,
  cross.area=2,
  category=c( "Drd1a (dSPNs)\n(150)", "Drd2 (iSPNs)\n(51)"  ),
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

dev.copy(device=cairo_pdf, file="~/Dropbox/Projects/Broad/PD_mouse/Manuscript/figures_0929/venn_dopDep_DOWN.pdf", width=16, height=16 )
dev.off();
dev.off();



# chronic low, up-regulated
venn.plot <- draw.pairwise.venn(
  area1=1307,
  area2=50,
  cross.area=24,
  category=c( "Drd1a (dSPNs)\n(1307)", "Drd2\n(iSPNs)\n  (50)"  ),
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
dev.copy(device=cairo_pdf, file="~/Dropbox/Projects/Broad/PD_mouse/Manuscript/figures_0929/venn_chronicLow_UP.pdf", width=20, height=16 )
dev.off();
dev.off();



# chronic low, down-regulated
venn.plot <- draw.pairwise.venn(
  area1=1683,
  area2=26,
  cross.area=9,
  category=c( "Drd1a (dSPNs)\n(1683)", "Drd2\n(iSPNs)\n(26)"  ),
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

dev.copy(device=cairo_pdf, file="~/Dropbox/Projects/Broad/PD_mouse/Manuscript/figures_0929/venn_chronicLow_DOWN.pdf", width=20, height=16 )
dev.off();
dev.off();



# chronic High, up-regulated
venn.plot <- draw.pairwise.venn(
  area1=1871,
  area2=243,
  cross.area=105,
  category=c( "Drd1a (dSPNs)\n(1871)", "Drd2 (iSPNs)\n(243)"  ),  
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

dev.copy(device=cairo_pdf, file="~/Dropbox/Projects/Broad/PD_mouse/Manuscript/figures_0929/venn_chronicHigh_UP.pdf", width=16, height=16 )
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
  area1=2668,
  area2=166,
  cross.area=67,
  category=c( "Drd1a (dSPNs)\n(2668)", "Drd2\n(iSPNs)\n(166)"  ),
  
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


dev.copy(device=cairo_pdf, file="~/Dropbox/Projects/Broad/PD_mouse/Manuscript/figures_0929/venn_chronicHigh_DOWN.pdf", width=16, height=16 )
dev.off();
dev.off();


