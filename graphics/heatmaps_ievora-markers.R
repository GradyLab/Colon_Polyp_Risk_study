# ievora heatmaps

library(ComplexHeatmap)


hmcpg <- roc.epic[roc.epic$auc>=0.7,]$cpgid # prefilter top iEVORA probes (from EPIC validation)

gdisc <- gfilt_ievora; gval <- gfilt.epic


g1 <- gdisc[hmcpg,order(gdisc$Risk_Group.index)]; g1 <- g1[order(match(rownames(g1),hmcpg)),] 
g2 <- gtcgatest[hmcpg,order(gtcgatest$sampletype)]; g2 <- g2[order(match(rownames(g2),hmcpg)),] 
g3 <- gval[hmcpg,order(gval$Risk_Group.index)]; g3 <- g3[order(match(rownames(g3),hmcpg)),] 

identical(hmcpg,rownames(g1)); identical(hmcpg,rownames(g2)); identical(hmcpg,rownames(g3))

drg <- g1$Risk_Group.index; drg[is.na(drg)] <- "0"
hmdisc.col <- HeatmapAnnotation(show_legend = TRUE,
                                df = data.frame(risk_group=drg,location=g1$Anatomic_Location2), 
                                col = list(risk_group = c("0" =  "green","NA"="green","1" = "pink", "2"="red","3"="purple"),
                                           location=c("right_colon"="burlywood1","left_colon"="aquamarine2")), 
                                name = "Risk Group")

tcgaloc <- g2$Anatomic_Loc; tcgaloc <- ifelse(tcgaloc %in% c("Ascending Colon","Cecum","Hepatic Flexure","Transverse Colon"),
                                              "right_colon",ifelse(tcgaloc %in% c("[Discrepancy]","[Not Available]"),"not_available",
                                                                   "left_colon"))
table(tcgaloc,g2$Anatomic_Loc)

hmfilt.col <- HeatmapAnnotation(show_legend = TRUE,
                                df = data.frame(tissue_type=g2$sampletype,location=tcgaloc), 
                                col = list(tissue_type = c("01" =  "purple","11"="green"),
                                           location=c("right_colon"="burlywood1","left_colon"="aquamarine2","not_available"="gray")), 
                                name = "Tissue Type")

hmval.col <- HeatmapAnnotation(show_legend = FALSE,
                                df = data.frame(risk_group=g3$Risk_Group.index,location=g3$Anatomic_Location2), 
                                col = list(risk_group = c("0" =  "green","NA"="green","1" = "pink", "2A"="red","2B"="red","3"="purple"),
                                           location=c("left_colon"="aquamarine2")), 
                                name = "Risk Group")
library(circlize)
breaks=seq(0,1,0.05)
hmcol = colorRamp2(breaks,colorRampPalette(c("darkblue","yellow"))(n=length(breaks)))

hmdisc <- Heatmap(getBeta(g1),col=hmcol,
                  top_annotation=hmdisc.col,
                  heatmap_legend_param = list(color_bar="continuous"),
                  name="DVMC Methylation\n(Beta_value)",
                  row_title=paste0("Validated DVMC (CpG ID, N=",nrow(g1),")"),
                  column_title=paste0("HM450K Discovery\nSamples (N=",ncol(g1),")"),
                  show_row_names=FALSE,
                  show_column_names=FALSE,
                  show_column_dend = FALSE,
                  show_heatmap_legend = TRUE,
                  show_row_dend = FALSE)

hmval <- Heatmap(getBeta(g3),col=hmcol,
                 top_annotation=hmval.col,
                 row_title="Validated DVMC (CpG ID)",
                 column_title=paste0("EPIC Validation\nSamples (N=",ncol(g3),")"),
                 show_row_names=FALSE,
                 show_column_names=FALSE,
                 show_column_dend=FALSE,
                 show_heatmap_legend = FALSE,
                 show_row_dend = FALSE)

hmfilt <- Heatmap(getBeta(g2),col=hmcol,
                  top_annotation=hmfilt.col,
                  row_title=paste0("Validated DVMC (CpG ID, N=",nrow(g2),")"),
                  column_title=paste0("TCGA Filter Samples (N=",ncol(g2),")"),
                  name="DVMC Methylation\n(Beta_value)",
                  show_row_names=FALSE,
                  show_column_names=FALSE,
                  show_column_dend = FALSE,
                  show_heatmap_legend = TRUE,
                  show_row_dend = FALSE)


jpeg("HM-ievora-valDVMC.jpg",10,5,units="in",res=400)
draw(hmdisc+hmval)
dev.off()

jpeg("HM-ievora-tcgaFilt-valDVMC.jpg",7,4,units="in",res=400)
draw(hmfilt)
dev.off()
