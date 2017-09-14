# ievora graphics

#



#==================================
# four markers validated in epic
valmarkers <- intersect(roc.epic$cpgid,old.markers)


# preamble = colors
gx <- gdisc; gx <- gx[,order(gx$Risk_Group.index)] 
rgid <- gx$Risk_Group.index; rgid[is.na(rgid)] <- "0"
dv <- rgid; dv <- ifelse(dv==3,"purple",ifelse(dv %in% c("2","2A","2B"),"red",ifelse(dv=="1","pink","green")))
table(dv,gx$Risk_Group.index)

gx <- gtcgatest[markeri,order(gtcgatest$sampletype)]; 
tcgatissues <- gx$sampletype; tcgatissues <- ifelse(tcgatissues=="01","T","NM")
tcgacol <- gx$sampletype; tcgacol <- ifelse(tcgacol=="01","blue","green")

gx <- gval;
gx <- gx[markeri,]; gx <- gx[,order(gx$Risk_Group.index)] 
rgiv <- gx$Risk_Group.index; rgiv[is.na(rgiv)] <- "0"
predvar <- rgiv; predvar<-ifelse(predvar %in% c("0","1"),0,1) 
dv <- rgiv; dv <- ifelse(dv==3,"purple",ifelse(dv %in% c("2","2A","2B"),"red",ifelse(dv=="1","pink","green")))
table(dv,gx$Risk_Group.index)


jpeg("plot-summary_ievora-epicvalmarkers.jpg",12,12,units="in",res=400)
par(mfrow=c(4,4),oma=c(0,12,0,0))

bpcex1 = 1.5; bpcex2 = 1.5

for(i in 1:4){
  markeri <- valmarkers[i]
  gx <- gdisc; gx <- gx[,order(gx$Risk_Group.index)]
  titlex <- paste0(markeri,"; HM450 Discovery Samples")
  colgrpx <- c("green","pink","red","purple")
  
  if(i==1){
    # boxplot 1
    boxplot(as.numeric(getBeta(gx[markeri,]))~rgid,col=colgrpx,
            ylim=c(0,1),ylab="Beta-value",cex.axis=bpcex1,cex.lab=bpcex2)
    mtext("HM450\nDiscovery",las=1,side=2,line=2,at=grconvertY(0.5,"npc","nic"),outer=TRUE)
    mtext(paste0(markeri),las=1,side=3,line=-3,at=grconvertX(0.5,"npc","nic"),outer=TRUE)
  } else{
    # boxplot 1
    boxplot(as.numeric(getBeta(gx[markeri,]))~rgid,col=colgrpx,
            ylim=c(0,1),cex.axis=bpcex1,cex.lab=bpcex2)
    mtext(paste0(markeri),las=1,side=3,line=-3,at=grconvertX(0.5,"npc","nic"),outer=TRUE)
  }
  
}

for(i in 1:4){
  markeri <- valmarkers[i]
  gx <- gval; gx <- gx[,order(gx$Risk_Group.index)]
  colgrpx <- c("green","pink","red","red","purple")
  # boxplot 2
  
  if(i==1){
    boxplot(as.numeric(getBeta(gx[markeri,]))~rgiv,col=colgrpx,
            ylim=c(0,1),ylab="Beta-value",cex.axis=bpcex1,cex.lab=bpcex2)
    mtext("EPIC\nValidation",las=1,side=2,line=2,at=grconvertY(0.5,"npc","nic"),outer=TRUE)
  } else{
    boxplot(as.numeric(getBeta(gx[markeri,]))~rgiv,col=colgrpx,
            ylim=c(0,1),cex.axis=bpcex1,cex.lab=bpcex2)
  }
  
}
  
for(i in 1:4){
  markeri <- valmarkers[i]
  gx <- gtcgatest[markeri,order(gtcgatest$sampletype)];
  colgrpx <- c("blue","green")
  # boxplot 3
  
  if(i==1){
    boxplot(as.numeric(getBeta(gx))~tcgatissues,col=colgrpx,
            ylim=c(0,1),ylab="Beta-value",cex.axis=bpcex1,cex.lab=bpcex2)
    mtext("TCGA Filter",las=1,side=2,line=2,at=grconvertY(0.5,"npc","nic"),outer=TRUE)
  } else{
    boxplot(as.numeric(getBeta(gx))~tcgatissues,col=colgrpx,
            ylim=c(0,1),cex.axis=bpcex1,cex.lab=bpcex2)
  }
  
  
}
  
for(i in 1:4){
  markeri <- valmarkers[i]
  # roc curve
  gx <- gval; gx <- gx[,order(gx$Risk_Group.index)]
  if(i==1){
    plot.roc(roc(predvar~as.numeric(getBeta(gx[markeri,]))),
             auc.polygon=TRUE,
             print.auc=TRUE,
             legacy.axes = TRUE,
             print.auc.cex=1.5,cex.axis=bpcex1,cex.lab=bpcex2)
    mtext("EPIC Samples\nROC Curve",las=1,side=2,line=2,at=grconvertY(0.5,"npc","nic"),outer=TRUE)
  } else{
    plot.roc(roc(predvar~as.numeric(getBeta(gx[markeri,]))),
             auc.polygon=TRUE,
             print.auc=TRUE,
             legacy.axes = TRUE,
             print.auc.cex=1.5,cex.axis=bpcex1,cex.lab=bpcex2)
  }
  
} 

dev.off()
