# ievora EPIC array validation

dvmc.toval <- rownames(dmptcga)


# subset and preprocess epic arrays
#rgepic <- rg.epic.sub[getAnnotation(rg.epic.sub)$Name %in% getAnnotation(rgsub.rl.450)$Name,]
rgepic <- rg.epic.sub[getAnnotation(rg.epic.sub)$Methyl450_Loci==TRUE,]; dim(rgepic)
rgepic <- rgepic[,pData(rgepic)$Anatomic_Location2=="left_colon" & 
                   !(pData(rgepic)$Participant %in% pData(rgsub.rl.450)$Participant)]
rgepic <- rgepic[,!colnames(rgepic)=="NA"]
dim(rgepic);table(rgepic$Risk_Group.index)
summary(as.data.frame(table(rgepic$Participant))[,2])

library(methyPre)
epicmpre <- epicmpre <- methyPre(rgepic,
                                 workflow=c("norm","pfilt","map","minfifilt","crxcgfilt"),
                                 normfun=c("illumina","SWAN"),
                                 crxcgFilt = "epic")
save(epicmpre,"epicmpre.rda")

gfilt.epic <- epicmpre[[6]]; dim(gfilt.epic)

###

#============================
# ROC curve/ AUC evaluation
library(pROC)

gval <- gfilt.epic[rownames(gfilt.epic) %in% rownames(dmptcga),]
roc.epic <- data.frame(cpgid=rownames(gval),
                       auc=rep(NA,nrow(gval)),
                       stringsAsFactors = FALSE)

valvar <- gval$Risk_Group.index; valvar <- ifelse(valvar %in% c("3","2A","2B"),1,0)
table(valvar,gval$Risk_Group.index)

rocvard <- as.numeric(valvar); betaval <- getBeta(gval)
for(i in 1:nrow(roc.epic)){
  bi <- betaval[i,]
  roc.epic$auc[i] <- round(as.numeric(auc(roc(rocvard~bi))),4)
}; roc.epic <- roc.epic[rev(order(roc.epic$auc)),]
