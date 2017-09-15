#===============================
# ievora EPIC array validation
#===============================
# author: Sean Maden

#=========================================
# Validate post-filter discovery markers
#=========================================

gval <- gfilt.epic[rownames(gfilt.epic) %in% rownames(dmptcga),] # all left colon


# recategorize NCCN risk group to be LR = 0/HR = 1
valvar <- gval$Risk_Group.index; valvar <- ifelse(valvar %in% c("3","2A","2B"),1,0)

table(valvar,gval$Risk_Group.index)
# valvar  0  1 2A 2B  3
# 0  7  9  0  0  0
# 1  0  0  1  1 27

# prefilter on methylation differences (LR<HR)
bepic <- getBeta(gval); 
blr <- rowMeans(bepic[,which(valvar==0)]); bhr <- rowMeans(bepic[,which(valvar==1)])
bepic <- bepic[which(blr<bhr),]; dim(bepic)

# ROC curve/ AUC evaluation
gval <- gval[rownames(gval) %in% rownames(bepic),]

library(pROC)

roc.epic <- data.frame(cpgid=rownames(gval),
                       auc=rep(NA,nrow(gval)),
                       stringsAsFactors = FALSE)

rocvard <- as.numeric(valvar); betaval <- getBeta(gval)
for(i in 1:nrow(roc.epic)){
  bi <- betaval[i,]
  roc.epic$auc[i] <- round(as.numeric(auc(roc(rocvard~bi))),4)
}; roc.epic <- roc.epic[rev(order(roc.epic$auc)),]

valdvmc <- roc.epic[roc.epic$auc>=0.7,]$cpgid # N = 122 CpGs

valdvmc

###
