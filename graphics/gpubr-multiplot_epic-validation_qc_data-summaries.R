#===================================================
# gpubr multiplot summary validation epic arrays
#===================================================
# Author: Sean Maden

library(ggpubr);library(ggfortify)
bhmepic <- getBeta(gfilt.epic)
bvar <- apply(bhmepic,1,sd); bvar <- bvar[!is.na(bvar)]; bvar <- bvar[rev(order(bvar))]
bpca <- bhmepic[names(bvar[1:1000]),]

datepic <- as.data.frame(pData(gfilt.epic))
epicievora <- datepic$Risk_Group.index; epicievora <- ifelse(epicievora %in% c("2","3"),"red","blue")
table(epicievora)

# PCA 1
xpca <- prcomp(t(bpca))

p1 <- autoplot(xpca,
               data=datepic,
               colour=as.factor(epicievora))
