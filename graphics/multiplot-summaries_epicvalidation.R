#===================================================
# multiplot summaries epic validation arrays
#===================================================
# Author: Sean Maden

library(ggpubr);library(ggfortify)
bhmepic <- getBeta(gfilt.epic)
bvar <- apply(bhmepic,1,sd); bvar <- bvar[!is.na(bvar)]; bvar <- bvar[rev(order(bvar))]
bpca <- bhmepic[names(bvar[1:1000]),]

datepic <- as.data.frame(pData(gfilt.epic))
epicievora <- datepic$Risk_Group.index; epicievora <- ifelse(epicievora %in% c("2","2A","2B","3"),"red","blue")
table(epicievora)

# PCA 1
xpca <- prcomp(t(bpca))

p1 <- autoplot(xpca,
               data=datepic,
               colour=as.factor(epicievora))


# Age Summary Table
stable.age <- matrix(c("LR",
                       "HR",
                       "all",
                       length(epicievora[epicievora=="blue"]),
                       length(epicievora[epicievora=="red"]),
                       length(epicievora),
                       round(mean(datepic[which(epicievora=="blue"),]$Age),0),
                       round(mean(datepic[which(epicievora=="red"),]$Age),0),
                       round(mean(datepic$Age),0),
                       round(sd(datepic[which(epicievora=="blue"),]$Age),0),
                       round(sd(datepic[which(epicievora=="red"),]$Age),0),
                       round(sd(datepic$Age),0)),
                     ncol=4,nrow=3)
colnames(stable.age) <- c("risk_group","n_samples","mean_age","sd_age")
stable.age <- ggtexttable(stable.age,rows=NULL,theme=ttheme("mBlue"))
stable.age

# Gender Summary Table
stable.gender <- matrix(c("LR",
                       "HR",
                       "all",
                       nrow(datepic[datepic$Risk_Group.index %in% c("0","1"),]),
                       nrow(datepic[datepic$Risk_Group.index %in% c("2A","2B","3"),]),
                       length(epicievora),
                       nrow(datepic[datepic$Risk_Group.index %in% c("0","1") & datepic$Gender=="M",]),
                       nrow(datepic[datepic$Risk_Group.index %in% c("2A","2B","3") & datepic$Gender=="M",]),
                       nrow(datepic[datepic$Gender=="M",]),
                       nrow(datepic[datepic$Risk_Group.index %in% c("0","1") & datepic$Gender=="F",]),
                       nrow(datepic[datepic$Risk_Group.index %in% c("2A","2B","3") & datepic$Gender=="F",]),
                       nrow(datepic[datepic$Gender=="F",]),
                       100*round(nrow(datepic[datepic$Risk_Group.index %in% c("0","1") & datepic$Gender=="F",])/nrow(datepic[datepic$Risk_Group.index %in% c("0","1"),]),2),
                       100*round(nrow(datepic[datepic$Risk_Group.index %in% c("2A","2B","3") & datepic$Gender=="F",])/nrow(datepic[datepic$Risk_Group.index %in% c("2A","2B","3"),]),2),
                       100*round(nrow(datepic[datepic$Gender=="F",])/nrow(datepic),2)),
                     ncol=5,nrow=3)
colnames(stable.gender) <- c("risk_group","n_samples","n_male","n_female","female(%)")
stable.gender <- ggtexttable(stable.gender,rows=NULL,theme=ttheme("mBlue"))
stable.gender

textt <- paste("EPIC validation set summary",
               "\n(left) PCA Analyses (1K Most Variable CpGs) colored on Risk group (LR = blue/HR = red)",
               "\n(right) Summary tables for, from top to bottom, 1. Age, 2. Gender, by risk group", sep=" ")
text.p <- ggparagraph(text = textt, size = 11, color = "black")

arrange1 <- ggarrange(p1,ncol=1,nrow=1)
arrange2 <- ggarrange(stable.age,stable.gender,text.p,ncol=1,nrow=3)

jpeg("epicval-multiplot-preprocess-summaries.jpg",10,15,units="in",res=400)
ggarrange(arrange1,arrange2,ncol = 2, nrow = 1)
dev.off()
