#==============================================
# gpubr multiplot, HM450K discovery set data
#=============================================
# Author: Sean Maden

library(ggpubr);library(ggfortify)


bhm450 <- getBeta(gfilt_ievora)
bvar <- apply(bhm450,1,sd); bvar <- bvar[!is.na(bvar)]; bvar <- bvar[rev(order(bvar))]
bpca <- bhm450[names(bvar),]

datrl <- as.data.frame(pData(gfilt_ievora))
rgievora <- datrl$Risk_Group.index; rgievora <- ifelse(rgievora %in% c("2","3"),"red","blue")
table(rgievora)

# PCA 1
xpca <- prcomp(t(bpca))

p1 <- autoplot(xpca,
               data=datrl,
               colour=as.factor(rgievora))

# PCA 2
col2 <- datrl$Run_Date; 
clist2 <- sample(colours(),length(unique(col2)));names(clist2) <- unique(col2)
for(i in 1:length(unique(col2))){col2[col2==names(clist2)[i]] <- clist2[i]}

p2 <- autoplot(xpca,
               data=datrl,
               colour=col2,legend=TRUE)+ theme(legend.position = c(0.8, 0.2)); p2

# PCA 3
col3 <- datrl$Anatomic_Location2; 
clist3 <- sample(colours(),length(unique(col3)));names(clist3) <- unique(col3)
for(i in 1:length(unique(col3))){col3[col3==names(clist3)[i]] <- clist3[i]}

p3 <- autoplot(xpca,
               data=datrl,
               colour=col3,legend=TRUE)+ theme(legend.position = c(0.8, 0.2)); p3

# PCA 4
col4 <- col3; 
clist4 <- clist3

bdvc <- getBeta(gfilt_ievora[rldvc.list,])
dvcvar <- apply(bdvc,1,sd); dvcvar <- dvcvar[!is.na(dvcvar)]; dvcvar <- dvcvar[rev(order(dvcvar))]
bpca.dvc <- bdvc[names(dvcvar[1:1000]),]; summary(apply(bpca.dvc,1,sd))

p4 <- autoplot(prcomp(t(bpca.dvc)),
               data=datrl,
               colour=col4,legend=TRUE)+ theme(legend.position = c(0.8, 0.2)); p4

# Age summary table
stable.age <- matrix(c("right_colon",
                       "left_colon",
                       "all_colon",
                       nrow(datrl[datrl$Anatomic_Location2=="right_colon",]),
                       nrow(datrl[datrl$Anatomic_Location2=="left_colon",]),
                       nrow(datrl),
                       round(mean(datrl[datrl$Anatomic_Location2=="right_colon",]$Age),0),
                       round(mean(datrl[datrl$Anatomic_Location2=="left_colon",]$Age),0),
                       round(mean(datrl$Age),0),
                       round(sd(datrl[datrl$Anatomic_Location2=="right_colon",]$Age),0),
                       round(sd(datrl[datrl$Anatomic_Location2=="left_colon",]$Age),0),
                       round(sd(datrl$Age),0)),
                     ncol=4,nrow=3)
colnames(stable.age) <- c("anatomic_location","n_samples","mean_age","sd_age")
stable.age <- ggtexttable(stable.age,rows=NULL,theme=ttheme("mBlue"))
stable.age

# Anatomic Locations Table
stable.loc <- ggtexttable(t(table(datrl$Anatomic_Location2,datrl$Anatomic_Location1)),
                          theme=ttheme("mBlue",base_size=8))

# stable gender
stable.gender <- matrix(c("right_colon",
                          "left_colon",
                          "all_colon",
                          nrow(datrl[datrl$Anatomic_Location2=="right_colon",]),
                          nrow(datrl[datrl$Anatomic_Location2=="left_colon",]),
                          nrow(datrl),
                          nrow(datrl[datrl$Anatomic_Location2=="right_colon" & datrl$Gender=="M",]),
                          nrow(datrl[datrl$Anatomic_Location2=="left_colon" & datrl$Gender=="M",]),
                          nrow(datrl[datrl$Gender=="M",]),
                          nrow(datrl[datrl$Anatomic_Location2=="right_colon" & datrl$Gender=="F",]),
                          nrow(datrl[datrl$Anatomic_Location2=="left_colon" & datrl$Gender=="F",]),
                          nrow(datrl[datrl$Gender=="F",]),
                          100*round(nrow(datrl[datrl$Anatomic_Location2=="right_colon" & datrl$Gender=="F",])/nrow(datrl[datrl$Anatomic_Location2=="right_colon",]),2),
                          100*round(nrow(datrl[datrl$Anatomic_Location2=="left_colon" & datrl$Gender=="F",])/nrow(datrl[datrl$Anatomic_Location2=="left_colon",]),2),
                          100*round(nrow(datrl[datrl$Gender=="F",])/nrow(datrl),2)),
                        ncol=5,nrow=3)
colnames(stable.gender) <- c("anatomic_location","n_samples","n_males","n_females","female (%)")
stable.gender <- ggtexttable(t(stable.gender),theme=ttheme("mBlue",base_size = 8))
stable.gender


textt <- paste("HM450K Discovery Set Summary",
               "\n(left) PCA Analyses (1K Most Variable CpGs) colorized by (top to bottom):",
               "\n1. iEVORA Risk Group, LR = blue / HR = red, all CpGs",
               "\n2. Run Date/Batch, at all CpGs",
               "\n3. Colon Location at all CpGs,",
               "\n4. Colon Location, at most variable DVCs",
               "\n(right) Summary tables of (top to bottom):",
               "1.Age,\n2.Gender,\n3.Colon Location",
               sep=" ")
text.p <- ggparagraph(text = textt, size = 11, color = "black")

arrange1 <- ggarrange(p1,p2,p3,p4,ncol=1,nrow=4)
arrange2 <- ggarrange(stable.age,stable.gender,stable.loc,text.p,ncol=1,nrow=4)

jpeg("hm450disc-multiplot-preprocess-summaries.jpg",10,15,units="in",res=400)
ggarrange(arrange1,arrange2,ncol = 2, nrow = 1)
dev.off()
