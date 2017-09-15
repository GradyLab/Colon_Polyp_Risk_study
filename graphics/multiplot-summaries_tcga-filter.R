#====================================
# gtcga multiplot, iEVORA analysis
#==================================
# Author: Sean Maden

library(ggpubr);library(ggfortify)

#gtcgatest <- gtcga[,colnames(gtcga) %in% colnames(gtcgatest)]
#save(gtcgatest,file="gfilt-tcga_ievora_allcpg.rda")

btcga <- getBeta(gtcgatest)
bvar <- apply(btcga,1,sd); btcga <- btcga[rev(order(bvar)),]
bpca <- btcga[1:1000,]

xpca <- prcomp(t(bpca))

ggplot(xpca,aes(x=apply(xpca,2,var)))+geom_hist()

p1 <- autoplot(prcomp(t(bpca)),
         data=as.data.frame(pData(gtcgatest)),
         colour=as.factor(gtcgatest$sampletype))

col2 <- substr(gtcgatest$TCGA_Barcode,22,25); 
clist2 <- sample(colours(),length(unique(col2)));names(clist2) <- unique(col2)
for(i in 1:length(unique(col2))){col2[col2==names(clist2)[i]] <- clist2[i]}

p2 <- autoplot(prcomp(t(bpca)),
               data=as.data.frame(pData(gtcgatest)),
               colour=col2,legend=TRUE)+ theme(legend.position = c(0.8, 0.2)); p2

# Table of descriptive stats and data summary
xtcga <- as.data.frame(pData(gtcgatest)); xtcga$Age <- as.numeric(xtcga$Age)

stable <- desc_statby(xtcga, measure.var = "Age",
                      grps = "sampletype")
stable <- stable[, c("sampletype", "length", "mean", "sd")]
stable.age <- ggtexttable(stable, rows = NULL, 
                        theme = ttheme("mBlue"))

stable.gender <- matrix(c("01","11",
                          nrow(xtcga[xtcga$sampletype=="01" & xtcga$Gender=="MALE",]),
                          nrow(xtcga[xtcga$sampletype=="11" & xtcga$Gender=="MALE",]),
                          nrow(xtcga[xtcga$sampletype=="01" & xtcga$Gender=="FEMALE",]),
                          nrow(xtcga[xtcga$sampletype=="11" & xtcga$Gender=="FEMALE",])),ncol=3)
colnames(stable.gender) <- c("sampletype","male","female")
stable.gender <- ggtexttable(stable.gender,rows=NULL,theme=ttheme("mBlue"))


textt <- paste("TCGA Dataset Preprocessing and summary.",
               "(topleft) PCA of top 1k most variable probes,",
               "samples colored by tissue (black=T / green=NM).",
               "(bottomleft) PCA of top 1k most variable probes, samples colored by plate.",
               "(top- and midright) Descriptive statistics summary tables for Age and Gender",sep=" ")
text.p <- ggparagraph(text = textt, size = 11, color = "black")

arrange1 <- ggarrange(p1,p2,ncol=1,nrow=2)
arrange2 <- ggarrange(stable.age,stable.gender,text.p,ncol=1,nrow=3)
ggarrange(arrange1,arrange2,ncol = 2, nrow = 1)
