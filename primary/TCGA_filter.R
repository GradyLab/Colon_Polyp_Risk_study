# ievora tcga filters

table(substr(gtcga$TCGA_Barcode,22,25)) # batches/plates

mdsPlot(getBeta(gtcga),
        sampGroups=substr(gtcga$TCGA_Barcode,22,25),
        sampNames=gtcga$sampletype) 

# exclude matched tumors
gtcgatest <- gtcga[,!(gtcga$Patient_ID %in% gtcga[,gtcga$sampletype=="11"]$Patient_ID &
                     gtcga$sampletype=="01")]
gtcgatest <- gtcgatest[,!duplicated(gtcgatest$Patient_ID)]

summary(as.data.frame(table(gtcga$Patient_ID)))
summary(as.data.frame(table(gtcgatest$Patient_ID)))

table(gtcgatest$Patient_ID,gtcgatest$sampletype)

# grab dmps where T>NM
x1 <- rowMeans(getBeta(gtcgatest[,gtcgatest$sampletype=="11"]))
x2 <- rowMeans(getBeta(gtcgatest[,gtcgatest$sampletype=="01"]))
hyptcg <- which(x2-x1 >=0.2)

gtcgatest <- gtcgatest[intersect(names(hyptcg),rownames(rlsig.dvmc.df)),] # only perform DMP search among cgs where T>NM
dmptcga <- dmpFinder(getBeta(gtcgatest),
                     type="categorical",
                     pheno=gtcgatest$sampletype)


save(dmptcga,file="dvmcs-tcgadmp-postfilter-q05.rda")
###
