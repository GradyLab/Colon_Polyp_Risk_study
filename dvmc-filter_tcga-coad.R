#======================
# ievora tcga filters
#======================
# author: Sean Maden

# background: This test filters discovered risk markers on methylation in tumor (T) vs. normal matched (NM) TCGA colon samples

#==========================
# Data summaries, subsets
#==========================
table(substr(gtcga$TCGA_Barcode,22,25)) # batches/plates

mdsPlot(getBeta(gtcga),
        sampGroups=substr(gtcga$TCGA_Barcode,22,25),
        sampNames=gtcga$sampletype) 

# exclude tumor sampels with matched normal tissues
gtcgatest <- gtcga[,!(gtcga$Patient_ID %in% gtcga[,gtcga$sampletype=="11"]$Patient_ID &
                        gtcga$sampletype=="01")] 
gtcgatest <- gtcgatest[,!duplicated(gtcgatest$Patient_ID)]

summary(as.data.frame(table(gtcga$Patient_ID))) # up to four patient samples
summary(as.data.frame(table(gtcgatest$Patient_ID))) # only one unique sample per patient

table(gtcgatest$Patient_ID,gtcgatest$sampletype)
table(gtcgatest$sampletype)

dim(gtcgatest)
gtcgatest <- gtcgatest[,gtcgatest$sampletype %in% c("01","11")]; dim(gtcgatest)

#===================================
# Differential Methylation Testing
#===================================

# DMPs where T > NM (T - NM >= 0.2)
x1 <- rowMeans(getBeta(gtcgatest[,gtcgatest$sampletype=="11"]))
x2 <- rowMeans(getBeta(gtcgatest[,gtcgatest$sampletype=="01"]))
hyptcg <- which(x2-x1 >=0.2)

# DMP search among T-hypermethylated CpGs that are also risk markers from discovery
# ievora.ddf a df of DVMCs from testing all HM450k arrays
gtcgatest <- gtcgatest[intersect(names(hyptcg),rownames(ievora.ddf)),] 
dmptcga <- dmpFinder(getBeta(gtcgatest),
                     type="categorical",
                     pheno=gtcgatest$sampletype)

summary(p.adjust(dmptcga$pval,method="BH")) # padj < 2.5e-5
dim(dmptcga) # 671 post-filter markers

#save(dmptcga,file="dvmcs-tcgadmp-postfilter-pbh05.rda")

###
