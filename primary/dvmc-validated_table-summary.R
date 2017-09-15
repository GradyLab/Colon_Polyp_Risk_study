#==============================
# validated dvmc summary stats
#==============================

x <- roc.epic[roc.epic$auc>=0.7,]
x <- x[rev(order(x$auc)),]; head(x)

ge <- gfilt.epic[rownames(gfilt.epic) %in% x$cpgid,];
ge <- ge[order(match(rownames(ge),x$cpgid)),]; identical(rownames(ge),x$cpgid)

gt <- gtcgatest[rownames(gtcgatest) %in% x$cpgid,]
gt <- gt[order(match(rownames(gt),x$cpgid)),]; identical(rownames(gt),x$cpgid)

gd <- gfilt.ievora[rownames(gfilt.ievora) %in% x$cpgid,]
gd <- gd[order(match(rownames(gd),x$cpgid)),]; identical(rownames(gd),x$cpgid)

ad <- as.data.frame(getAnnotation(gd))

summarydf <- data.frame(dvmc.id=x$cpgid,
                        auc.epic=x$auc,
                        memean.epic.lr=rowMeans(getBeta(ge[,ge$Risk_Group.index %in% c("0","1")])),
                        memean.epic.hr=rowMeans(getBeta(ge[,ge$Risk_Group.index %in% c("2A","2B","3")])),
                        memean.tcga.t=rowMeans(getBeta(gt[,gt$sampletype=="01"])),
                        memean.tcga.nm=rowMeans(getBeta(gt[,gt$sampletype=="11"])),
                        memean.hm450.lr=rowMeans(getBeta(gd[,gd$Risk_Group.index %in% c("0","1")])),
                        memean.hm450.hr=rowMeans(getBeta(gd[,gd$Risk_Group.index %in% c("2","3")])),
                        cgisland.id=ad$Islands_Name,
                        cgisland.relation=ad$Relation_to_Island,
                        refgene.group=ad$UCSC_RefGene_Group,
                        refgene.name=ad$UCSC_RefGene_Name,
                        chr=ad$chr,
                        pos=ad$pos,
                        strand=ad$strand,
                        stringsAsFactors = FALSE)

write.csv(summarydf,file="dvmc-validated_cgsummaries.csv")
