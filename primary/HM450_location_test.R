# location-indpendent and dependent ievora markers

gfilt_ievora <- gfilt_ievora[,!duplicated(gfilt_ievora$PatientID_cohort)]; dim(gfilt_ievora)
table(gfilt_ievora$Risk_Group.index,gfilt_ievora$Anatomic_Location2)
save(gfilt_ievora,file="gfilt-ievora-RL-test-uniquepat.rda")

ievora.var <- gfilt_ievora$Risk_Group.index
ievora.var <- ifelse(ievora.var %in% c("3","2"),1,0)

rsample <- gfilt_ievora[,gfilt_ievora$Anatomic_Location2=="right_colon"]
rvar <- ievora.var[which(gfilt_ievora$Anatomic_Location2=="right_colon")]

lsample <- gfilt_ievora[,gfilt_ievora$Anatomic_Location2=="left_colon"]
lvar <- ievora.var[which(gfilt_ievora$Anatomic_Location2=="left_colon")]


# ievora probes

rdvmc <- iEVORA(getBeta(rsample),rvar,thDV=0.05,thDM=1)
ldvmc <- iEVORA(getBeta(lsample),lvar,thDV=0.05,thDM=1)
rldvmc <- intersect(rownames(ldvmc),rownames(rdvmc))

old.markers <- c("cg00727675",
                 "cg06422261",
                 "cg24148882",
                 "cg11291003",
                 "cg11111460",
                 "cg19008720",
                 "cg19980151",
                 "cg03756778",
                 "cg12052265")

length(intersect(old.markers,rldvmc))


# alternate testing

rdv <- t(apply(getBeta(rsample),1,doDV,rvar));

# R colon tests
rDVC <- t(apply(getBeta(rsample),1,doDV,rvar));
rdvcdf <- as.data.frame(rDVC)
rdvcdf$padj.bh <- p.adjust(rdvcdf[,2],method="BH")

# L colon tests
lDVC <- t(apply(getBeta(lsample),1,doDV,lvar));
ldvcdf <- as.data.frame(lDVC)
ldvcdf$padj.bh <- p.adjust(ldvcdf[,2],method="BH")

# get intersect of significantly variable probes

rldvc <- intersect(rownames(rdvcdf[rdvcdf$padj.bh<0.05,]),rownames(ldvcdf[ldvcdf$padj.bh<0.05,])) # 60k markers

save(rldvc,file="rlcolon-common-DVClist-padjBH.rda")
#=====================





