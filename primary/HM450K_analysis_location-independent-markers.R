#===================================================
# location-indpendent and dependent ievora markers
#===================================================
# Author: Sean Maden

ievora.var <- gfilt_ievora$Risk_Group.index
ievora.var <- ifelse(ievora.var %in% c("3","2"),1,0)

rsample <- gfilt_ievora[,gfilt_ievora$Anatomic_Location2=="right_colon"]
rvar <- ievora.var[which(gfilt_ievora$Anatomic_Location2=="right_colon")]

lsample <- gfilt_ievora[,gfilt_ievora$Anatomic_Location2=="left_colon"]
lvar <- ievora.var[which(gfilt_ievora$Anatomic_Location2=="left_colon")]


# ievora probes

rbt <- apply(getBeta(rsample),1,function(x){return(bartlett.test(x=x,g=rvar))})
rbtp <- c(); 
for(i in 1:length(rbt)){
  rbtp[i] <- rbt[[i]]$p.value;
  message(paste0(round(100*i/length(rbt),3),"%"))
  }
rbtpbh <- p.adjust(rbtp,method="BH"); length(rbtpbh[rbtpbh<=0.05])
rdvclist <- names(rbt[which(rbtpbh<=0.05)])

lbt <- apply(getBeta(lsample),1,function(x){return(bartlett.test(x=x,g=lvar))})
lbtp <- c(); 
for(i in 1:length(lbt)){
  lbtp[i] <- lbt[[i]]$p.value;
  message(paste0(round(100*i/length(lbt),3),"%"))
}
lbtpbh <- p.adjust(lbtp,method="BH"); length(lbtpbh[lbtpbh<=0.05])
ldvclist <- names(lbt[which(lbtpbh<=0.05)])

rldvc.list <- intersect(ldvclist,rdvclist)
save(rldvc.list,file="DVC-list_RLColon-Overlap_BT-PBH-05.rda")

markersoi <- c("cg00727675","cg06422261","cg24148882","cg11291003","cg11111460","cg19008720","cg19980151","cg03756778","cg12052265")
length(intersect(markersoi,rldvc.list))

###
