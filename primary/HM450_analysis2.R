# ievora 091117

# rg.all N=490

rg.rl <- rg.all[,rg.all$Anatomic_Location2!="rectum"]
# left_colon right_colon 
# 199          72

rg.rl.450 <- rg.rl[,rg.rl$Array_platform=="HM450K"]
table(rg.rl.450$Anatomic_Location2)
# left_colon right_colon 
# 145          71 


# resections
patall.ccar <- colnames(rg.rl.450[,rg.rl.450$Project=="CCAR" & 
                           rg.rl.450$Run_Date %in% c("20120412","20150103","20161108")])

resec <- read.csv("resections_colocares_senplate2.csv",stringsAsFactors = FALSE)
resecpat <- c();patid <- unique(resec$ParticipantID)
for(i in 1:length(patid)){
  rx <- resec[resec$ParticipantID==patid[i],]
  if("FZN Distal Adjacent" %in% rx$SpecimenType){
    resecpat[i] <- rx[rx$SpecimenType=="FZN Distal Adjacent",]$DNA.Database.ID
  } else{
    resecpat[i] <- rx[!duplicated(rx$ParticipantID),]$DNA.Database.ID
  }
}

array.ccar <- c(patall.ccar,colnames(rg.rl.450[,rg.rl.450$Solution_ID %in% resecpat]))
array.all <- c(colnames(rg.rl.450[,rg.rl.450$Project %in% c("GICA","AKNC")]), array.ccar)
rgsub.rl.450 <- rg.rl.450[,array.all]

rgsub.rl.450 <- rgsub.rl.450[,rgsub.rl.450$Institution_tissue %in% c("UPMC","UWAS")]
table(rgsub.rl.450$Anatomic_Location2,
      rgsub.rl.450$Risk_Group.index)

save(rgsub.rl.450,file="rgsub_rlcolon_ievora.rda")

summary(as.data.frame(table(rg.rl.450$PatientID_cohort))[,2])


# preprocess

library(methyPre)

ievora.var <- rgsub.rl.450$Risk_Group.index
ievora.var <- ifelse(ievora.var %in% c("2","3"),"1","0")
data.frame(rgsub.rl.450$Risk_Group.index,ievora.var)

mpre <- methyPre(rgsub.rl.450,
                 workflow=c("norm","pfilt","map","minfifilt","crxcgfilt"),
                 normfun=c("SWAN","illumina"),
                 pfilt=TRUE,detPcutoff = 0.05,
                 minfiFilt=TRUE,uFilt=FALSE,crxcgFilt = "hm450")

batchpre <- methyPre(mpre[[6]],
                 workflow=c("batch"),
                 batchCorrect=TRUE,
                 bcCovariate=ievora.var,
                 bcBatch=rgsub.rl.450$Run_Date)


save(mpre,file="esets_preprocessing_ievora-RL-test.rda")
