# ievora 2017 documentation working script

#=====================
# ADD RESECTION DATA
#=====================

# read in resection location info
resec <- read.csv("resections_colocares_senplate2.csv",stringsAsFactors = FALSE)
hset <- read.csv("harmonized_set_colon-rectum.csv",stringsAsFactors=FALSE)

length(intersect(resec[,3],hset$Solution_ID))

# add dna.soln ids
p.gr <- as.data.frame(pData(gr)) # only right colon samples
summary(as.data.frame(table(p.gr$patID))[,2]) # includes replicates from resections, Colocares
p.gl <- as.data.frame(pData(gl)) # only left colon samples

p.gr$DNA.soln <- ""
for(i in 1:nrow(p.gr)){
  p.gr$DNA.soln[i] <- hset[hset$Array_ID==rownames(p.gr)[i],]$Solution_ID 
}

p.gl$DNA.soln <- ""
for(i in 1:nrow(p.gl)){
  p.gl$DNA.soln[i] <- hset[hset$Array_ID==rownames(p.gl)[i],]$Solution_ID 
}

# add resection info
p.gr$resec.loc <- "";p.gl$resec.loc <- ""
for(i in 1:nrow(resec)){
  if(resec$DNA.Database.ID[i] %in% p.gr$DNA.soln){
    p.gr[p.gr$DNA.soln==resec$DNA.Database.ID[i],]$resec.loc <- resec$SpecimenType[i]
  }
  if(resec$DNA.Database.ID[i] %in% p.gl$DNA.soln){
    p.gl[p.gl$DNA.soln==resec$DNA.Database.ID[i],]$resec.loc <- resec$SpecimenType[i]
  }
}

summary(as.data.frame(table(p.gr[p.gr$resec.loc=="",]$patID))[,2]) # no duplicates without resection
summary(as.data.frame(table(p.gl[p.gl$resec.loc=="",]$patID))[,2]) # no duplicates without resection

##########################################
#
# PT1. Grab Location-independent Probes
#
#

#===============================
# SAMPLE SUBGROUPING AND iEVORA
#===============================
library(qvalue)


#=========================
# right colon iEVORA
identical(rownames(p.gr),colnames(gr))
rvar <- rep(0,ncol(gr)); rvar[which(gr$risk_index %in% c(2,3))] <- 1 # make ievora marker
p.gr$rvar <- rvar

# subset, remove duplicate resection samples
# select only unique resections (Location = FZN Distal Adjacent); only 1 SAMPLE PER PATIENT
p.grsub <- p.gr[p.gr$resec.loc %in% c("","FZN Distal Adjacent"),]; 
dim(p.grsub); summary(as.data.frame(table(p.grsub$patID))[,2])
grsub <- gr[,rownames(p.grsub)]; dim(grsub)
#p.grsub <- pData(gsuball.combat[,gsuball.combat$loc2=="right_colon"])
#grsub <- gsuball.combat[,rownames(p.grsub)]

# run iEVORA on R samples
ievora.disc.grsub <- iEVORA(data.m=as.matrix(getBeta(grsub)),
                                    pheno.v=p.grsub$var,
                                    thDV=0.001,
                                    thDM=0.05)

#===========================
# left colon iEVORA
identical(rownames(p.gl),colnames(gl))
lvar <- rep(0,ncol(gl)); lvar[which(gl$risk_index %in% c(2,3))] <- 1 # make ievora marker
p.gl$lvar <- lvar

# subset, remove duplicate resection samples
# select only unique resections (Location = FZN Distal Adjacent); only 1 SAMPLE PER PATIENT
p.glsub <- p.gl[p.gl$resec.loc %in% c("","FZN Distal Adjacent"),]; 
dim(p.glsub); summary(as.data.frame(table(p.glsub$patID))[,2])
glsub <- gl[,rownames(p.glsub)]; dim(glsub)

# run iEVORA on R samples
ievora.disc.glsub <- iEVORA(data.m=as.matrix(getBeta(glsub)),
                                         pheno.v=p.glsub$lvar,
                                         thDV=0.001,
                                         thDM=0.05)

#===============================
# location-independent markers

locind.ievoracpg <- intersect(rownames(ievora.disc.grsub),rownames(ievora.disc.glsub))

###################################
# Whole colon iEVORA

p.allsub <- rbind(p.grsub[,1:13],p.glsub[,1:13])
summary(as.data.frame(table(p.allsub$patID))[,2])

beta.allsub <- cbind(getBeta(grsub[locind.ievoracpg,]),getBeta(glsub[locind.ievoracpg,]))

# discovery subset iEVORA probes
allsub.var <- rep(0,nrow(p.allsub)); allsub.var[which(p.allsub$risk_index %in% c(2,3))] <- 1
table(allsub.var)
p.allsub$var <- allsub.var

disc0 <- sample(rownames(p.allsub[which(allsub.var==0),]),round(2*(89/3),0))
disc1 <- sample(rownames(p.allsub[which(allsub.var==1),]),round(2*(40/3),0))

palld <- p.allsub[rownames(p.allsub) %in% c(disc0,disc1),]; table(palld$var)
pallv <- p.allsub[!rownames(p.allsub) %in% rownames(palld),]; table(pallv$var)

gd <- gsuball.combat[,rownames(palld)]
gv <- gsuball.combat[,rownames(pallv)]


# ievora with discovery samples


##########################
# Batch Effects Testing

beta.allsub <- cbind(getBeta(grsub),getBeta(glsub))
x <- apply(beta.allsub,1,sd)
bmds <- beta.allsub[names(x[x>=quantile(x,seq(0,1,0.01))[100]]),]

mdsPlot(bmds,sampNames=p.allsub$risk_index,sampGroups=p.allsub$batch,legendPos = "topright",legendNCol = 2)

# try combat correction
library(sva)
combat_M <- ComBat(logit2(beta.allsub),batch=p.allsub$batch,mod=model.matrix(~p.allsub$var))

gsuball.combat <- GenomicRatioSet(gr=granges(grsub),
                           Beta=ilogit2(combat_M),
                           M=combat_M,
                           annotation=annotation(grsub),
                           preprocessMethod=preprocessMethod(grsub),
                           pData=p.allsub)


mdsPlot(getBeta(gsuball.combat),
        sampNames=gsuball.combat$risk_index,
        sampGroups=gsuball.combat$batch,
        legendNCol=2,legendPos = "bottomright")

save(gsuball.combat,file="gcombat_ievora_all-loc.rda")




