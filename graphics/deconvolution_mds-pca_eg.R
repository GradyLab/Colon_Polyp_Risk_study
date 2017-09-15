#=============================================================
# iEVORA deconvolution analysis, at most variable CpGs/DVCs
#=============================================================
# Author: Sean Maden

# Using all CpGs: MVPs
x <- apply(getBeta(gfilt.ievora),1,var)
x2 <- x[!is.na(x)]
x2 <- x2[rev(order(x2))]
summary(x2)
summary(x2[1:1000])
gpca <- gfilt.ievora[names(x2[1:1000]),]

# mds plot
mdsPlot(getBeta(gpca),sampGroups=gpca$Anatomic_Location2)
# principal component analysis plot
autoplot(prcomp(t(getBeta(gpca))),#data=as.data.frame(pData(gpca)),
         colour=ifelse(colx=="right_colon","orange","green"))

# Using location-independent DVCs
x <- apply(getBeta(gfilt.ievora[rldvc.list,]),1,var)
x2 <- x[!is.na(x)]
x2 <- x2[rev(order(x2))]
summary(x2)
summary(x2[1:1000])
gpca <- gfilt.ievora[names(x2[1:1000]),]

# mds plot
mdsPlot(getBeta(gpca),sampGroups=gpca$Anatomic_Location2,main="MDS at Most Variable DVCs")
# principal component analysis plot
autoplot(prcomp(t(getBeta(gpca))),#data=as.data.frame(pData(gpca)),
         colour=ifelse(colx=="right_colon","orange","green"),main="PCA at Most Variable DVCs")

###
