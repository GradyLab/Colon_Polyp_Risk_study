#==============================
# main iEVORA marker discovery
#==============================
# Author: Sean Maden

# use all HM450k arrays at all colon locations
g <- gfilt.ievora
ievoravar <- g$Risk_Group.index
ievoravar <- ifelse(ievoravar %in% c("2","3"),1,0)
table(ievoravar); dim(g)

gdisc <- g[rownames(g) %in% rldvc.list,]

ievora.disc <- iEVORA(getBeta(gdisc),ievoravar,thDV=1,thDM=1)
ievora.ddf <- as.data.frame(ievora.disc)

###
