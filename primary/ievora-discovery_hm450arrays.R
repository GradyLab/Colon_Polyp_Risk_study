#=================================
# iEVORA marker discovery tests
#=================================
# Author: Sean Maden

g <- gfilt.ievora
ievoravar <- g$Risk_Group.index
ievoravar <- ifelse(ievoravar %in% c("2","3"),1,0)
table(ievoravar); dim(g)

gdisc <- g[rownames(g) %in% rldvc.list,]

ievora.disc <- iEVORA(getBeta(gdisc),ievoravar,thDV=1,thDM=1)
