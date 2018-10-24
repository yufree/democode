library('VennDiagram')
library('scales')

list <- list(v1,v2,v3,v4,v5,v6,v7)
venn.diagram(list,'venn.png',fill=hue_pal()(5),margin=c(.05,.05,.05,.05),imagetype = 'png')