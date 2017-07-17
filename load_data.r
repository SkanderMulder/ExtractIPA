#thisplots objtoman

setwd("C:/Users/smulder/Dropbox/tmpwork/Paper1/Integration/CreateIPAfiles/IPApathway_output/final"); RESURSIVE =TRUE
library('openxlsx');source('C:/Users/smulder/Dropbox/tmpwork/Paper1/Integration/CreateIPAfiles/IPApathway_output/readIPAfiles.r');set.seed(1)
set.seed(1)
library('hellno'); require(googleVis)

#prepare data
names(objtoman)
metselector<-names(objtoman)[c(10,12)]
names(objtoman[[metselector[1]]])

toplot1<-cbind(objtoman[[metselector[1]]][["Canonical Pathways for My Projects"  ]][,c(1,4)],go='Early DKD')
toplot2<-cbind(objtoman[[metselector[2]]][["Canonical Pathways for My Projects"  ]][,c(1,4)],go='Late DKD'); names(toplot2)<-names(toplot2)[3:1]

toplot<-rbind(toplot1,toplot2)
toplot[,1]<-(as.character(toplot[,1]))
toplot[,2]<-as.numeric(as.character(toplot[,2]))
toplot[,3]<-(as.character(toplot[,3]))


###scripts 

pvalue<- -log(0.05); ratio<- 0.1 ;overlapfactor=3
ord<-1:length(allselector)
fixNum<-function(x)as.numeric(as.character(unlist(x)))

