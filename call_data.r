source("searchfunv2.r")
source('load_data.r')
setwd("C:/Users/smulder/Dropbox/tmpwork/Paper1/Integration/CreateIPAfiles/IPApathway_output/ExtractIPA");source('functionSankey.r')

savedObjtoman<-(objtoman)

x=1:3

#all
allselector<-c( "allMet_renal30conf0HR.xls" ,'allProt_renal30conf0HR.xls','dkd_glom_ellu_egfr.xls','dkd_tub_ellu_egfr.xls')
allselector %in% names(objtoman)
allselectorSYNC<-c('Serum Metabolites','Urinary Proteins','Glomerulus Transcriptome','Tubule Transcriptome')
CanonicalPathway_all<-CanonicalPathwayPlot(turnround2=TRUE,coloverload= rep(rgb(200/255,200/255,200/255),10000),turnaround=FALSE,overlapfactor=2,p=0.05,transR=FALSE,OLAP=0)#,coloverload=c('gray','red','blue','black')) 
#early
allselector<-c( "earlyMet_renal30conf0HR.xls" ,'earlyProt_renal30conf0HR.xls','dkd_glom_ellu_egfr.xls','dkd_tub_ellu_egfr.xls')
allselector %in% names(objtoman)
allselectorSYNC<-c('Serum Metabolites','Urinary Proteins','Glomerulus Transcriptome','Tubule Transcriptome')
CanonicalPathway_early<-CanonicalPathwayPlot(turnround2=TRUE,coloverload= rep(rgb(200/255,200/255,200/255),10000),turnaround=FALSE,overlapfactor=1,p=0.05,transR=FALSE,OLAP=0)

#late
allselector<-c( "lateMet_renal30conf0HR.xls" ,'lateProt_renal30conf0HR.xls','dkd_glom_ellu_egfr.xls','dkd_tub_ellu_egfr.xls')
allselector %in% names(objtoman)
allselectorSYNC<-c('Serum Metabolites','Urinary Proteins','Glomerulus Transcriptome','Tubule Transcriptome')
CanonicalPathway_late<-CanonicalPathwayPlot(turnround2=TRUE,coloverload= rep(rgb(200/255,200/255,200/255),10000),turnaround=FALSE,overlapfactor=1,p=0.05,transR=FALSE,OLAP=0)

#baselinestate
allselector<-c( "allMet_egfr_est.xls" ,'allProt_egfr_est.xls','dkd_glom_ellu_egfr.xls','dkd_tub_ellu_egfr.xls')
allselector %in% names(objtoman)
allselectorSYNC<-c('Serum Metabolites','Urinary Proteins','Glomerulus Transcriptome','Tubule Transcriptome')
CanonicalPathway_bl<-CanonicalPathwayPlot(turnround2=TRUE,coloverload= rep(rgb(200/255,200/255,200/255),10000),turnaround=FALSE,overlapfactor=1,p=0.05,transR=FALSE,OLAP=0)

#Combendpoints
allselector<-c( '000grps/earlyMetEGRGRP.xls',"000grps/earlyprotEGFRGRP.xls","000grps/lateprotEGFRGRP.xls","000grps/lateMetEGFRGRP.xls" )#,'dkd_glom_ellu_egfr.xls','dkd_tub_ellu_egfr.xls')
allselector %in% names(objtoman)
# allselectorSYNC<-c('Direct Metabolites','Direct Proteins','Sun-Macro Metabolites','Sun-Macro Proteins')#,'Glomerulus Transcriptome','Tubule Transcriptome')
allselectorSYNC<-c('Direct','Direct','Sun-Macro','Sun-Macro')#,'Glomerulus Transcriptome','Tubule Transcriptome')
CanonicalPathway_bl<-CanonicalPathwayPlot(CC=ccc,overlapfactor=,p=0.05 )

#COMBINEDENDPTS
allselector<-c( "lateMet_renal30conf0HR.xls","lateProt_renal30conf0HR.xls","000grps/lateMetEGFRGRP.xls","000grps/lateprotEGFRGRP.xls" )#,'dkd_glom_ellu_egfr.xls','dkd_tub_ellu_egfr.xls')
allselector %in% names(objtoman)
allselectorSYNC<-c('Sun-Macro1','Sun-Macro1','Sun-Macro2','Sun-Macro2')#,'Glomerulus Transcriptome','Tubule Transcriptome')

CP1<-CanonicalPathwayPlot(CC ='1',plot=TRUE,overlapfactor=1,p=0.05,turnaround=FALSE,coloverload= c('red','green') )

allselector<-c( "earlyMet_renal30conf0HR.xls","earlyProt_renal30conf0HR.xls","000grps/earlyMetEGRGRP.xls","000grps/earlyprotEGFRGRP.xls" )#,'dkd_glom_ellu_egfr.xls','dkd_tub_ellu_egfr.xls')
allselector %in% names(objtoman)
allselectorSYNC<-c('Direct1','Direct1','Direct2','Direct2')#,'Glomerulus Transcriptome','Tubule Transcriptome')

CP2<-CanonicalPathwayPlot(CC ='2',plot=FALSE,overlapfactor=1,p=0.05,turnaround=FALSE,coloverload= c('red','green') )
# names(CP2); names(CP2)<-c('Go','Pathway','Ratio')

TOPLOT<-rbind(CP2,CP1)
sum(table(TOPLOT[,2])==4) 





allselectorSYNC<-1:10

library('openxlsx')

wb2<-createWorkbook()
addWorksheet(wb2, "allmet", gridLines = TRUE)
addWorksheet(wb2, "earlymet", gridLines = TRUE)
addWorksheet(wb2, "latemet", gridLines = TRUE)
addWorksheet(wb2, "allprot", gridLines = TRUE)
addWorksheet(wb2, "earlyprot", gridLines = TRUE)
addWorksheet(wb2, "lateprot", gridLines = TRUE)

allselector<-c( "allMet_renal30conf0HR.xls" ,'000grps/allMetEGFRGRP.xls');allselector %in% names(objtoman)
writeData(wb2,"allmet",unique(CanonicalPathwayPlot(CC='',overlapfactor=1,p=0.05 ,plot=FALSE)$Go))

allselector<-c( "earlyMet_renal30conf0HR.xls" ,'000grps/earlyMetEGRGRP.xls');allselector %in% names(objtoman)
writeData(wb2,"earlymet",unique(CanonicalPathwayPlot(CC='',overlapfactor=1,p=0.05 ,plot=FALSE)$Go))

allselector<-c( "lateMet_renal30conf0HR.xls" ,'000grps/lateMetEGFRGRP.xls');allselector %in% names(objtoman)
writeData(wb2,"latemet",unique(CanonicalPathwayPlot(CC='',overlapfactor=1,p=0.05 ,plot=FALSE)$Go))

allselector<-c( "allProt_renal30conf0HR.xls" ,'000grps/allprotEGFRGRP.xls');allselector %in% names(objtoman)
writeData(wb2,"allprot",unique(CanonicalPathwayPlot(CC='',overlapfactor=1,p=0.05 ,plot=FALSE)$Go))

allselector<-c( "earlyProt_renal30conf0HR.xls" ,'000grps/earlyprotEGFRGRP.xls');allselector %in% names(objtoman)
writeData(wb2,"earlyprot",unique(CanonicalPathwayPlot(CC='',overlapfactor=1,p=0.05 ,plot=FALSE)$Go))

allselector<-c( "lateProt_renal30conf0HR.xls" ,'000grps/lateprotEGFRGRP.xls');allselector %in% names(objtoman)
writeData(wb2,"lateprot",unique(CanonicalPathwayPlot(CC='',overlapfactor=1,p=0.05 ,plot=FALSE)$Go))
# saveWorkbook(wb2, file = paste(trialname,name), overwrite = TRUE)

saveWorkbook(wb2, 'ZZZoput.xlsx',overwrite=TRUE)



###NEW

#early
allselector<-c('000grps/earlyMetEGRGRP.xls',"000grps/earlyprotEGFRGRP.xls","earlyMet_renal30conf0HR.xls" ,'earlyProt_renal30conf0HR.xls','000grps/lateMetEGFRGRP.xls',"000grps/lateprotEGFRGRP.xls", "lateMet_renal30conf0HR.xls" ,'lateProt_renal30conf0HR.xls','dkd_glom_ellu_egfr.xls','dkd_tub_ellu_egfr.xls')
allselector %in% names(objtoman)
allselectorSYNC<-c(rep('Early DKD',4),rep('Late DKD',4),rep('Glom Transcriptome',1),rep('Tub Transcriptome',1))


colstovl<-rep(rgb(200/255,200/255,200/255),10000)

COLOBJ<-matrix(c('Early DKD','green','Late DKD','red'),ncol=2,byrow=TRUE);


colstovl<-rep(rgb(200/255,200/255,200/255),10000)
colstovl[999]<-'black'; colstovl[10]<-'red';colstovl[1]<-'blue'

CanonicalPathway_eaaa<-CanonicalPathwayPlot(turnaround=FALSE,overlapfactor=1,p=0.05,coloverload =colstovl,OLAP=3,COLOBJ=COLOBJ)#,coloverload=rev(c('green','orange','gray')) )




wb2<-createWorkbook()
for(i in 0:3)
{
print(i)
addWorksheet(wb2, paste0("olap",i), gridLines = TRUE)
CanonicalPathway_eaaa<-CanonicalPathwayPlot(plots=TRUE,turnaround=FALSE,overlapfactor=1,p=0.05,coloverload =colstovl,OLAP=1)#,coloverload=rev(c('green','orange','gray')) )
print(paste(i, dim(CanonicalPathway_eaaa)))
wholelist=unique(c(CanonicalPathway_eaaa$Go, CanonicalPathway_eaaa$Pathway))
writeData(wb2,paste0("olap",i),data.frame(wholelist))
}

setwd("C:/Users/smulder/Dropbox/tmpwork/Paper1/Integration/CreateIPAfiles/IPApathway_output/ELLU"); RESURSIVE =FALSE
saveWorkbook(wb2, 'wholelist.xlsx',overwrite=TRUE)
