#thisplots objtoman


setwd("C:/Users/smulder/Dropbox/tmpwork/Paper1/Integration/CreateIPAfiles/IPApathway_output/final"); RESURSIVE =TRUE
library('openxlsx');source('C:/Users/smulder/Dropbox/tmpwork/Paper1/Integration/CreateIPAfiles/IPApathway_output/readIPAfiles.r');set.seed(1)
;set.seed(1)
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

CanonicalPathwayPlot<-function(turnround2=TRUE,dimension = c(1000,1000),CC='cc',turnaround=TRUE,overlapfactor=1,p=0.05, plots=TRUE, coloverload =NA,transR=TRUE,OLAP=1,COLOBJ=NA)
{
#STARTSELHERE
pvalue<- -log(p)

# allselectorSYNC=allselector<-NA
if(turnaround)
{
allselector<-rev(allselector)
allselectorSYNC<-rev(allselectorSYNC)
}

#loop to create object
dfr<-data.frame(Pathway=NULL,Logratio =NULL, Zscore=NULL,Ratio =NULL, Molecules=NULL,Go=NULL)
ites_sav_nrow<-NULL
for(i in 1:length(allselector))
{


v1<-objtoman[[allselector[i]]][[1]]
v1<-subset(v1, fixNum(v1[" -log(p-value)"])  > pvalue ) #& fixNum(v1['Ratio']) <ratio)
names(v1)<-c('Pathway','Logratio','Zscore','Ratio','Molecules')

v1$Go <- allselectorSYNC[i]
#if(i %% 2 == 0)#isevennumber
nforranscrp<-grep('ranscri',allselectorSYNC)


if(transR)
if(!i %in% nforranscrp)
{
	v1$Go <- v1$Pathway
	v1$Pathway <- allselectorSYNC[i]
}
ites_sav_nrow<-c(ites_sav_nrow,nrow(v1))
dfr<-rbind(dfr,v1)
}

#preparingdata

TOPLOT<-(dfr[,c(1,6,4)])
TOPLOT[,2]<-as.character(as.factor(TOPLOT[,2]))
TOPLOT[,3]<-as.numeric (as.character(as.factor(TOPLOT[,3])))
TOPLOT[,1]<-as.character(as.factor(TOPLOT[,1]))

funs<-function(x) {return(sum(c(TOPLOT[,2],TOPLOT[,1]) ==x)) }
##ALLTHING
name<-data.frame(a=c(TOPLOT[,2],TOPLOT[,1]))
tomerg<-data.frame(table(name))
savelookup<-merge(x=name,by.x='a',y=tomerg,by.y='name')
return_number<-function(x)
{
if(length(x)>1)return(sapply(x,return_number))
rval<-(savelookup[which(savelookup[,1] == x),2])[1]
;return(rval)
}

###
overlapfactors<-pmin(return_number(TOPLOT[,2]),return_number(TOPLOT[,1]))  >=overlapfactor  
overlapfactors13 <- 1:length(overlapfactors) %in% fun_search('+Transcripto',names(overlapfactors))
TOPLOTs<-subset(TOPLOT, overlapfactors13 | overlapfactors);if(!plots)return(TOPLOT);TOPLOTs
#endpreparedta

if(turnround2)names(TOPLOT)<-names(TOPLOT[,c(2,1,3)])
#musthaveoverlapWtranscriptome

olap1<-TOPLOTs[,1] %in%  allselectorSYNC  |  TOPLOTs[,1] %in% TOPLOTs[,2] #checkifserumin trans
olap2<-TOPLOTs[,2] %in%  allselectorSYNC  | TOPLOTs[,2] %in%  TOPLOTs[,1]#checks if trans in serum


if(!  OLAP==0)
{
if(OLAP==1)TOPLOTs<-TOPLOTs[which(olap1),]
if(OLAP==2)TOPLOTs<-TOPLOTs[which(olap2),]
if(OLAP==3) TOPLOTs<-TOPLOTs[which(olap1 &olap2),]
}
whatTranscript<-grep('ranscr',TOPLOTs[,2])

whhatNONTranscript<-which(!1:nrow(TOPLOTs) %in% whatTranscript)
devfact<-sum(TOPLOTs$Ratio [whatTranscript]) / sum(TOPLOTs$Ratio [whhatNONTranscript])
TOPLOTs$Ratio [whhatNONTranscript]<-TOPLOTs$Ratio [whhatNONTranscript]*devfact

TOPLOT$Ratio[whhatNONTranscript]<-TOPLOT$Ratio[whhatNONTranscript]

conv_listobj_invec<-function(listobj)
{
appender<-NULL
for(i in 1:length(listobj))
{
appender<-c(appender,rep(i,length(listobj[[i]])))
}
return(appender)}

#assignanobjtodoesthelinkcolors
listobj<-lapply(paste0('+',allselectorSYNC), fun_search,db=apply(data.frame(a=names(return_number(TOPLOTs[,1])) , b=names(return_number(TOPLOTs[,2]))),1,toString))
coltosel<-conv_listobj_invec(listobj)


# SLEFT
to=TOPLOTs[,2];fr=TOPLOTs[,1]





Sleft<-function(Y=1,ztk=NULL,SEARCHID='NA')
{
x=NA

if(!is.character(Y))x=fr[Y]
if(is.character(Y)) x=Y
print(paste('search=', x))

loc<-unlist(lapply(x, function(x)which( x==fr)))
# return(loc)
ztk<-c(ztk,loc)
tos<-unlist(lapply(to[loc],function(x)which(fr %in% x)))
if(length(tos)==0)return(ztk)
if( SEARCHID %in% to[loc])return('yes')
Sleft(tos,ztk)
}


Sright<-function(Y=1,ztk=NULL,SEARCHID='NA')
{
x=NA

if(!is.character(Y))x=to[Y]
if(is.character(Y)) x=Y
print(paste('search=', x))

loc<-unlist(lapply(x, function(x)which( x==to)))
# return(loc)
ztk<-c(ztk,loc)
tos<-unlist(lapply(to[loc],function(x)which(to %in% x)))
if(length(tos)==0)return(ztk)
if( SEARCHID %in% to[loc])return('yes')
Sleft(tos,ztk)
}



funwrapsleft<-function(x=1)
{
tosearch<-Sleft(x)
search_result<-funwrapsleft(Sleft(to[tosearch]))
if(length(search_result ==0)) x=x+1
}
# funwrapsleft()



inputtocolnodes<-unique(as.vector(apply(TOPLOTs[,1:2],1,byrowdf<-function(OBJ){QQ<-as.vector(unlist(OBJ));colnames(QQ)<-NULL; if(length(QQ) ==0)return(0) ;return(QQ)})))
obj_pw_col<-rbind(obj_pw_col,cbind( unique(TOPLOTs[,1]),rep('FFFFFF',length(unique(TOPLOTs[,1])))))

colors_node <-sapply(inputtocolnodes ,retcolmerger)

##HEREGETORDER
names_pahtwayorder<-unlist(data.frame(t(TOPLOT[,1:2])))
names_pahtwayorder<-names_pahtwayorder[!duplicated(names_pahtwayorder)]
names(names_pahtwayorder)<-NULL; names_pahtwayorder

if(!is.na(COLOBJ))for(i in 1:nrow(COLOBJ))
{
colors_linkz[which(COLOBJ[,i] == names_pahtwayorder )] <- COLOBJ[,i]
}
#END

gry<-rgb(200/255,200/255,200/255);bluelight<-'#4e9afc'; purplelight<-'#ff3db4' ;greenlight<-'#90f77b'
colors_linkz<-coloverload

colors_linkz[which(is.na(colors_linkz) | is.na(colors_linkz) =='NA')]<- '000000'
colors_node[which(is.na(colors_node) | is.na(colors_node) =='NA')]<- '000000'

colors_link_array <- paste0("[", paste0("'", colors_linkz,"'", collapse = ','), "]")
colors_node_array <- paste0("[", paste0("'", colors_node,"'", collapse = ','), "]")
# print(colors_link_array)


##HEREGETORDER
names_pahtwayorder<-unlist(data.frame(t(TOPLOT[,1:2])))
names_pahtwayorder<-names_pahtwayorder[!duplicated(names_pahtwayorder)]
names(names_pahtwayorder)<-NULL; names_pahtwayorder




opts <- paste0("{
        link: { colorMode: 'gradient',
                colors: ", colors_link_array ," },
        node: { colors: ", colors_node_array ,",width:10,labelPadding: -0, nodePadding: 3,   label: {
      fontName: 'Times-Roman',
      fontSize: 25,
      color: '#000',
      bold: true,
      italic: false
    }}}")
print(opts)
	
	if(plots)
plot(
  gvisSankey(TOPLOTs[,1:3] ,from="Go", 
             to="Pathway",# weight='Ratio',
             options=list(
               height=dimension[1],width=dimension[2],
               sankey=opts,
			   title="Women",titleTextStyle="{color:'red', fontName:'S', fontSize:99}"
               ))
)



#TOPLOTs$Pathway<-paste0(CC,TOPLOTs$Pathway)
TOPLOTs$Pathway<-paste0(TOPLOTs$Pathway)
return(TOPLOTs)

}

getallpathway<-function()
{
appvec<-NULL;for(i in 1:14)
{
# print(i)
ap<-unlist(as.character(objtoman[[i]][[1]]$"Ingenuity Canonical Pathways"))

appvec<-c(appvec,ap)
}
return(appvec)
}

fun<-function(x)
{
ccc<-substr(rgb(x[1]/7, x[2]/7, x[3]/7, 1, names = NULL, maxColorValue = 1),2,7)
return(ccc)
}

v1<-unique(as.vector(getallpathway()))
v2<-sample(apply((expand.grid(0:7,0:7,0:7)),1,fun),length(v1)) ###?
obj_pw_col<-cbind(v1,v2)
obj_pw_col

pathway="Protein Kinase A Signaling" 
retcolmerger<-function(pathway)
{
robj<-obj_pw_col[,2][which(obj_pw_col[,1]==pathway)][1]
if(length(robj)==0)return('black')
return(robj)
}





# function(x)
#endscripts


savedObjtoman<-(objtoman)


x=1:3

#all
allselector<-c( "allMet_renal30conf0HR.xls" ,'allProt_renal30conf0HR.xls','dkd_glom_ellu_egfr.xls','dkd_tub_ellu_egfr.xls')
allselector %in% names(objtoman)
allselectorSYNC<-c('Serum Metabolites','Urinary Proteins','Glomerulus Transcriptome','Tubule Transcriptome')
CanonicalPathway_all<-CanonicalPathwayPlot(turnround2=TRUE,coloverload= rep(rgb(200/255,200/255,200/255),10000),turnaround=FALSE,overlapfactor=1,p=0.05,transR=FALSE,OLAP=0)#,coloverload=c('gray','red','blue','black')) 

#early
allselector<-c( "earlyMet_renal30conf0HR.xls" ,'earlyProt_renal30conf0HR.xls','dkd_glom_ellu_egfr.xls','dkd_tub_ellu_egfr.xls')
allselector %in% names(objtoman)
allselectorSYNC<-c('Serum Metabolites','Urinary Proteins','Glomerulus Transcriptome','Tubule Transcriptome')
CanonicalPathway_early<-CanonicalPathwayPlot(CC='early DKD',overlapfactor=2,p=0.05 )

#late
allselector<-c( "lateMet_renal30conf0HR.xls" ,'lateProt_renal30conf0HR.xls','dkd_glom_ellu_egfr.xls','dkd_tub_ellu_egfr.xls')
allselector %in% names(objtoman)
allselectorSYNC<-c('Serum Metabolites','Urinary Proteins','Glomerulus Transcriptome','Tubule Transcriptome')
CanonicalPathway_late<-CanonicalPathwayPlot(CC='late DKD',overlapfactor=2,p=0.05 )

#baselinestate
allselector<-c( "allMet_egfr_est.xls" ,'allProt_egfr_est.xls','dkd_glom_ellu_egfr.xls','dkd_tub_ellu_egfr.xls')
allselector %in% names(objtoman)
allselectorSYNC<-c('Serum Metabolites','Urinary Proteins','Glomerulus Transcriptome','Tubule Transcriptome')
CanonicalPathway_bl<-CanonicalPathwayPlot(CC=ccc,overlapfactor=,p=0.05 )

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
