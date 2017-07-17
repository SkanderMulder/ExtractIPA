library('xlsx');library('hellno')
list.files()

# x<-read.xlsx('allMet_alb_est.xls',1,startRow=2)
source("searchfunv2.r")

#forloop to get the listnames
convertRawIPAoutput<-function(x,foldername='+.>Skander_shared')
{
indexloc<-fun_search(foldername,as.character(x[,1]))
# print('names of listobj:');
obj=vector('list',length(indexloc))
listnames<-NULL;for(i in 1:length(indexloc))
{
pos = regexpr('->', x[indexloc,1][i])[1]
# print(substr(x[indexloc,1][i],1,pos-1))
listnames<-c(listnames,substr(x[indexloc,1][i],1,pos-1))
}
#created list object
names(obj)<-listnames
#save settings for loop #2
headrloc<-indexloc+1
startloc<-indexloc+2
endloc<-c(indexloc[(2:length(indexloc))]-1, c(dim(x)[1]))
#endset
for(y in 1:length(listnames))
{
if((endloc[y]-indexloc[y]) >2)
{
# ;print(y)
newnames<-as.character(unlist(x[headrloc[y],][fun_search('+.-NA',x[headrloc[y],])]))
tmp<-x[startloc[y]:endloc[y],fun_search('+.-NA',x[headrloc[y],])]
colnames(tmp)<-(newnames);head(tmp)
obj[[y]]<-tmp
}
}
return(obj)
}


detach('package:openxlsx',)
iteror<-list.files(recursive = RESURSIVE);iteror<-iteror[fun_search('+xls-xlsx',iteror)];vts<-vector('list',length(iteror))
for(zz in 1:length(iteror))
{
vts[[zz]]<-convertRawIPAoutput(read.xlsx(iteror[zz],1,startRow=2))
print(paste0('ZZ is ',zz,' for :',iteror[zz]));flush.console()
}
names(vts)<-iteror;objtoman<-vts

print(objtoman)

names(objtoman)
objtoman$allMet_alb_est.xls['Analysis Ready Molecules for My Projects']
names(objtoman$allMet_alb_est.xls)

names(objtoman)
names(objtoman[[1]])
objtoman[[1]][[1]]


tosel<-2;names(objtoman[[1]])[[tosel]];head(objtoman[[1]][[tosel]])

