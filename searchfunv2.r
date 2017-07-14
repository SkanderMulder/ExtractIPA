#f
# s<-"+gfr"
# names_of_db<-names(dataset)
# s<-"+/././././"

# names_of_db[fun_search("+PC.aa-uacr",names_of_db)]

fun_search<-function(s,db)#fun to select names from db
{
cuts<-unlist(strsplit(s,"[|]"))
result_of_search<-NULL


for(zz in 1:length(cuts))
{
# print(zz)
  sav_min<- sav_plus<-vector()
  abc = cuts[zz]

  pluspile<-minpile<-NULL;curtype<-NA
  for(i in 1:nchar(abc))
  {
    CUR<-substring(abc,i,i)
    if(CUR=="-"){curtype <- "-"}
    if(CUR=="+"){curtype <- "+"}
    
    # print(paste(i,CUR))
    
    if(CUR != "-" && CUR != "+")
    {
      if(curtype == "-")minpile<-paste0(minpile,(CUR))
      if(curtype == "+")pluspile<-paste0(pluspile,(CUR))
    }
    
    if(CUR == "-" || CUR == "+")
    {
      if(curtype == "-")minpile<-paste0(minpile,("+-"))
      if(curtype == "+")pluspile<-paste0(pluspile,("+-"))
      # print("hit")
    }
  }
  
  if(length(pluspile)!=0)
  {
  pluspile<-unlist(strsplit(pluspile, "[+-]"));pluspile<-pluspile[(""!=pluspile)]
    for(i in 1:length(pluspile))
  {
    sav_plus<-   c(sav_plus,grep(pluspile[i],(db)))
  }
  }
  
  if(length(minpile)!=0)
  {
  minpile<-unlist(strsplit(minpile, "[+-]"));minpile<-minpile[(""!=minpile)]
    for(i in 1:length(minpile))
  {
    sav_min<-   c(sav_min,grep(minpile[i],(db)))
  }
  }
	#saving results

  sav_plus<-as.numeric(names(table(sav_plus))[which(table(sav_plus) >=  length(pluspile))]) 
  found_in_iter<-sav_plus[which(!sav_plus %in% sav_min)]  #this

  # print(paste("Found:", db[found_in_iter]))
  result_of_search<-c(result_of_search,unlist(found_in_iter)) 
 } 
 
result_of_search<-unique(result_of_search)
#print(paste("Found:", db[sav_plus[which(!sav_plus %in% sav_min)]]))
return(result_of_search)
}


# abc ="-abba+koei-gs-sd+gks"
# db<-c("koeks", "sddgs","gs", "koeiii")
# db[fun("",db)]


