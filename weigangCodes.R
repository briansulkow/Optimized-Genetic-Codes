#GOAL- Merge Dr Qiu's codon list with bins and  paths of amino acids

###SECTION I 
##Retrieve Files= 

##Weiga 
weigangCode<-read.csv("C:/Users/brian/OneDrive/Documents/AANeighbors/weigangCodonList.csv", sep=" ",header=FALSE)
weigangCode<-data.frame(weigangCode[,13],weigangCode[,1])
ncol(weigangCode)
x<-order(weigangCode[,1])
weigangCode[order(weigangCode[,1]),]
weigangCode

#get the titvx3 file names
fileList<-list()
for (i in 1:length(aamapsreOrder)){
 name<-paste("hydro",i,sep="_")
 name<-paste(name,'txt', sep=".")
 path<-paste("C:/Users/brian/OneDrive/Documents/AANeighbors/parameterTesting/oneParameter/oneParameterHydro",name,sep="/")
 fileList[[i]]<-path
 }

subset(weigangCode,weigangCode[,1]==8)


#read them into tables

tiFiles<-lapply(fileList,read.table,header=TRUE)

#
#Function to Take the aa path out data set
unique(tiFiles[[2]][,2])
uniqueAAPath<-function(data){
  a<-unique(data[,2]) 
  return(a)
}
weigangCode

wc<-subset(weigangCode,weigangCode[,1]!=8)
nrow(wc)
#apply function to all files
tiCol<-lapply(tiFiles,uniqueAAPath)
tiCol
#function to alligncode
tiCol[[1]]
weigangCode[1,]
a<-insertXInAamapsreorder(aamapsreOrder[[1]])
a[weigangCode[3,1],]
#AAweigangCodedf

AAWeigangCodedf<-function(data) {
  AminoAcids<-c()
  for (i in 1:64) {
    AminoAcids[i]<-data[weigangCode[i,1],3]
  } 
  df<-data.frame(weigangCode[,2],AminoAcids)
  return(df)
}
fullaaList[[1]]

completeCode<-lapply(fullaaList,AAWeigangCodedf)
completeCode[[1]]
b<-AAWeigangCodedf(a)
subset(b, b[,2]=="X")
#lapply  allign AAWeigangCodedf

weigangHydroVol<-lapply(aamapsreOrder,AAWeigangCodedf)



#Version for 20x20 codes
wc<-subset(weigangCode,weigangCode[,1]!=8)  #take out stop codons and aa's




AAWeigangCodedf2020<-function(data) {
  AminoAcids<-c()
  for (i in 1:61) {
    AminoAcids[i]<-data[wc[i,1],3]
  } 
  df<-data.frame(wc[,2],AminoAcids)
  return(df)
}


aamapsreOrder[[1]][,3]

#lapply  allign AAWeigangCodedf2020

weigangHydroVol<-lapply(aamapsreOrder,AAWeigangCodedf2020)

weigangHydroVol[[1]]
length(weigangHydroVol)


# check if all paths contain 21 amino acids
for ( i in 1:length(weigangHydroVol)) {

  if(length(unique(weigangHydroVol[[i]][,2]))!=21) {
    print(i)
  }
}

weigangHydroVol[[1]]
#check if all codes have 64 codons
# check if all paths contain 21 amino acids
for ( i in 1:length(weigangHydroVol)) {
  
  if(length(unique(weigangHydroVol[[i]][,1]))!=64) {
    print(i)
  }
}

length(weigangHydroVol)
weigangHydroVol[[1]]
#File names dataframe_*.txt  
naa<-c()
for (i in 1:length(completeCode)){
  naa[i]<-paste("dataframe",i,sep="_")
  naa[i]<-paste(naa[i],"txt",sep=".")
  naa[i]<-paste("C:/Users/brian/OneDrive/Documents/AANeighbors/2020Paths/2020PolVolIsoHydro",naa[i],sep="/")}
length(naa)

for (i in 1:length(completeCode)) {
  write.table(completeCode[[i]], naa[i], sep="\t")
}
nrow(fixXcodes[[1]])
