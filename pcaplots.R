#Plots  for paths, and Distance Density Plots with and without PCA coordinates.
#All dataframs contain codon paths determined by the path in weigangCode
#ggplot2 is used for the plots


#PART I -  Plots with Polarity v Volume Dimension.  No PCA

library(factoextra) # PCA library gets you coordinates

#Section I- Load code.  aadistance, minPath, weigangCode

# weigang's code  

weigangCode<-read.csv("C:/Users/brian/OneDrive/Documents/AANeighbors/weigangCodonList.csv", sep=" ",header=FALSE)
weigangCode<-data.frame(weigangCode[,13],weigangCode[,1])
ncol(weigangCode)
x<-order(weigangCode[,1])
weigangCode[order(weigangCode[,1]),]


#Amino Acid Coordinates Dataframe


aadistances<- read.csv("C:/Users/brian/OneDrive/Documents/MonteCarloClub/Neural Nets/TSP/aadistances.csv",header = TRUE)



#minPath- a dataframe with a tsp path optimized for polarity and Volume error.
minPath<- read.table("C:/Users/brian/OneDrive/Documents/AANeighbors/parameterTesting/twoParameters/WPolVol/Wtitvx4_1.txt",header = TRUE)



#Section II- Attaching Coordinates to Paths

# attachCoordinates- function to attach coordinates (Polarity and Volume) to MinPath

attachCoordinates<-function(data) {
  minPathCoordinate<-data
  polCol<-c()
  volCol<-c()
  for (i in 1:64) {
    x<-subset(aadistances, aadistances[,1]==data[i,2])
    
    polCol[i]<-as.numeric(x[2])
    volCol[i]<-as.numeric(x[4])
  }
  minPathCoordinate<-cbind(minPathCoordinate,polCol)
  minPathCoordinate<-cbind(minPathCoordinate,volCol)
  return(minPathCoordinate)
}

#minPathnew-tsp dataframe with coordinates. 
minPathnew<-attachCoordinates(minPath)

head(minPathnew)
#plot of amino acids
plot(minPathnew[,3],minPathnew[,4])

computeDistance<- function (data) {
  l<-data[,3:4]
  totaldistance<-0
  disk<-as.matrix(dist(l, method="euclidean"))
  for (i in 1:63) {
    totaldistance<- totaldistance+ disk[i,i+1]
  }
  return(totaldistance)
}

computeDistance(minPathnew)


#fileList- list of tsp file paths
fileList<-list()
for (i in 1:52){
  name<-paste("Wtitvx4",i,sep="_")
  name<-paste(name,'txt', sep=".")
  path<-paste("C:/Users/brian/OneDrive/Documents/AANeighbors/parameterTesting/twoParameters/WPolVol",name,sep="/")
  fileList[[i]]<-path
}

#Section II- Attaching Coordinates to Paths
fileList

# attachCoordinates- function to attach PCA coordinates  to tsp path (minPath)


#attachCoordinatesPCA- data-input file name
attachCoordinatesFiles<-function(data) {
  minPath<- read.table(data)
  minPathCoordinate<-minPath
  polCol<-c()
  volCol<-c()
  for (i in 1:64) {
    x<-subset(aadistances, aadistances[,1]==minPath[i,2])
    
    polCol[i]<-as.numeric(x[1])
    volCol[i]<-as.numeric(x[2])
  }
  minPathCoordinate<-cbind(minPathCoordinate,polCol)
  minPathCoordinate<-cbind(minPathCoordinate,volCol)
  return(minPathCoordinate)
}
length(fileList)

#minPathList- list of paths of AA with 2D PCA coodinate

minPathFiles<-lapply(fileList,attachCoordinatesFiles)
minPathlength<-lapply(minPathFiles,computeDistance)

minPathlengthvector<-c()
for (i in 1:52) {
  minPathlengthvector[i]<-minPathlength[[i]]
}

cor(exp(errordata[[2]]),log(minPathlengthvector))
cor((errordata[[3]]),minPathlengthvector)
lm1<-lm(log(minPathlengthvector)~exp(errordata[[2]])+exp(errordata[[3]])+exp(errordata[[1]]))

summary(lm1)
#Check to see if codon path of minPathnew is alligned with code of weigangCode

for (i in 1:64) {
  if (as.character(minPathnew[i,1])!=as.character(weigangCode[i,2])) {
    print("no")
  }
  
}

# Constructing the SGC dataframe- dd  .
#Note- this dataframe will be alligned with Weigang's Codon Path

codon<-c("ATT, ATC, ATA","CTT, CTC, CTA, CTG, TTA, TTG","GTT, GTC, GTA, GTG","TTT, TTC","ATG","TGT, TGC",
         "GCT, GCC, GCA, GCG","GGT, GGC, GGA, GGG","CCT, CCC, CCA, CCG","ACT, ACC, ACA, ACG","TCT, TCC, TCA, TCG, AGT, AGC",
         "TAT, TAC","TGG","CAA, CAG","AAT, AAC","CAT, CAC","GAA, GAG","GAT, GAC","AAA, AAG","CGT, CGC, CGA, CGG, AGA, AGG","TAA, TAG, TGA")

slc<-c("I","L","V","F","M","C","A","G","P","T","S","Y","W","Q","N","H","E","D","K","R","X")


codon.list<-strsplit(codon,", ")

codon.list

#dd- a  (codon,aa) data frame.  

codonaa<-data.frame(slc[[1]],codon.list[[1]])
dataframenames<-list()
for (i in 1:21) {
  dataframenames[[length(dataframenames)+1]]<-data.frame(slc[[i]],codon.list[[i]])
}

for (i in 1:21) {
  colnames(dataframenames[[i]])<-c("AA","Codons")
}
dataframenames[[2]]
dd<-dataframenames[[1]]
for (i in 1:20) {
  dd<-rbind(dd,dataframenames[[i+1]])
}



#Use Weigang's ordering to find the code.  

#SGCCoordinates-function whose output is dataframe with a path order (data=WeigangCode) and coordinates.

SGCCoordinates<-function(data) {
  polCol<-c()
  volCol<-c()
  aan<-data.frame(matrix(ncol = 2, nrow = 0))
  v <- c("AA", "Codons")
  colnames(aan) <- v
  
  
  for (i in 1:64) {                     #reorders dd according to WeigangCode
    
    x<-subset(dd, data[i,2]==dd[,2])
    aan<-rbind(aan,x)
  
  }
  
  
  for ( j in 1:64) {                   #attaches coordinates to the dataframe
    z<-subset(aadistances, aadistances[,1]==aan[j,1])
    polCol[j]<-as.numeric(z[2])
    volCol[j]<-as.numeric(z[4])
  }

  aan<-cbind(aan,polCol)
  aan<-cbind(aan,volCol)
  return(aan)
  
}


#sgcweigangCode-dataframe
sgcWeigangCode<-SGCCoordinates(weigangCode)




#check to see if order of sgcweigangCode df follows  
for (i in 1:64) {
  if (as.character(randomCode[i,2])!=as.character(weigangCode[i,2])) {
    print("no")
  }
}
tail(minPathnew)
#Generate Random code by shuffling blocks

# Use SGC dataframe, but shuffle slc  before you do the frame. 

codon<-c("ATT, ATC, ATA","CTT, CTC, CTA, CTG, TTA, TTG","GTT, GTC, GTA, GTG","TTT, TTC","ATG","TGT, TGC",
         "GCT, GCC, GCA, GCG","GGT, GGC, GGA, GGG","CCT, CCC, CCA, CCG","ACT, ACC, ACA, ACG","TCT, TCC, TCA, TCG, AGT, AGC",
         "TAT, TAC","TGG","CAA, CAG","AAT, AAC","CAT, CAC","GAA, GAG","GAT, GAC","AAA, AAG","CGT, CGC, CGA, CGG, AGA, AGG","TAA, TAG, TGA")

slc<-c("I","L","V","F","M","C","A","G","P","T","S","Y","W","Q","N","H","E","D","K","R","X")



slc<-sample(slc,21)  #ALERT: shuffles the blocks. Only use this to generate random codes
slc
codon.list<-strsplit(codon,", ")

codon.list

#make a codon aa data frame.  Start with a list of data frames the rbind elts of the list to the dataframe.

codonaa<-data.frame(slc[[1]],codon.list[[1]])
dataframenames<-list()
for (i in 1:21) {
  dataframenames[[length(dataframenames)+1]]<-data.frame(slc[[i]],codon.list[[i]])
}

for (i in 1:21) {
  colnames(dataframenames[[i]])<-c("AA","Codons")
}
dataframenames[[2]]
dd<-dataframenames[[1]]
for (i in 1:20) {
  dd<-rbind(dd,dataframenames[[i+1]])
}
dd

#randomCoordinates- function to add coordinates to random paths, that use WeigangCode, ordering
randomCoordinates<-function(data) {
  polCol<-c()
  volCol<-c()
  aan<-data.frame(matrix(ncol = 2, nrow = 0))
  v <- c("AA", "Codons")
  colnames(aan) <- v
  
  
  for (i in 1:64) {
    
    x<-subset(dd, data[i,2]==dd[,2])
    aan<-rbind(aan,x)
    
  }
  
  
  for ( j in 1:64) {
    z<-subset(aadistances, aadistances[,1]==aan[j,1])
    polCol[j]<-as.numeric(z[2])
    volCol[j]<-as.numeric(z[4])
  }
  
  aan<-cbind(aan,polCol)
  aan<-cbind(aan,volCol)
  return(aan)
  
}




#Section 4-- Making the ggplot

#Give everyone the same colnames
colnames(minPathnew)<-colnames(randomCode)
  
#Make a data frame using all three codes.  Mark each dataframe by the terms "SGC" 
randomCode<-randomCoordinates(weigangCode)
df1<-rbind(sgcWeigangCode,randomCode)
df1<-rbind(df1,minPathnew)
df1$type<-c(rep("SGC",64),rep("Random",64),rep("Optimized Code",64))
randomCode$type=c(rep("Amino Acids",64))

#use ggplot2 to plot
p<-ggplot(data=df1,aes(polCol,volCol))
qplot(polCol,volCol,data=df1) +geom_path(aes(colour=type),arrow=arrow(angle = 20, length = unit(0.15, "inches"), ends = "first", type = "closed"))+
facet_wrap(~type)+labs(x='Polarity',y='Volume')+
  #geom_point(data=subset(randomCode,type=="Amino Acids"),aes(polCol,volCol,label=AA))+
  geom_text(data=randomCode,aes(polCol,volCol,label=AA))


#PART II Do the same as partI but with PCA coodinates


#pcadf-  a dataframe (PCAdim1, PCAdim2,AA)


aadists<-aadistances[,2:5]
colnames(aadists)





res.pca <- prcomp(aadists)
res.ind <- get_pca_ind(res.pca)
res.ind$coord    # Coordinates for the 21 amino acids


pcadf<-data.frame(matrix(ncol=2,nrow=0))

pcadf<-data.frame(res.ind$coord[,1],res.ind$coord[,2])
colnames(pcadf)<-c("pcad1","pcad2")
pcadf$AA<-aadistances[,1]

plot(pcadf[,1],pcadf[,2])

#Section II -call the tsp files

#fileList- list of tsp file paths
fileList<-list()
for (i in 1:52){
  name<-paste("Wtitvx4",i,sep="_")
  name<-paste(name,'txt', sep=".")
  path<-paste("C:/Users/brian/OneDrive/Documents/AANeighbors/parameterTesting/twoParameters/WPolVol",name,sep="/")
  fileList[[i]]<-path
}

#Section II- Attaching Coordinates to Paths

# attachCoordinates- function to attach PCA coordinates  to tsp path (minPath)


#attachCoordinatesPCA- data-input file name
attachCoordinatesPCA<-function(data) {
  minPath<- read.table(data)
  minPathCoordinate<-minPath
  polCol<-c()
  volCol<-c()
  for (i in 1:64) {
    x<-subset(pcadf, pcadf[,3]==minPath[i,2])
    
    polCol[i]<-as.numeric(x[1])
    volCol[i]<-as.numeric(x[2])
  }
  minPathCoordinate<-cbind(minPathCoordinate,polCol)
  minPathCoordinate<-cbind(minPathCoordinate,volCol)
  return(minPathCoordinate)
}
length(fileList)

#minPathList- list of paths of AA with 2D PCA coodinate

minPathList<-lapply(fileList,attachCoordinatesPCA)



#computeDistance-function  input a path with coodinates get the length of path.  
computeDistance<- function (data) {
  l<-data[,3:4]
  totaldistance<-0
  disk<-as.matrix(dist(l, method="euclidean"))
  for (i in 1:63) {
    totaldistance<- totaldistance+ disk[i,i+1]
  }
  return(totaldistance)
}
tspDistances<-lapply(minPathList,computeDistance)
tspDist<-c()
for (i in 1:52) {
  tspDist[i]<-tspDistances[[i]]
}



#check that listing of codons alligns with Qiu's codon path.

minPathnewPCA<-attachCoordinatesPCA(minPath)
head(minPathnewPCA)

plot(minPathnew[,3],minPathnew[,4])
for (i in 1:64) {
  if (as.character(minPathnew[i,1])!=as.character(weigangCode[i,2])) {
    print("no")
  }
  
}


# SGC dataframe ordered according to QIU's codon path- 
#Section I- dd Construction- dd is the SGC represented in a dataframe (AA, Codons) 

codon<-c("ATT, ATC, ATA","CTT, CTC, CTA, CTG, TTA, TTG","GTT, GTC, GTA, GTG","TTT, TTC","ATG","TGT, TGC",
         "GCT, GCC, GCA, GCG","GGT, GGC, GGA, GGG","CCT, CCC, CCA, CCG","ACT, ACC, ACA, ACG","TCT, TCC, TCA, TCG, AGT, AGC",
         "TAT, TAC","TGG","CAA, CAG","AAT, AAC","CAT, CAC","GAA, GAG","GAT, GAC","AAA, AAG","CGT, CGC, CGA, CGG, AGA, AGG","TAA, TAG, TGA")

slc<-c("I","L","V","F","M","C","A","G","P","T","S","Y","W","Q","N","H","E","D","K","R","X")


codon.list<-strsplit(codon,", ")

codon.list

#make a codon aa data frame.  Start with a list of data frames the rbind elts of the list to the dataframe.

codonaa<-data.frame(slc[[1]],codon.list[[1]])
dataframenames<-list()
for (i in 1:21) {
  dataframenames[[length(dataframenames)+1]]<-data.frame(slc[[i]],codon.list[[i]])
}

for (i in 1:21) {
  colnames(dataframenames[[i]])<-c("AA","Codons")
}
dataframenames[[2]]
dd<-dataframenames[[1]]
for (i in 1:20) {
  dd<-rbind(dd,dataframenames[[i+1]])
}

#Section II
#Use Qiu's path to get a dataframe

# Q- code dataframe (AA number, Codons) 
weigangCode<-read.csv("C:/Users/brian/OneDrive/Documents/AANeighbors/weigangCodonList.csv", sep=" ",header=FALSE)
weigangCode<-data.frame(weigangCode[,13],weigangCode[,1])
ncol(weigangCode)
weigangCode
 

#SGCCoordinates-function to transform ordering of codons in SGC Dataframe into Weigang's path.

SGCCoordinatesPCA<-function(data) {
  polCol<-c()
  volCol<-c()
  aan<-data.frame(matrix(ncol = 2, nrow = 0))
  v <- c("AA", "Codons")
  colnames(aan) <- v
  
  
  for (i in 1:64) {
    
    x<-subset(dd, data[i,2]==dd[,2])
    aan<-rbind(aan,x)
    
  }
  
  
  for ( j in 1:64) {
    z<-subset(pcadf, pcadf[,3]==aan[j,1])
    polCol[j]<-as.numeric(z[1])
    volCol[j]<-as.numeric(z[2])
  }
  
  aan<-cbind(aan,polCol)
  aan<-cbind(aan,volCol)
  return(aan)
  
}

sgcWeigangCodePCA<-SGCCoordinatesPCA(weigangCode)
sgcWeigangCodePCA




for (i in 1:64) {
  if (as.character(randomCode[i,2])!=as.character(weigangCode[i,2])) {
    print("no")
  }
}

#Section 3

#Generate Random codes by shuffling blocks

#create dd -Use SGC dataframe, but shuffle slc  before you do the frame. 

slcPermList<-list()
for( i in 1:100) {
codon<-c("ATT, ATC, ATA","CTT, CTC, CTA, CTG, TTA, TTG","GTT, GTC, GTA, GTG","TTT, TTC","ATG","TGT, TGC",
         "GCT, GCC, GCA, GCG","GGT, GGC, GGA, GGG","CCT, CCC, CCA, CCG","ACT, ACC, ACA, ACG","TCT, TCC, TCA, TCG, AGT, AGC",
         "TAT, TAC","TGG","CAA, CAG","AAT, AAC","CAT, CAC","GAA, GAG","GAT, GAC","AAA, AAG","CGT, CGC, CGA, CGG, AGA, AGG","TAA, TAG, TGA")

slc<-c("I","L","V","F","M","C","A","G","P","T","S","Y","W","Q","N","H","E","D","K","R","X")
slc<-sample(slc,21)
codon.list<-strsplit(codon,", ")
#make a codon aa data frame.  Start with a list of data frames the rbind elts of the list to the dataframe.
codonaa<-data.frame(slc[[1]],codon.list[[1]])
dataframenames<-list()
for (i in 1:21) {
  dataframenames[[length(dataframenames)+1]]<-data.frame(slc[[i]],codon.list[[i]])
}
for (i in 1:21) {
  colnames(dataframenames[[i]])<-c("AA","Codons")
}
dd<-dataframenames[[1]]
for (i in 1:20) {
  dd<-rbind(dd,dataframenames[[i+1]])
}
slcPermList[[length(slcPermList)+1]]<-dd
}

#make a list of 100 permutations of AA-
slcPermList[[2]]




#Now use this dd to make a code
randomCoordinatesPCA<-function(data) {
  polCol<-c()
  volCol<-c()
  aan<-data.frame(matrix(ncol = 2, nrow = 0))
  v <- c("AA", "Codons")
  colnames(aan) <- v
  
  
  for (i in 1:64) {
    
    x<-subset(data, weigangCode[i,2]==data[,2])
    aan<-rbind(aan,x)
    
  }
  
  
  for ( j in 1:64) {
    z<-subset(pcadf, pcadf[,3]==aan[j,1])
    polCol[j]<-as.numeric(z[1])
    volCol[j]<-as.numeric(z[2])
  }
  
  aan<-cbind(aan,polCol)
  aan<-cbind(aan,volCol)
  return(aan)
  
}

randomCodeList<-lapply(slcPermList,randomCoordinatesPCA)
randomCodeList[[1]]
#Computing the length of random codes.
a<-randomCodeList[[1]][,3:4]
dis<-as.matrix(dist(a,method="euclidean"))
dis[1,4]

computeDistance<- function (data) {
  l<-data[,3:4]
  totaldistance<-0
  disk<-as.matrix(dist(l, method="euclidean"))
  for (i in 1:63) {
    totaldistance<- totaldistance+ disk[i,i+1]
  }
  return(totaldistance)
}
nrow(randomCodeList[[1]])
computeDistance(randomCodeList[[1]])

distanceofRandomPaths<-lapply(randomCodeList,computeDistance)
distanceofRandomPaths
randomDistance<-c()
for (i in 1:100) {
  randomDistance[i]<-distanceofRandomPaths[[i]]

}
plot(density(randomDistance))
#Section IV
#ggplot2
#Make a data frame using all three codes.

colnames(minPathnewPCA)<-colnames(randomCodePCA)
randomCodePCA<-randomCoordinatesPCA(weigangCode)
df1<-rbind(sgcWeigangCodePCA,randomCodePCA)
df1<-rbind(df1,minPathnewPCA)
df1$type<-c(rep("SGC",64),rep("Random Code",64),rep("Optimized Code",64))
head(df1)
randomCodePCA$type=c(rep("Amino Acids",64)) # randomcodePCA this will be used for geom_txt amino acid plot
head(randomCodePCA)


data.frame(errordata[[1]],errordata[[3]],tspDist)





#ggplot2 to plot

qplot(polCol,volCol,data=df1) +geom_path(aes(colour=type),arrow=arrow(angle = 20, length = unit(0.15, "inches"), ends = "first", type = "closed"))+
  facet_wrap(~type)+labs(x='PCA dim1',y='PCA dim2')+
  #geom_point(data=subset(randomCodePCA,type=="Amino Acids"),aes(polCol,volCol,label=AA))+
  geom_text(data=randomCodePCA,aes(polCol,volCol,label=AA))










#Section V -  Computing distances

#computeDistance (function(data))-input path with coodinates(nrow-64) Use lists of random and  tsp paths

computeDistance<- function (data) {
  l<-data[,3:4]
  totaldistance<-0
  disk<-as.matrix(dist(l, method="euclidean"))
  for (i in 1:63) {
    totaldistance<- totaldistance+ disk[i,i+1]
  }
  return(totaldistance)
}

#sgcdist- length of sgc path.
sgcdist<-computeDistance(sgcWeigangCodePCA)


#Section VIII computing distance density and histogram
#computing 
mean(tspDist)
mean(randomDistance)
sd(tspDist)
sd(randomDistance)

plot(density(tspDist))
plot(density(randomDistance))


#Section IX GGPLOT of Distance

#pathdf- dataframe with tspDistances, and Random Code distances.  

pathDist<-c(tspDist,randomDistance)
type<-c(rep("Optimized Code",52),rep("RandomCode",100))
pathdf<-data.frame(pathDist,type)
sgcFrame<-data.frame(sgcdist,0)


#ggplot code
p<-ggplot(aes(pathDist,colour=type),data=pathdf)
p+geom_density()+geom_vline(aes(xintercept=sgcdist),colour="red",data=sgcFrame)+
  labs(x='Path Distance',y='Frequency')+
  scale_color_manual(values=c("black","blue"))

