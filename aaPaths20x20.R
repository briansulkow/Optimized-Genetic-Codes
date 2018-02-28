##TSP Paths, via Hopfield Networks,   through  with No Stop Codon:
#  are using the Amino Acids(no stop codons) 20x20 case for Weigang:  take out "X" row from  aadistance   

#Goal: USE Traveling Salesman Problem algorithm To find a path through the amino acids using two of more of the following coodinates as coordinates (polarity, hydropathy, volume, isoelectricity).

#The TSP algorithm is found in Miller Bout paper. 


#Sections of code are organized by Roman Numerals.  Subsections of code are organized by Capital Letters


####Section I- Data Wrangling- Construct a distance Matrix Using W Qiu's normalized amino acid measurements.  

###I-A
##Import Dr Qiu's normalized amino acid measurements (aadistances). Make matrix (aadist) that includes columns of aadistances to be optimized, ex (pol, vol).


aadistances<- read.csv("C:/Users/brian/Onedrive/Documents/MonteCarloClub/Neural Nets/TSP/aadistances.csv",header = TRUE)
aadistances
aadistances<-aadistances[-1,]
##boxplot of 4 dimensions.  Refer to 12/2/17 for discussion. 
boxplot(aadistances[,-1], xlabs=colnames(aadistances[,-1]), main="Boxplot of Four Types of Amino Acid Measures")

## aadist- new df containing columns of aadistances you want to optimize.
aadist<-aadistances[,-1]
aadist[1,]


aadist<-aadist[,-4]
aadist<-aadist[,-3]
aadist<-aadist[,-2]
colnames(aadist)
aadist

###I-B Make the Distance matrix (NormaaMatrix) for TSP hopfield network

## aaMatrix- Use aadist for distances
aaMatrix<-dist(aadist, method="euclidean")


##NormaaMatrix-The matrix used for analysis. max value=1.

NormaaMatrix<-aaMatrix/max(aaMatrix)
NormaaMatrix<-as.matrix(NormaaMatrix)

#STOP Section I

####Section II
###In section II we make matrices that will be used for tsp.  This approach, of making matrices, differs from Olivers approach but is equivalent. 


###Section (II-A)
##Filter cols- This matrix is associated with summation that has the dp parameter in the Exi and energy terms of the algorithm.

#Size of the matrix is 21*21. This is the size of the set of all pairs (amino acid, time). 

filtercols<-matrix(c(rep(0,(20*20)^2)),nrow=20*20)
for (i in 1:20){
  for (j in 1:20){
    for(k in 1:20) {
      if(i!=k) { 
        filtercols[[20*(i-1)+j,20*(k-1)+j]]<-1
      }
    }
  }
}

#check of values of matrix. Should be 1's
filtercols[1,]
for(i in 1:20) {
  print(filtercols[22,20*(i-1)+2])
}



###Section (II-B)

##distpen-This Matrix corresponds to the second summation in the EXi and Energy equations. It is the term where distance appears.   
#(I'll provide a lower dimensional example of this matrix.  for you and Oliver to look at in the next day or two.)



distpen<-matrix(c(rep(0,(20*20)^2)),nrow=20*20)

#the first and last for loops in this construction ensure that you won't have a cyclic path, rather there will be distinct start and end points. 

for(i in 1:20) {
  for(k in 1:20) {
    distpen[[20*(i-1)+1,20*(k-1)+2]]<-NormaaMatrix[i,k]
  }
}

for(i in 1:20) {
  for(l in 2:19) {
    for(k in 1:20) {
      
      distpen[[20*(i-1)+l,20*(k-1)+(l-1)]]<-NormaaMatrix[i,k]
      distpen[[20*(i-1)+l,20*(k-1)+l+1]]<-NormaaMatrix[i,k]
    }
  }
}

for(i in 1:20) {
  for(k in 1:20) {
    distpen[[20*(i-1)+20,20*(k-1)+19]]<-NormaaMatrix[i,k]
  }
}

##STOP Section II


####Section III


###Section III-A
## Initial Conditions

#V-initial condition- The initial condition is a vector V of size 21*21, one entry for each pair (amino acid, time). Hopfield and Tank ('85 pg147) recommend that the initial condition have random values close to .5

V<-c(rep(1,20*20))
for (i in 1:(20*20)) {
  V[i]<-.5+runif(1,-.1,.1)  #random entry around .5
}
length(V)
V



### Section III-B
## Vectors and Dataframes to keep track of Energy, valid paths.  
energytemps<-c()
voltageout<-data.frame (0,0,0)
voltageout
colnames(voltageout)<-c("temperature","dparameter","Valid")
rvoltsls<-list()

###Section III-C
##THis section give a program for the algorithm by Miller, and van der Bout p 131.  ## Some comments on the program will refer to specific lines from the Miller, Van der Bout algorithm. Those  comments will start with the word "LINE".

#Note- there are TWO outer loops. These are not part of the algorithm on pg 131.  The  loop using "m" let's me run the algorithm for various values of the dp parameter (most outerloop). The loop using "k" lets me run the algorithm for various values of T. That's related to the annealing.  Rather than taking the temperature down exponentially as Oliver's algorithm does, the temperature is taken down linearly.  
#Note- Since I found values that worked well, I held dp and T constant in the program below . 
# With dp=.7 and T=.1 I was able to get 52 valid paths from 1500 runs (dim's were pol and volume). If you run this algorithm, try 50 runs.  You'll get 1 or 2 valid paths.

for (m in 0:50)
{
  for (k in 1:1) {
    T<-.1 #-k*1e-1
    dp<-.7 # .1*(m-1)
    voltvector<-V     #voltvector =initial condition V
    energyresults<-list() #grabs info on energy.
    j<-700    #parameter of while loop given below  
    samplepoint<-c()   # vector to catch the amino acids that were randomly sampled used.
    while(j!=0) {  #LINE- the while loop is used for "do until (a fixed point is found)"
      samplecodon<-sample(c(1:20),1) #LINE-"Select a city at random"
      samplepoint[length(samplepoint)+1]<-samplecodon
      meanfield<-c(rep(1,20))
      sumterms<-0
      for(i in 1:20){  #LINE-for (i in i<=n)..
        meanfield[i]<-dp*filtercols[20*(samplecodon-1)+i,]%*%voltvector +
          distpen[20*(samplecodon-1)+i,]%*%voltvector #LINE- EXi equation.
        sumterms<-sumterms+exp(-meanfield[i]/T)  #LINE- "sum<=sum +e^-EXi/T"
        
      }
      
      for(i in 1:20) { ##LINE=for(i<=1;...) 
        voltvector[20*(samplecodon-1)+i]<- exp(-meanfield[i]/T)/sumterms
      }  #LINE-     "VXi=exp^(-EXi/T)/sum" This loop adjusts the voltvector, changing values corresponding to city=X time=i for all i
      #end loop  
      
      ##sumtest is meant to break the loop if any voltage goes to infinity.
      sumtest<-c(rep(1,20))
      for(i in 1:20){
        sumtest[i]<-voltvector[20*(samplecodon-1)+i]
      }
      
      st<-sum(sumtest)  
      if (is.nan(st)==TRUE) {
        print("bad")
        break}
      #LINE- "E=dp/2...."  THis is the energy equation  
      energy<-(dp/2)*voltvector%*%filtercols%*%voltvector+
        .5*voltvector%*%distpen%*%voltvector
      
      energyresults[[length(energyresults)+1]]<-energy
      j<-j-1
      
    }
    #End While Loop
    
    
    #the next few lines of the program is help guard against non convergence.  
    energytemps[length(energytemps)+1]<-energy
    if(is.nan(st)==TRUE) 
    { voltageout<-rbind(voltageout,c(T,dp,"NaN"))
    
    rvoltsls[[length(rvoltsls)+1]]<-"noconvergence" 
    next }
    volts<-t(matrix(voltvector,nrow=20))  # volts- is the matrix form of  voltvector. Notice we have to take the transpose. 
    
    #rvolts will create our matrix of 1;s and 0;s using volts matrix.  
    rvolts<-volts
    for(p in 1:20) {
      for (q in 1:20) {
        if(rvolts[p,q]==max(volts[p,])) 
        {rvolts[p,q]<-1} else {rvolts[p,q]<-0}
      }}
    
    #Useful output-
    #rvoltsls-a list of rvolts matrices.          #voltage- a dataframe that's used in processing the output. 
    rowvolts<-c()   
    
    
    for (i in 1:20) {
      rowvolts[length(rowvolts)+1]<-sum(rvolts[i,])
    }
    sumrv<-sum(rowvolts)
    
    if(sum(rowvolts)==20)
    {voltageout<-rbind(voltageout,c(T,dp,sumrv))
    rvoltsls[[length(rvoltsls)+1]]<-rvolts
    }
    if(sum(rowvolts)<20)
    {voltageout<-rbind(voltageout,c(T,dp,sumrv))
    rvoltsls[[length(rvoltsls)+1]]<-rvolts
    
    }
  } 
}



#END PROGRAM






plot(energytemps)
####Section IV
###FINDING ALL VALID PATHS AND GETTING LISTS OF MATRICES OF VALID PATHS

####Section IV-A
##We first need to get column and rowsums.
#rvoltsround- a list of the rounded matrices of rvoltsls.This can be stream lines
rvoltsround<-lapply(rvoltsls,round)
length(rvoltsround)


#Column sums and Row sums of outputs.
sumsrvoltsround<-lapply(rvoltsround,colSums)
sumsrvoltsround
rowsumsrvoltsround<-lapply(rvoltsround,rowSums)
rowsumsrvoltsround

###Section IV-B
##Finding valid paths.
#roundvector- tells you which elements of rvoltsround represent valid paths.
#goodvector-a list of matrices representing good paths.


roundvector<-c()          
#roundvectors find valid paths
for (i in 1:length(sumsrvoltsround)) {
  transformsum<-c()
  for (j in 1:20) {
    transformsum[j]<-sumsrvoltsround[[i]][j]/sumsrvoltsround[[i]][j]
  }
  for (k in 1:20) {
    if (is.nan(transformsum[k])) {
      transformsum[k]<-0
    }
  }
  if (sum(rowsumsrvoltsround[[i]])==20 & sum(transformsum)==20) {
    roundvector[length(roundvector)+1]<- i
  }
}
roundvector
length(roundvector)

goodvolts<-list()
for (i in 1:length(roundvector)) {
  goodvolts[[length(goodvolts)+1]]<- rvoltsls[[roundvector[i]]]
}

roundvolts<-lapply(goodvolts,round)
length(roundvolts)
length(roundvolts[[1]][1,])

aamapsreOrder[[1]]


####Section V
###Output Paths in the form (AminoAcid Names, Paths) 
####NOTE YOU MUST RUN FUNCTION IN SECTION VI BEFORE YOU CAN GO THROUGH THIS SECTION.


AminoMaps<-lapply(roundvolts,maps)
AminoMaps[[1]]

length(AminoMaps)


aamapsreOrder<-lapply(AminoMaps,reorderAA)
aamapsreOrder[[1]]

#InsertX= function  inserts X into the paths

insertXInAamapsreorder<-function(data) {
  row8<-c(0,0,'X')
  datar<-rbind(data[1:7,],row8,data[8:20,])
  return(datar)
}

fullaaList<-lapply(aamapsreOrder,insertXInAamapsreorder)

fullaaList[[1]]

#AAWeigangCOdedf- Maps Codon Path (weigangCode) to AA Path. Input a full path through the AA;s 
# that includes "X"

AAWeigangCodedf<-function(data) {
  AminoAcids<-c()
  for (i in 1:64) {
    AminoAcids[i]<-data[weigangCode[i,1],3]
  } 
  df<-data.frame(weigangCode[,2],AminoAcids)
  return(df)
}
weigangCode

#completeCode- List of (AA,Codon) paths. Ready to save and put through Code-Stats.pl

completeCode<-lapply(fullaaList,AAWeigangCodedf)
weigangCode
completeCode[[1]]

#Write Paths-put complete codes into a path. 
#File names dataframe_*.txt  
naa<-c()
for (i in 1:length(completeCode)){
  naa[i]<-paste("dataframe",i,sep="_")
  naa[i]<-paste(naa[i],"txt",sep=".")
  naa[i]<-paste("C:/Users/brian/OneDrive/Documents/AANeighbors/2020Paths/2020Pol4",naa[i],sep="/")}
length(naa)

for (i in 1:length(completeCode)) {
  write.table(completeCode[[i]], naa[i], sep="\t")
}
###SECTION VI
##IGNORE, GO TO SECTION VII


####SECTION VI
#FUNCTIONS FOR PROCESSING OUTPUT INTO 21 PAIRS (AMINO ACID NAMES, TIME)
#RUN this section before running section V


#maps-Outputs 20 pairs, (Amino Acid Number,Time) .  Notice we get the amino acid number.  The next path with get the pairs.  

maps<-function(data) {
  #the matrix is 21x21  we want to 
  
  
  mapper<-data.frame()
  
  for (l in 1:20){
    for (m in 1:20){ 
      if (data[m,l]!=0 )
      {mapper<-rbind(mapper,c(l,m))}
    }
  }
  colnames(mapper)<-c("AA","codonnum")
  return(mapper)
  
}



#make an  AA names vector (length=20), without "STOP" AA

#change the factor levels of aadistances. THis is a necessary step to avoid confusion. 
aadistances<-aadistances[-1,]
newAALevels<-factor(aadistances[,1], levels=(aadistances[,1]))
length(newAALevels)

###reorderAA- Takes the output of maps (aminoacid number, time) and outputs (aminoacid NAME, time)  
reorderAA<-function(data) {
  AAPerm<-c()
  for (i  in 1:20){
    AAPerm[i]<- as.character(newAALevels[data[i,2]])
    
  }
  data[,3]<-AAPerm
  return(data)
}

#weigangCode - dataframe (block number,codon).
#this is path through the codons, that  includes sgc block structure, that we will map to the Amino Acids.
weigangCode<-read.csv("C:/Users/brian/OneDrive/Documents/AANeighbors/weigangCodonList.csv", sep=" ",header=FALSE)
weigangCode<-data.frame(weigangCode[,13],weigangCode[,1])
ncol(weigangCode)
x<-order(weigangCode[,1])
weigangCode[order(weigangCode[,1]),]
weigangCode






