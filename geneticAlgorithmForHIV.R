#Goal-Using genetic algorithms we seek to evolve sets of HIV alleles, represented by binary sequences, so that minimum distance between a 
# any pair of alleles is maximized.   
#This question is relevant to immunologists looking for similarities between HIV strains that aren't related by recombinations. 
#It's also relevant to computer science folks looking for to maximize minimum distance between codewords with a prescribed length. 


#S1 Functions, and Priliminaries- population, fitness functions, mutation function

#population- A function to make initial population for genetic alg. The parameters are
#nloci- number of loci, nalleles-number of alleles, and popsize- population size.  
#An individual in this population is a matrix dim(nalleles*nloci) 
#In order to efficiently perform the genetic algorthm the matrix-individual- 
#is flattened into a vector of size-nalleles*nloci


population<-function(nloci,nalleles,popsize) {
  
  pop<-list()
  for (j in 1:popsize) {
    a<-sample(c(0,1),nloci*nalleles,replace = TRUE)
     pop[[length(pop)+1]]<-a
      }
  return(pop)
}


#"matdis" fitness function 1- outputs the hamming distance between each pair of alleles.
matdis<-function(data,nloci,nalleles) {
  suma<-c()
  for (i in 0:(nalleles-2)){
    for (j in (i+1):(nalleles-1)) {
      a<-c()
      for (k in 1:nloci) {
       a[length(a)+1]<-abs(data[nloci*i+k]- data[nloci*j+k])
      }
      
      suma[length(suma)+1]<-sum(a)
    }
  }
  return(suma)
}


#mutate-function causes an individuals alleles to change, with a given probability of .2
mutate<-function(data,nloci,nalleles) {
  for (j in 1:nalleles*nloci) {          #number of positions to mutate
    coin<-rbinom(1,1,.2)
    if (coin==1){
      if(data[j]==0) {
        data[j]<-1 }
      else {data[j]<-0}
    }
  }
  return(data)
}



#S2  TWO GENETIC ALGORITHMS
#genalg<- The genetic algorithm function.  Use population list as input.
##score and sumscore outputs give a mode for determining fitness by 
#1) taking individuals with sum of pw differences  (scores)
#or 2) taking individuals that maximize the minimum differences between scores (sumScore) 

 
genAlgTotalPWDiff<-function(data,nloci,nalleles,generations) {

avgdist<-c()
maxdist<-c()
mindist<-c()
p4<-data
for (i in 1:generations) {    #loop on generations
  
  p4<-lapply(p4,mutate,nloci,nalleles)      #mutate
  score<-lapply(p4,matdis,nloci,nalleles)  #score-list of pairwise differences btwn alleles for each individual 
  
  sumScore<-lapply(score,sum)  #sum score-  sums pw differences btwn alleles (scores) for each individual.
  
  sumScore<-unlist(sumScore)
  
  #we collect this data, which is based on the goal #1 of the function. 
  avgdist<-append(avgdist,mean(unlist(sumScore)))
  maxdist<-append(maxdist,max(unlist(sumScore)))
  mindist<-append(mindist,min(unlist(sumScore)))
  
  #fitness function determined by sumScore- which is the total pw difference of an individual. 
 #scoreVector-a vector that weights fit individuals based on total pw difference btwn pairs of alleles
  
  scoreVector<-c()
  for (i in 1:length(score)) {
    if(sumScore[i]==0) {
      next
    }
    scoreVector<-append(scoreVector,c(rep(i,sumScore[i])))
  }
  
  sc<-sample(scoreVector,length(p4))    #sample,  without replacement, from the scorevector to get number of copies of individuals   
  newList<-list()
  for (j in 1:length(sc)) {
    newList[[length(newList)+1]]<-p4[[sc[j]]]
  }
  p4<-newList  
  
}
 
 resultsdf<-data.frame(avgdist,maxdist,mindist)
  return(resultsdf)
}



#genAlgMaxMin- the fitness function is derived by taking the min of pw diff btwn alleles of an individual

scoreGTZero<-function(data) {
  positiveSet<-subset(data,data>0)
  return(positiveSet)
}


genAlgMaxMin<-function(data,nloci,nalleles,generations) {
  keys<-c(1:length(data))
  
  maxmindist<-c()
  maxdist<-c()
  p4<-data
  for (i in 1:generations) {    #loop on generations
    
    p4<-lapply(p4,mutate,nloci,nalleles)      #mutate
    score<-lapply(p4,matdis,nloci,nalleles)  #score-list of pairwise differences btwn alleles for each individual 
    scoreNotZero<-lapply(score,scoreGTZero)  #take out zero pw differences btwn alleles for individuals
    sumScore<-lapply(score,sum)  
  #make a df for the minScores gt zero
    minScore<-lapply(scoreNotZero,min)  #min score-takes min pw differences btwn alleles (scores) for each individual.
    minScore<-unlist(minScore)
    
    minScoredf<-data.frame(keys,minScore)
    minScoredf<-minScoredf[order(minScoredf[,2],decreasing = TRUE),]
    minScoredf<-minScoredf[1:floor(length(data)/2),]
  
    #we collect this data, which is based on the goal #1 of the function. 
    maxmindist<-append(maxmindist,max(minScore))
    maxdist<-append(maxdist,max(unlist(sumScore)))
    
    #fitness function determined by sumScore- which is the total pw difference of an individual. 
    #scoreVector-a vector that weights fit individuals based on total pw difference btwn pairs of alleles
    
    scoreVector<-c()
    for (i in 1:nrow(minScoredf)) {
      if(minScoredf[i,2]==0) {
        next
      }
      scoreVector<-append(scoreVector,c(rep(minScoredf[i,1],nalleles*minScoredf[i,2])))
    }
    
    sc<-sample(scoreVector,length(p4))    #sample,  without replacement, from the scorevector to get number of copies of individuals   
    newList<-list()
    for (j in 1:length(sc)) {
      newList[[length(newList)+1]]<-p4[[sc[j]]]
    }
    p4<-newList  
    
  }
  
  resultsdf<-data.frame(maxmindist,maxdist)
  return(resultsdf)
}


#genAlgMaxMin Test  - 
#Hold nalleles fixed and increase number of loci.
#output the max of the set of total distances 
#and the max of the minimum pw distances between alleles

maxmindistance<-c()
max2distance<-c()
lociSeq<-seq(10,30,10)
lociSeq
for (i in lociSeq)  {
  p1<-population(i,16,100)
  p1Evolve<-genAlgMaxMin(p1,i,16,10) 
  print(p1Evolve)
  maxmindistance[length(maxmindistance)+1]<-max(p1Evolve[,1])
  max2distance[length(max2distance)+1]<-max(p1Evolve[,2])
  
  }


plot(maxmindistance)
plot(max2distance)

#genAlgTotalPWDistance Test 
#Hold nalleles fixed and increase number of loci.
#output the max of the set of total distances 

maxdistance<-c()
lociSeq<-seq(10,30,10)

for (i in lociSeq)  {
  p2<-population(i,16,100)
  p2Evolve<-genAlgTotalPWDiff(p2,i,16,10) 
  print(p2Evolve)
  maxdistance[length(maxdistance)+1]<-max(p2Evolve[,2])
  
}


maxdistance

#Compare output of genetic algorithms

plot(maxdistance,max2distance)


