###Goal = Statistics and Graphs on Error Results-

#Note each list has errors computed 4 different ways-(polarity, hydro, vol, iso)

#Import sgc error values, and make a list.  
# sgcerror<- error values of sgc are in  the following order(polarity, hydro, vol, iso)
sgc<- read.csv("C:/Users/brian//OneDrive/Documents/MonteCarloClub/Neural Nets/tspCodonResults/graphingResults/sgcpolhydovoliso.csv",header = TRUE)
ncol(sgc)
sgc
#list of error values
sgcerror<-list()
for (i in 1:4) {
  sgcerror[[i]]<-sgc[1,i]
}
sgcerror

plot(density(aadistances[,2]))
boxplot(aadistances[,2:4])
gaerrordata<-list()
for (i in 1:4) { 
  a<-paste("gbinerror",i,sep="_")
  b<-paste(a,"csv",sep=".")
  c<-paste("C:/Users/brian/OneDrive/Documents/AANeighbors/ParameterHydro",b,sep="/")
  
  gaerrordata[[i]]<-read.csv(c,header = FALSE)
  gaerrordata[[i]]<-gaerrordata[[i]][,1]
}

gaerrordata

gaerrors<-c()
for (i in 1:4) {
  gaerrors[i]<-gaerrordata[[i]][1]
}

sgcerrors<-c()
for (i in 1:4) {
  sgcerrors[i]<-sgcerror[[i]][1]
}
sgcerrors
gaerrors


errordata<-list()
for (i in 1:4) { 
  a<-paste("binerror",i,sep="_")
  b<-paste(a,"csv",sep=".")
  c<-paste("C:/Users/brian/OneDrive/Documents/AANeighbors/2020Paths/2020Pol4",b,sep="/")
  errordata[[i]]<-read.csv(c,header = FALSE)
  errordata[[i]]<-errordata[[i]][,1]
}

min(errordata[[1]])
errordata[[1]]



####Finding paths that beat sgc in all dimensions. 

optPathList<-list()
for (i in 1:4){
    optPathList[[i]]<-which(errordata[[i]]< sgcerror[[i]])
}
optPathList  
a<-intersect(optPathList[[1]],optPathList[[2]])  #sswitch second to 2
a
b<-intersect(optPathList[[3]],a)
b
paths<-intersect(optPathList[[4]],b)
paths
optPathValues<-list()

for ( i in 1:length(paths)) {
  errs<-c()
  for ( j in 1:4) {
    errs[j]<-errordata[[j]][paths[i]]
  }
  print(errs)
  optPathValues[[i]]<-errs
}
length(optPathValues)
optPathValues   

sgcerror  
  perrs<-c()
  for (i in 1:4) {
    perrs[i]<-perrordata[[i]][1]
  }
8/143  
perrs
  #import 4  random shuffle error files
#serrordata<-list of error vectors from random shuffle. Order is as follows (polarity, hydro, vol, iso)
length(errordata)

serrordata<-list()
for (i in 1:4) { 
  a<-paste("nserror",i,sep="_")
  b<-paste(a,"csv",sep=".")
  c<-paste("C:/Users/brian/OneDrive/Documents/MonteCarloClub/Neural Nets/tspCodonResults",b,sep="/")
  
  serrordata[[i]]<-read.csv(c,header = FALSE)
  serrordata[[i]]<-serrordata[[i]][,1]
}
length(serrordata[[4]])
#Scaling Error Data -lists are center...errordata


"Computing ratios of errors"

spolhydroslopes<-c()
for (i in 1:length(errordata[[2]])) {
  spolhydroslopes[i]<-abs(serrordata[[1]][i]/serrordata[[2]][i])
  
}
hist(spolhydroslopes)

spolvolslopes<-c()
for (i in 1:length(serrordata[[2]])) {
  spolvolslopes[i]<-abs(serrordata[[1]][i]/serrordata[[3]][i])
  
}
par(mfrow=c(1,1))


plot(density(spolhydroslopes),main="Pol/Hydro Slopes")
lines(density(polhydroslopes),main="pol/hydro slopes", col="red")
plot(density(spolvolslopes),main="Pol/Vol Slopes")
lines(density(polvolslopes),col="red")

#Center Data


centerserrordata<-list()
for (i in 1:4)  {
a<-c(rep(mean(serrordata[[i]]),length(serrordata)))
centerserrordata[[i]]<-(serrordata[[i]]-a)/(sd(serrordata[[i]]))
}

centererrordata<-list()
for ( i in 1:4) {
  a<-c(rep(mean(serrordata[[i]]),length(errordata[[i]])))
  centererrordata[[i]]<-(1/sd(serrordata[[i]]))*(errordata[[i]]-a)
}
length((centererrordata[[1]]))

#computing stats on error data and random error data

df.sd<-data.frame("pol","hydro","vol","iso")
sderrvector<-c()
for (i in 1:4){
  sderrvector[i]<-sd(errordata[[i]])
}

sdserrvector<-c()
for (i in 1:4){
  sdserrvector[i]<-sd(serrordata[[i]])
}
sderrvector
sdserrvector  
df.sd<-rbind(sderrvector,sdserrvector)

df.sd




#TSP DENSITY  derrordata<- This list has the density associate with each density of error data set in the list errordata
derrordata<-lapply(errordata,density)

#scaled density
dcentererrordata<-lapply(centererrordata,density)


#Random shuffle densities, dserrordata<- This list has the density associate with each error random shuffle data set in the list serrordata
dserrordata<-lapply(serrordata,density)

#scaled density
dcenterserrordata<-lapply(centerserrordata,density)

#Comparing density plots

#names of plots
testlist<-list("Polarity Error Distributions", "Hydrophobicity Error Distributions","Volume Error Distributions", "Isoelectric Error Distribution")
par(mfrow=c(2,2))
randomTestError<-c("TSP","Random")
#Ploting for SCALED Data
for ( i in 1:4)  {
plot(dcentererrordata[[i]], pch=19, col="red",xlim=c(9/10*min(dcentererrordata[[i]][[1]],dcenterserrordata[[i]][[1]]),max(dcenterserrordata[[i]][[1]],dcentererrordata[[i]][[1]],((sgcerror[[i]]-mean(serrordata[[i]]))/sd(serrordata[[i]])))*12/10),ylim=c(0,12/10*max(max(dcentererrordata[[i]][[2]]),max(dcenterserrordata[[i]][[2]]))), xlab="Error", ylab="Freq",main=testlist[[i]])
# Add a line
lines(dcenterserrordata[[i]], pch=18, col="blue", lty=2)
abline(v=((sgcerror[[i]]-mean(serrordata[[i]]))/sd(serrordata[[i]])),col="green") #sgc line
# Add a legend
legend(.80*max(max(dcentererrordata[[i]][[1]]),10.5/10*max(dcenterserrordata[[i]][[1]])),5/6*max(max(dcentererrordata[[i]][[2]]),max(dcenterserrordata[[i]][[2]])), legend=c("TSP","Random","SGC"),
       col=c("red", "blue","green"),lty=1:3, cex=.5)
}
#Ploting for Data
for ( i in 1:4)  {
  plot(derrordata[[i]], pch=19, col="red",xlim=c(9/10*min(derrordata[[i]][[1]],dserrordata[[i]][[1]]),max(dserrordata[[i]][[1]],derrordata[[i]][[1]])*12/10),ylim=c(0,12/10*max(max(derrordata[[i]][[2]]),max(dserrordata[[i]][[2]]))), xlab="Error", ylab="Freq",main=testlist[[i]])
  # Add a line
  lines(dserrordata[[i]], pch=18, col="blue", lty=2)
  abline(v=sgcerror[[i]],col="green") #sgc line
  # Add a legend
  legend(.80*max(max(dcentererrordata[[i]][[1]]),10.5/10*max(dcenterserrordata[[i]][[1]])),5/6*max(max(dcentererrordata[[i]][[2]]),max(dcenterserrordata[[i]][[2]])), legend=c("TSP","Random","SGC"),
         col=c("red", "blue","green"),lty=1:3, cex=.5)
}
?subset
#Where does the sgc error fit into the empirical distribution of each data set?   
#sgc.quantile.df  is a data frame in which the quantile for sgc error is found relative to the type of measurement and associated tsp and random empirical distributions.  
errordata
measuretype<-c("pol","hydro","vol","iso")
tspquant<-c()
randomquant<-c()
for (i in 1:4) {
  tspquant[i]<-ecdf(errordata[[i]])(sgcerror[[i]])
  randomquant[i]<-ecdf(serrordata[[i]])(sgcerror[[i]])
}
sgc.quantile.df<-data.frame(measuretype,tspquant,randomquant)
sgc.quantile.df

# Are the random distributions and TSP distributions the same?
# Mann Whitney test for distributions (the non parametric analogue of the t-test say "no".)  
mannwhitney<-list()

for (j in 1:4) {
  mannwhitney[[j]]<-wilcox.test(errordata[[j]],serrordata[[j]],alternative="less", conf.int = .95) # mann whitney test on two group
}

mannwhitney
p.values<-c()
NullHyp<-c(rep("reject",4))

for (i in 1:4) {
p.values[i]<-mannwhitney[[i]][[3]]
}
mannwhitney.df<-data.frame(p.values,NullHyp)

mannwhitney.df

#location difference
locationDiff<-list()
for (i in 1:4){
  locationDiff[[i]]<-mannwhitney[[i]][[9]]
}

#pvalue and location different list
MannWhitneyResults<-list()
for (i in 1:4){
  a<-c(locationDiff[[i]],mannwhitney.df[i,1])
  MannWhitneyResults[[i]]<-a
}
MannWhitneyResults

### engineering approach function
engineerApproach(errordata,serrordata)



#Engineering Approach function




engineerApproach<- function(mydata,randomdata) {
  minerror<-list()

  for (i in 1:4) {
  minerror[[i]]<-min(mydata[[i]])
  }

  meanserror<-list()
  for (i in 1:4) {
    meanserror[[i]]<-mean(randomdata[[i]])
    }

  engineApproach<-list()
  for (i in 1:4) {
    engineApproach[[i]]<-(meanserror[[i]]-sgcerror[[i]])/(meanserror[[i]]-minerror[[i]])
  }
  return(engineApproach)
}

mean(aadistances[,2])
sd(aadistances[,2])
boxplot(aadistances[,2],aadistances[,3],aadistances[,4],aadistances[,5])
?par
par(mfrow=c(2,2))
for (i in 2:5) {
  hist(aadistances[,i])
}

