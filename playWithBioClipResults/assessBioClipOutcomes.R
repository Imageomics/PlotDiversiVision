setwd("~/Documents/research/imageomics/floraPalooza/playWithBioClipResults")
rm(list=ls())

#Text data packages 
#install.packages("tm", dependencies = TRUE)
library(tm)

#**** CHANGE FILE NAMES
#**filename for predictions
#filePred<-"CPER_009_neon_grid3.csv"
#filePred<-"SCBI_005_neon_grid3.csv"
filePred<-"SCBI_005_neon_grid2.csv"
filePred<-"SCBI_005_neon_grid4.csv"



#**filename for species labels
#filePredLabels<-"CPER_NEON_labels.csv"
filePredLabels<-"SCBI_plot_labels.csv"

#**filename for ground truth
#fileTrueLabels<-"CPER_009_subplot_labels.csv"
fileTrueLabels<-"SCBI_005_subplot_labels.csv"



#**predictions from BioClip
pathPred<-paste("../PlotDiversiVision/temp_results/", filePred, sep="")
pred<-read.csv(pathPred, sep=",", header=TRUE)
matrix(names(pred))
names(pred) <- gsub("\\.", " ", names(pred)) #change periods in column headers to spaces
names(pred)<-tolower(names(pred))
names(pred)
#define which columns in pred are meta-data and whch columns are probabilities for each speaces
metaCols<-1:which(names(pred)=="species_count")
probCols<-(which(names(pred)=="species_count")+1):dim(pred)[[2]]


#species labesl 
pathPredLabels<-paste("../PlotDiversiVision/assets/species_list/", filePredLabels,sep="")
labelPred<-read.csv(pathPredLabels, sep=",", header=TRUE)
dim(labelPred)
names(labelPred)
labs<-tolower(labelPred$resolved_labels)
N<-length(labs) #number of species in the region


#**Ground truth; i.e., species listed by NEON technician
pathTrueLabels<-paste("../PlotDiversiVision/assets/test_labels/", fileTrueLabels,sep="")
gt<-read.csv(pathTrueLabels, sep=",", header=TRUE) #look at column: resolved_labels
dim(gt)
head(gt)
names(gt)
#This is a "short fat"  dataset.  I need to make this tall and skinny.  So, one row for every named species.  This plot ID repeats. 

sepTruths<-function(cell){
  cell<-tolower(cell)
  cell<-strsplit(cell,";")
  return(unlist(cell))
}
sepTruths(gt$resolved_labels[1])

#create data for first row in gt
gt2<-cbind(plotID=gt$plotID[1], subplotID=gt$subplotID[1], 
           trueCnt=gt$label_count[1], 
           trueSpecies=sepTruths(gt$resolved_labels[1]))

#this can get more efficient
for (i in 2:dim(gt)[[1]]){
  gt2<-as.data.frame(rbind(gt2,cbind(plotID=gt$plotID[i], 
                   subplotID=gt$subplotID[i], 
                   trueCnt=gt$label_count[i], 
                   trueSpecies=sepTruths(gt$resolved_labels[i]))))
}
gt2$trueCnt<-as.numeric(gt2$trueCnt)
head(gt2)

#**** NOW! See if bioClip matches gt2!
#**PROCESS:
#**  1. per subplot_id, list species that have prob > thresh across all grids. 
#**  2. Compare this list to gt2 list for each subplot_id

#*1. Function to list of  predicted species
nameAboveThreshPerGrid<-function(probVec, labs,t=0.05){
  return(labs[probVec>t])
}
nameAboveThreshPerSubplot<-function(predSubPlot, labs,t=0.05){
  probMatrix<-predSubPlot[,-metaCols]
  out<-apply(probMatrix,1,nameAboveThreshPerGrid,labs, t)
  return(unique(unlist(out)))
}
#check functions
#spId<-"31_1_1"
#nameAboveThreshPerGrid(probVec=pred[1,-metaCols], names(pred)[-metaCols], t=.1)
#nameAboveThreshPerSubplot(predSubPlot=pred[pred$subplot_id==spId,], labs=names(pred)[-metaCols], t=.1)

#*1. For all subplots, list predictd species
predList<-split(pred, f=pred$subplot_id)
predSpec<-lapply(predList,nameAboveThreshPerSubplot, names(pred)[-metaCols], t=.1)
#*2. For all subplot, cmpare predicted species to true species
sens<- unique(gt2[,1:3])
sens$predCnt<-NA
sens$tpCnt<-NA
sens$fnCnt<-NA
sens$fpCnt<-NA
sens$tnCnt<-NA

for (i in 1:dim(sens)[[1]]){
  trueSpecSubPl<-gt2[gt2$subplotID==sens$subplotID[i],]
  if (sens$subplotID[i] %in% names(predSpec)){
    predSpecSubPL<-predSpec[[which(names(predSpec)==sens$subplotID[i])]]
    sens$predCnt[i]<-length(predSpecSubPL)
    sens$tpCnt[i]<- length(which(trueSpecSubPl$trueSpecies%in% predSpecSubPL))
    sens$fnCnt[i]<- length(which(!(trueSpecSubPl$trueSpecies%in% predSpecSubPL)))
    sens$fpCnt[i]<- length(which(!(predSpecSubPL%in% trueSpecSubPl$trueSpecies)))
    sens$tnCnt[i]<- N-sens$trueCnt[i]-sens$fpCnt[i]
  }
}

sens$tpr<-sens$tpCnt/sens$trueCnt
sens$fnr<-sens$fnCnt/sens$trueCnt
sens$fpr<-sens$fpCnt/(N-sens$trueCnt)
sens$tnr<-sens$tnCnt/(N-sens$trueCnt)
sens$fdr<-sens$fpCnt/sens$predCnt #fals DISCOVERY rate
sens2<-na.omit(sens)
sens2
summary(sens2)
  
par(mfrow=c(2,2), mar=c(4,4,1,1))
plot(sens2$fpr, sens2$tpr, xlab="fpr", ylab="tpr", main=sens$plotID[1])
plot(sens2$fnr, sens2$tpr, xlab="fnr", ylab="tpr", main=sens$plotID[1])
plot(sens2$predCnt, sens2$trueCnt, xlab="# Predicted", ylab="True #", main=sens$plotID[1])
abline(a=0,b=1)


#*** create a tall, skinny data frame that includes true with predicted species and adds column

#*** ASSUMED FIRST subplotID in truth IS predicted by bioclip!!!

#Reformat predSpec
temp<- unique(gt2[,1:3])

#create data for first row in gt (this is why make above assumption)
predOnly<-cbind(plotID=temp$plotID[1], subplotID=temp$subplotID[1], 
                trueCnt=length(predSpec[[1]]), 
                trueSpecies=predSpec[[1]], source="pred")
trueOnly<-cbind(gt2[gt2$subplotID==gt2$subplotID[1],], source="true")


for (i in 2:dim(temp)[[1]]){
  if (temp$subplotID[i] %in% names(predSpec)){
    #pull from gt2
    trueOnly<-as.data.frame(rbind(trueOnly,
                                  cbind(gt2[gt2$subplotID==temp$subplotID[i],], source="true")))
    
    #predicted
    predSpecSubPL<-predSpec[[which(names(predSpec)==sens$subplotID[i])]]
  
    predOnly<-as.data.frame(rbind(predOnly,cbind(plotID=temp$plotID[i], 
                                                 subplotID=temp$subplotID[i], 
                                                 trueCnt=length(predSpecSubPL), 
                                                 trueSpecies=predSpecSubPL, source="pred")))
  }
}

predAndTrue<-rbind(trueOnly, predOnly)
#rename columns
names(predAndTrue)<-c("plotID","subplotID","cnt", "species", "source")
predAndTrue
