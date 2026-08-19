setwd("~/Documents/research/imageomics/floraPalooza/playWithBioClipResults")
rm(list=ls())

#Text data packages 
#install.packages("tm", dependencies = TRUE)
library(tm)

#*** This code 
#*A) aligns all prediction filesnames, with label filenames, and ground truth filenames.
#*B) Then, calculates  tpr       fnr         fpr       tnr fdr for ALL sites, subplots, grids, labels-ised, etc.


#*************************
#**** A) Align file names
#*************************
#Line up which files should go together.  E.g., 
#"SCBI_005_neon_grid4.csv" + "SCBI_plot_labels.csv"+"SCBI_005_subplot_labels.csv"
#"SCBI_005_neon_grid2.csv" + "SCBI_plot_labels.csv"+"SCBI_005_subplot_labels.csv"

#**** Get FILE NAMES
#**filename for predictions: PlotDiversiVision/temp_results/...
allFilesPred<-system("ls ../PlotDiversiVision/temp_results",intern = TRUE)
#delete toy...
allFilesPred<-allFilesPred[-which(substr(allFilesPred,1,3)=="toy")]
tempList<-strsplit(allFilesPred, "_")
getItemInList<-function(list,i=1){
  return(list[i])
}
#make data frame
allFilesPred<-as.data.frame(allFilesPred)
names(allFilesPred)<-"predFile"
allFilesPred$site<-substr(allFilesPred$predFile,1,4)
allFilesPred$plotId<-substr(allFilesPred$predFile,1,8)
allFilesPred$labels<-unlist(lapply(tempList,getItemInList,3))
allFilesPred$grid<-substr(unlist(lapply(tempList,getItemInList,4)),5,5)


#delete rows with colorado  and plot labels
allFilesPred<-allFilesPred[allFilesPred$labels!="colorado",]
table(allFilesPred$site)
table(allFilesPred$plotId)
table(allFilesPred$labels)


#**filename for species labels: PlotDiversiVision/assets/species_list/...
allFilesPredLabels<-system("ls ../PlotDiversiVision/assets/species_list",intern = TRUE)
tempList2<-strsplit(allFilesPredLabels, "_")
allFilesPredLabels<-as.data.frame(allFilesPredLabels)
names(allFilesPredLabels)<-"predLabelsFile"
allFilesPredLabels$site<-unlist(lapply(tempList2,getItemInList,1))
allFilesPredLabels$labels<-unlist(lapply(tempList2,getItemInList,2))
#clean plot labels to align with other files
#delete  plot labels
allFilesPredLabels<-allFilesPredLabels[allFilesPredLabels$labels!="plot",]
allFilesPredLabels$labels[allFilesPredLabels$labels=="BONAPlist"]<-"bonap" 
allFilesPredLabels$labels[allFilesPredLabels$labels=="conus"]<-"state" 
allFilesPredLabels$labels<-tolower(allFilesPredLabels$labels)
table(allFilesPredLabels$site)
table(allFilesPredLabels$labels)

#**filename for ground truth: ../PlotDiversiVision/assets/test_labels/
allFilesTrue<-system("ls ../PlotDiversiVision/assets/test_labels",intern = TRUE)
allFilesTrue<-as.data.frame(allFilesTrue)
names(allFilesTrue)<-"truthFile"
allFilesTrue$site<-substr(allFilesTrue$truthFile,1,4)
allFilesTrue$plotId<-substr(allFilesTrue$truthFile,1,8)
dim(allFilesTrue)
table(allFilesTrue$site)
table(allFilesTrue$plotId)

#assess what I did:
allFilesPred[1:4,]
allFilesTrue[1:4,]
allFilesPredLabels[1:4,]


#*** MERGE Prediction data with ground truth data
predTrue<-merge(allFilesPred,allFilesTrue, by=c("site", "plotId"))
dim(allFilesPred)
dim(allFilesTrue)
dim(predTrue)
predTrue[1:4,]

#*** MERGE Prediction data + ground truth data WITH Label files 
#*BY SITE AND labels
table(allFilesPredLabels$labels)
table(predTrue$labels) 
dim(allFilesPredLabels)
dim(predTrue)
predTrueLabel<-merge(predTrue,allFilesPredLabels, by=c("site", "labels"), all=TRUE)
dim(predTrueLabel)
predTrueLabel[1:10,]

#* ADD INDEX Column for processing Bioclip
predTrueLabel$index<-1:dim(predTrueLabel)[[1]]



#*****************************
#** B) Process bioclip results
#*****************************

processBioClip<-function(rowIndex){
  #**predictions from BioClip
  filePred<-predTrueLabel$predFile[rowIndex]
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
  filePredLabels<-predTrueLabel$predLabelsFile[rowIndex]
  pathPredLabels<-paste("../PlotDiversiVision/assets/species_list/", filePredLabels,sep="")
  labelPred<-read.csv(pathPredLabels, sep=",", header=TRUE)
  dim(labelPred)
  names(labelPred)
  labs<-tolower(labelPred$resolved_labels)
  N<-length(labs) #number of species in the region
  
  #**Ground truth; i.e., species listed by NEON technician
  fileTrueLabels<-predTrueLabel$truthFile[rowIndex]
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
  sens2$site<-predTrueLabel$site[rowIndex]
  sens2$labels<-predTrueLabel$labels[rowIndex]
  sens2$grid<-predTrueLabel$grid[rowIndex]

  return(sens2)
}

allSens<-processBioClip(1)
for (k in 2:dim(predTrueLabel)[[1]]){
#for (k in 2:10){
    temp<-processBioClip(k)
  allSens<-rbind(allSens, temp)
}


write.table(allSens, "tprFpr.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)

