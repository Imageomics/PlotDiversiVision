#begin Phylogenetic taxonomic data __ 
#Author Braedon Lineman 8/19/2026
######## REQUIRED************* Run "assessBioCLipOutcomes.r" first

Neon_plot_data<-read.csv("C:/Users/Braec/Desktop/Florapalooza/PlotDiversiVision/assets/NEON_plotData.csv")

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
names(predAndTrue)<-c("plotID","subplotID","n", "species","source")
predAndTrue

PhyloDF<- data.frame(scinam = Neon_plot_data$scientificName, family = Neon_plot_data$family)

print(predAndTrue)
print(PhyloDF)

#Phylodf creations
library(dplyr)
library(stringr)
library(tidyr)

library(V.PhyloMaker2)
library(ape) 
library(phytools) 
library(dplyr)
library(stringr)

# Extract just "genus species" from scinam and standardize case
PhyloDF_clean <- PhyloDF %>%
  mutate(
    genus_species = str_extract(scinam, "^\\S+\\s+\\S+"),      # first two words
    genus_species = str_to_lower(genus_species)
  ) %>%
  select(genus_species, family) %>%
  distinct(genus_species, .keep_all = TRUE)   # in case of duplicate name->family rows

# Standardize predAndTrue's species column the same way 
predAndTrue <- predAndTrue %>%
  mutate(species = str_to_lower(str_trim(species)))

# Join family onto predAndTrue
final_phylo_df <- predAndTrue %>%
  left_join(PhyloDF_clean, by = c("species" = "genus_species"))

# Check for anything that didn't match
unmatched <- final_phylo_df %>% filter(is.na(family)) %>% distinct(species)
print(unmatched)

# Look up each unmatched species' family 
# and fill in the correct family name below.
manual_family_lookup <- tibble::tribble(
  ~species,                    ~family,
  "persicaria perfoliata",     "Polygonaceae",
  "potentilla indica",         "Rosaceae",
  "persicaria longiseta",      "Polygonaceae",
  "persicaria virginiana",     "Polygonaceae",
  "cryptotaenia canadensis",   "Apiaceae"
)

final_phylo_df <- final_phylo_df %>%
  left_join(manual_family_lookup, by = "species", suffix = c("", "_manual")) %>%
  mutate(family = coalesce(family, family_manual)) %>%
  select(-family_manual)

# Confirm nothing is missing now
still_missing <- final_phylo_df %>% filter(is.na(family)) %>% distinct(species)
print(still_missing)

# split "genus species" into separate genus and species columns ---
final_phylo_df <- final_phylo_df %>%
  separate(species, into = c("genus", "species"), sep = " ", extra = "merge", remove = TRUE)

#Phylo analysis df for tree
TreeDF<- data_frame(species = final_phylo_df$species, genus = final_phylo_df$genus,family = final_phylo_df$family)

#phylo analysis
result_tree <- phylo.maker(TreeDF, tree     = GBOTB.extended.TPL,
                           nodes    = nodes.info.1.LCVP,
                           scenarios=c("S1","S2","S3"))


write.tree(result_tree$scenario.3, file = "C:/Users/Braec/Desktop/Florapalooza/PlotDiversiVision/assets/PlotVision_tree.tre")


### plot the phylogenies with node ages displayed.
par(mfrow = c(1, 3))
plot.phylo(result_tree$scenario.1, cex = 1.5, main = "scenario.1")

nodelabels(round(branching.times(result_tree$scenario.1), 1), cex = 1)

plot.phylo(result_tree$scenario.2, cex = 1.5, main = "scenario.2")

nodelabels(round(branching.times(result_tree$scenario.2), 1), cex = 1)

plot.phylo(result_tree$scenario.3, cex = 1.5, main = "scenario.3")

nodelabels(round(branching.times(result_tree$scenario.3), 1), cex = 1)
#plot

plot.phylo(
  result_tree$scenario.3,
  type           = "phylogram",
  show.tip.label = TRUE,
  tip.color      = "black",
  cex            = 0.80,
  no.margin      = TRUE,
  edge.width     = 1.5
)

# Geological time axis -- branch lengths are in millions of years (Ma)
axisPhylo(side = 1, cex.axis = 0.8)
mtext("Time (Ma before present)", side = 1, line = 2.5, cex = 0.9)
title("Time-Dated Phylogenetic Tree\n(V.PhyloMaker2, S3)")

