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

nodelabels(round(branching.times(tree), 1), cex = 1)
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

tree<- read.tree("C:/Users/Braec/Desktop/Florapalooza/PlotDiversiVision/assets/PlotVision_tree.tre")

pdist<-cophenetic(tree)

#############Begin recovered script

## ============================================================
## RECOVERED SCRIPT
## Reconstructed from R console history after accidental
## git-pull overwrite. Dead-end attempts, typos, and repeated
## re-definitions from the trace have been removed/consolidated.
## Please review before re-running -- a few ambiguous spots are
## flagged with comments below.
## ============================================================

library(ape)
library(ggplot2)
library(tidyr)

## ------------------------------------------------------------
## 1. Load tree & basic diagnostics
## ------------------------------------------------------------
tree <- read.tree("C:/Users/Braec/Desktop/Florapalooza/PlotDiversiVision/assets/PlotVision_tree.tre")

distance_matrix <- cophenetic(tree)
## ------------------------------------------------------------
## 2. Plot the time-dated phylogeny
## ------------------------------------------------------------
plot.phylo(
  tree,
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

# NOTE: original trace referenced `result_tree$scenario.3`, which was
# never defined in this session (threw "object 'result_tree' not found").
# The working fallback used in the trace was simply:
nodelabels(round(branching.times(tree), 1), cex = 1)

## ------------------------------------------------------------
## 3. Build long-format cophenetic distance data for heatmap
## ------------------------------------------------------------
Pdist_mat <- cophenetic(tree)

Pdist_df <- as.data.frame(Pdist_mat)
Pdist_df$sp1 <- gsub("_", " ", rownames(Pdist_mat))

Pdist_long <- tidyr::pivot_longer(
  Pdist_df,
  cols      = -sp1,
  names_to  = "sp2",
  values_to = "dist"
)
Pdist_long$sp2 <- gsub("_", " ", Pdist_long$sp2)

# Order species by hierarchical clustering for a nicer heatmap layout
sp_clust  <- hclust(as.dist(Pdist_mat))$order
sp_levels <- gsub("_", " ", rownames(Pdist_mat)[sp_clust])
Pdist_long$sp1 <- factor(Pdist_long$sp1, levels = sp_levels)
Pdist_long$sp2 <- factor(Pdist_long$sp2, levels = sp_levels)

## ------------------------------------------------------------
## 4. Custom theme
## ------------------------------------------------------------
theme_adm <- function() {
  theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 12),
      plot.subtitle    = element_text(size = 9, color = "grey40"),
      axis.title       = element_text(size = 10),
      axis.text        = element_text(size = 9),
      legend.title     = element_text(size = 9),
      legend.text      = element_text(size = 8),
      panel.grid.minor = element_blank()
    )
}

## ------------------------------------------------------------
## 5. Heatmap of pairwise phylogenetic (cophenetic) distance
## ------------------------------------------------------------

p2 <- ggplot(Pdist_long, aes(x = sp1, y = sp2, fill = dist)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("white", "#FDE68A", "#F97316", "#7C2D12"),
                       name = "dist") +
  labs(title    = "Phylogenetic Distance",
       subtitle = "Based on UPhyloMaker2",
       x = NULL, y = NULL) +
  theme_adm() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 7),
        axis.text.y = element_text(face = "italic", size = 7),
        panel.grid  = element_blank())

print(p2)

## ------------------------------------------------------------
## 6. Quick exploratory checks
## ------------------------------------------------------------
hist(Pdist_long$dist)
plot(Pdist_long$sp1, Pdist_long$dist)

# tprfpr <- read.csv("C:/Users/Braec/Desktop/Florapalooza/PlotDiversiVision/playWithBioClipResults/tprFpr.csv")
# View(tprfpr)
