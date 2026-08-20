#############################################################################################
#' @title create result figures of initial bioclip results
#'
#' @author
#' Dave T Barnett \email{dbarnettl@battelleecology.org} \cr
#'
#' @description creates figures of initial results by pulling in results summary file 
#'
#' @import dplyr ggplot2
#'
#' @param 
#'
#' @return Script returns figures of results, 3 for the species list constraints, 3 for the grid size
#'
#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007
#'
#' @export
#'
#' @examples
#' 



# changelog and author contributions / copyrights
#   Dave T Barnett (2026-08-20)
#     original creation
##############################################################################################


rm(list = ls())

library(dplyr)
library(ggplot2)

results <- read.csv('C:/Users/dbarnett/Documents/GitHub/PlotDiversiVision/playWithBioClipResults/tprFpr.csv', stringsAsFactors = FALSE)


#1. Setup and structure check

table(results$site, results$labels)

results <- results %>%
  mutate(listSize = tpCnt + fnCnt + fpCnt + tnCnt)

results %>%
  group_by(site, labels) %>%
  summarise(nSubplot = n(),
            listMin = min(listSize),
            listMax = max(listSize),
            nTrueZero = sum(trueCnt == 0),
            nPredZero = sum(predCnt == 0),
            .groups = 'drop')

#1b. grid conversion
gridOrder <- sort(unique(results$grid))
gridOrder

results <- results %>%
  mutate(nTile = grid ^ 2)

tileOrder <- sort(unique(results$nTile))
tileOrder

results <- results %>%
  mutate(tileF = factor(nTile, levels = tileOrder))

table(results$site, results$tileF)

#2. Macro summary, mean of per-subplot rates
macro <- results %>%
  group_by(site, labels) %>%
  summarise(nSubplot = n(),
            tprMean = mean(tpr, na.rm = TRUE),
            tprSd = sd(tpr, na.rm = TRUE),
            tprMedian = median(tpr, na.rm = TRUE),
            fnrMean = mean(fnr, na.rm = TRUE),
            fprMean = mean(fpr, na.rm = TRUE),
            tnrMean = mean(tnr, na.rm = TRUE),
            fdrMean = mean(fdr, na.rm = TRUE),
            fdrSd = sd(fdr, na.rm = TRUE),
            fdrMedian = median(fdr, na.rm = TRUE),
            .groups = 'drop') %>%
  arrange(labels, site)

macro

#3. Micro summary, rates pooled from summed counts
micro <- results %>%
  group_by(site, labels) %>%
  summarise(tpSum = sum(tpCnt),
            fnSum = sum(fnCnt),
            fpSum = sum(fpCnt),
            tnSum = sum(tnCnt),
            .groups = 'drop') %>%
  mutate(tpr = tpSum / (tpSum + fnSum),
         fnr = fnSum / (tpSum + fnSum),
         fpr = fpSum / (fpSum + tnSum),
         tnr = tnSum / (fpSum + tnSum),
         fdr = fpSum / (tpSum + fpSum)) %>%
  arrange(labels, site)

micro %>%
  select(site, labels, tpr, fdr)


#4. Paired comparison across label constraints
setA <- results %>%
  filter(labels == 'bonap') %>%
  select(site, plotID, subplotID, tprA = tpr, fdrA = fdr)

setB <- results %>%
  filter(labels == 'OTHER') %>%
  select(plotID, subplotID, tprB = tpr, fdrB = fdr)

paired <- setA %>%
  left_join(setB, by = c('plotID', 'subplotID')) %>%
  mutate(deltaTpr = tprA - tprB,
         deltaFdr = fdrA - fdrB)

paired %>%
  group_by(site) %>%
  summarise(nSubplot = n(),
            meanDeltaTpr = mean(deltaTpr, na.rm = TRUE),
            meanDeltaFdr = mean(deltaFdr, na.rm = TRUE),
            .groups = 'drop')


labelOrder <- c('neon', 'gbif', 'bonap', 'state')

# confirm the vector matches the values actually present
setdiff(unique(results$labels), labelOrder)
setdiff(labelOrder, unique(results$labels))

results <- results %>%
  mutate(labels = factor(labels, levels = labelOrder))

levels(results$labels)

####species list constraint figures 
#true positive rate
ggplot(results, aes(x = labels, y = tpr, fill = labels)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.6) +
  facet_wrap(~ site) +
  scale_fill_brewer(palette = 'YlGnBu') +
  labs(x = 'possible species list constraint', y = 'true positive rate') +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1))

#false discovery rate
ggplot(results, aes(x = labels, y = fdr, fill = labels)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.6) +
  facet_wrap(~ site) +
  scale_fill_brewer(palette = 'YlGnBu') +
  labs(x = 'possible species list constraint', y = 'false discovery rate') +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1))

#richness obs vs predicted
ggplot(results, aes(x = trueCnt, y = predCnt)) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'grey40') +
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.4, size = 1, color = 'steelblue') +
  geom_smooth(method = 'lm', se = FALSE, color = 'firebrick', linewidth = 0.6) +
  facet_grid(labels ~ site) +
  labs(x = 'observed richness', y = 'predicted richness') +
  theme_bw()


##what showing

calib <- results %>%
  group_by(site, labels) %>%
  summarise(nSubplot = n(),
            meanTrue = mean(trueCnt),
            meanPred = mean(predCnt),
            bias = mean(predCnt - trueCnt),
            sdPred = sd(predCnt),
            corr = cor(trueCnt, predCnt),
            slope = coef(lm(predCnt ~ trueCnt))[2],
            intercept = coef(lm(predCnt ~ trueCnt))[1],
            .groups = 'drop')

calib


####grid numbers
#true positive rate
ggplot(results, aes(x = tileF, y = tpr, fill = tileF)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.6) +
  facet_wrap(~ site) +
  scale_fill_brewer(palette = 'OrRd') +
  labs(x = 'tiles per subplot image', y = 'true positive rate') +
  theme_bw() +
  theme(legend.position = 'none')

#false discovery rate
ggplot(results, aes(x = tileF, y = fdr, fill = tileF)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.6) +
  facet_wrap(~ site) +
  scale_fill_brewer(palette = 'OrRd') +
  labs(x = 'tiles per subplot image', y = 'false discovery rate') +
  theme_bw() +
  theme(legend.position = 'none')

#richness obs vs predicted
ggplot(results, aes(x = trueCnt, y = predCnt)) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'grey40') +
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.4, size = 1, color = 'chocolate3') +
  geom_smooth(method = 'lm', se = FALSE, color = 'darkslateblue', linewidth = 0.6) +
  facet_grid(tileF ~ site) +
  labs(x = 'observed richness', y = 'predicted richness') +
  theme_bw()

