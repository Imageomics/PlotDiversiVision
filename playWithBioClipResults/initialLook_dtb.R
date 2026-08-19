
rm(list = ls())

library(plyr)
library(dplyr)
library(stringr)
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

#5. Distribution shape
ggplot(results, aes(x = labels, y = tpr, fill = labels)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.6) +
  facet_wrap(~ site) +
  labs(x = 'label constraint', y = 'true positive rate') +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(results, aes(x = labels, y = fdr, fill = labels)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.6) +
  facet_wrap(~ site) +
  labs(x = 'label constraint', y = 'false discovery rate') +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1))

#6. Richness calibration check
ggplot(results, aes(x = trueCnt, y = predCnt, color = labels)) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.5, size = 1) +
  facet_wrap(~ site) +
  labs(x = 'observed richness', y = 'predicted richness') +
  theme_bw()


