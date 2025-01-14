require(tidyverse)
require(lme4)
require(lmerTest)

# Set working directory and plot directory
setwd("/Volumes/bedfordlab/Nocturia/Data_Analysis")
plotdir = paste(getwd(), "/Plots/", sep = "") 

# Read in the qPCR data:
df<- read.csv("/Volumes/bedfordlab/Nocturia/Data_Analysis/Nocturia Lookup Tables - qPCR.csv")

# Make new grouping variables:
df<- df %>% mutate(AgeTime = paste(Age, Time),
                   TissueGene = paste(Tissue, Gene))

df$AgeTime<- factor(df$AgeTime, levels = c("Young AM", "Young PM",
                                           "Old AM", "Old PM"))

tissueGene<- unique(df$TissueGene)

# Plot fold change for each gene
for (i in tissueGene){
  fo <- df %>%
    filter(TissueGene == i)
  
  ggplot(data = fo, aes(x = AgeTime, y = log10(foldChange))) +
    facet_grid(.~ Sex) +
    geom_boxplot() +
    geom_point(colour="black", pch=21, position=position_jitter(0.1), size=1) +
    ggtitle(i) +
    theme_classic()
  ggsave(paste0(plotdir, "qPCR/", i, ".pdf"), width = 20, height = 10, units = "cm", useDingbats = FALSE)
  
}

# Run stats:
library(broom)

stats<- df %>% group_by(TissueGene) %>%
  do(broom::tidy(lm(log1p(foldChange) ~ Age*Time + Sex, .))) %>% 
  filter(!term == "(Intercept)") %>%
  mutate(P_round = round(p.value, digits = 4)) %>%
  ungroup()
write.csv(stats, "Nocturia_qPCR_stats_log1p.csv")

hist(log10(df$foldChange), breaks = 50)
hist(log1p(df$foldChange), breaks = 50)

# Sample size: 
ss = df %>% select(Age, Sex, Mouse_ID) %>% unique()
length(unique(ss$Mouse_ID))
table(ss$Sex, ss$Age)

## BLADDER MODELS:

# Bladder Bmal1:
summary(lm(log1p(foldChange) ~ Age*Time, data = df %>% filter(Sex == "F", TissueGene == "Bladder Bmal1"))) # NS
summary(lm(log1p(foldChange) ~ Age + Time, data = df %>% filter(Sex == "F", TissueGene == "Bladder Bmal1"))) # NS
summary(lm(log1p(foldChange) ~ Time, data = df %>% filter(Sex == "F", Age == "Young", TissueGene == "Bladder Bmal1"))) # NS
summary(lm(log1p(foldChange) ~ Time, data = df %>% filter(Sex == "F", Age == "Old", TissueGene == "Bladder Bmal1"))) # NS

summary(lm(log1p(foldChange) ~ Age*Time, data = df %>% filter(Sex == "M", TissueGene == "Bladder Bmal1"))) # NS
summary(lm(log1p(foldChange) ~ Age + Time, data = df %>% filter(Sex == "M", TissueGene == "Bladder Bmal1"))) # NS
summary(lm(log1p(foldChange) ~ Time, data = df %>% filter(Sex == "M", Age == "Young", TissueGene == "Bladder Bmal1"))) # 0.1410
summary(lm(log1p(foldChange) ~ Time, data = df %>% filter(Sex == "M", Age == "Old", TissueGene == "Bladder Bmal1"))) # NS

# Bladder Per2:
summary(lm(log1p(foldChange) ~ Age*Time, data = df %>% filter(Sex == "F", TissueGene == "Bladder Per2"))) # Age*Time: 0.00939 **
summary(lm(log1p(foldChange) ~ Time, data = df %>% filter(Sex == "F", Age == "Young", TissueGene == "Bladder Per2"))) # 0.00461 **
summary(lm(log1p(foldChange) ~ Time, data = df %>% filter(Sex == "F", Age == "Old", TissueGene == "Bladder Per2"))) # NS

summary(lm(log1p(foldChange) ~ Age*Time, data = df %>% filter(Sex == "M", TissueGene == "Bladder Per2"))) # NS
summary(lm(log1p(foldChange) ~ Age + Time, data = df %>% filter(Sex == "M", TissueGene == "Bladder Per2"))) # NS
summary(lm(log1p(foldChange) ~ Time, data = df %>% filter(Sex == "M", Age == "Young", TissueGene == "Bladder Per2"))) # 0.1180
summary(lm(log1p(foldChange) ~ Time, data = df %>% filter(Sex == "M", Age == "Old", TissueGene == "Bladder Per2"))) # NS

# Bladder Piezo:
summary(lm(log1p(foldChange) ~ Age*Time, data = df %>% filter(Sex == "F", TissueGene == "Bladder Piezo"))) # NS
summary(lm(log1p(foldChange) ~ Time, data = df %>% filter(Sex == "F", Age == "Young", TissueGene == "Bladder Piezo"))) # 0.0415 *
summary(lm(log1p(foldChange) ~ Time, data = df %>% filter(Sex == "F", Age == "Old", TissueGene == "Bladder Piezo"))) # NS

summary(lm(log1p(foldChange) ~ Age*Time, data = df %>% filter(Sex == "M", TissueGene == "Bladder Piezo"))) # Time: 0.0367 *; Age*Time: 0.0450 *
summary(lm(log1p(foldChange) ~ Age + Time, data = df %>% filter(Sex == "M", TissueGene == "Bladder Piezo"))) # Time: 0.0367 *; Age*Time: 0.0450 *
summary(lm(log1p(foldChange) ~ Time, data = df %>% filter(Sex == "M", Age == "Young", TissueGene == "Bladder Piezo"))) # NS
summary(lm(log1p(foldChange) ~ Time, data = df %>% filter(Sex == "M", Age == "Old", TissueGene == "Bladder Piezo"))) # NS

## KIDNEY MODELS:

# kidney Bmal1:
summary(lm(log1p(foldChange) ~ Age*Time, data = df %>% filter(Sex == "F", TissueGene == "kidney Bmal1"))) # NS
summary(lm(log1p(foldChange) ~ Time, data = df %>% filter(Sex == "F", Age == "Young", TissueGene == "kidney Bmal1"))) # 0.01979 *
summary(lm(log1p(foldChange) ~ Time, data = df %>% filter(Sex == "F", Age == "Old", TissueGene == "kidney Bmal1"))) # 0.09650 .

summary(lm(log1p(foldChange) ~ Age*Time, data = df %>% filter(Sex == "M", TissueGene == "kidney Bmal1"))) # NS
summary(lm(log1p(foldChange) ~ Time, data = df %>% filter(Sex == "M", Age == "Young", TissueGene == "kidney Bmal1"))) # 0.00196 **
summary(lm(log1p(foldChange) ~ Time, data = df %>% filter(Sex == "M", Age == "Old", TissueGene == "kidney Bmal1"))) # NS

# kidney Per2:
summary(lm(log1p(foldChange) ~ Age*Time, data = df %>% filter(Sex == "F", TissueGene == "kidney Per2"))) # 0.05591 .
summary(lm(log1p(foldChange) ~ Time, data = df %>% filter(Sex == "F", Age == "Young", TissueGene == "kidney Per2"))) # 0.00487 **
summary(lm(log1p(foldChange) ~ Time, data = df %>% filter(Sex == "F", Age == "Old", TissueGene == "kidney Per2"))) # NS

summary(lm(log1p(foldChange) ~ Age*Time, data = df %>% filter(Sex == "M", TissueGene == "kidney Per2"))) # NS
summary(lm(log1p(foldChange) ~ Time, data = df %>% filter(Sex == "M", Age == "Young", TissueGene == "kidney Per2"))) # 4.2e-05 ***
summary(lm(log1p(foldChange) ~ Time, data = df %>% filter(Sex == "M", Age == "Old", TissueGene == "kidney Per2"))) # NS

df<- df %>%
  mutate(Age = fct_relevel(Age, c("Young", "Old"))) %>%
  mutate(foldChange_unTrans = foldChange) %>%
  mutate(foldChange = log1p(foldChange)) # Express as log10 (fold change + 1)

# What is the distribution of fold changes? Are there outliers?
hist(df$foldChange_unTrans, breaks = 50)
hist(df$foldChange, breaks = 50)

## BLADDER PLOTS:
library(Rmisc)

## Effect size table for Piezo1 bladder:
ES <- summarySE(data = df %>% filter(TissueGene == "Bladder Piezo"),
                measurevar = "foldChange", groupvars = c("Sex", "Age", "Time"), na.rm=T)

# (dark - light) / light
round((ES[6,5] - ES[5,5]) / ES[5,5], digits = 3) # Young M 33.1%
round((ES[8,5] - ES[7,5]) / ES[7,5], digits = 3) # Old M -58.2%

# Percent Difference
cumulative.segment_max_wide_sum2$Diff_pct<- (cumulative.segment_max_wide_sum2$dark - cumulative.segment_max_wide_sum2$light) / cumulative.segment_max_wide_sum2$light
ES_totAct2<- summarySE(data=cumulative.segment_max_wide_sum2, measurevar = "Diff_pct", groupvars = c("Sex", "Age"))
round(ES_totAct2$Diff_pct, digits = 1)

## Bladder Bmal1:

df %>% filter(TissueGene == "Bladder Bmal1") %>% 
  dplyr::summarise(max = max(foldChange, na.rm = TRUE))

summaryF <- summarySE(data = df %>% filter(TissueGene == "Bladder Bmal1", Sex == "F"),
                      measurevar = "foldChange", groupvars = c("Age", "Time"), na.rm=T)
ggplot(data=df %>% filter(TissueGene == "Bladder Bmal1", Sex == "F"), 
       aes(x=Time, y=foldChange, fill = Time)) +
  facet_wrap(~Age) +
  geom_bar(data=summaryF, aes(fill=Time), stat="identity", width=0.85) +
  geom_errorbar(data=summaryF, aes(ymin=foldChange-se, ymax=foldChange+se, x=Time), size=0.5, width=0.25) + 
  geom_point(aes(fill=Time), colour="black", pch=21, position=position_jitter(width = 0.1), size=1) +
  scale_fill_manual(values=c("#f7941d", "#262262")) +
  coord_cartesian(ylim = c(0, 1.68)) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig4.Bladder Bmal1_foldChange_F.pdf", width=5.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)

summaryM <- summarySE(data = df %>% filter(TissueGene == "Bladder Bmal1", Sex == "M"),
                      measurevar = "foldChange", groupvars = c("Age", "Time"), na.rm=T)
ggplot(data=df %>% filter(TissueGene == "Bladder Bmal1", Sex == "M"), 
       aes(x=Time, y=foldChange, fill = Time)) +
  facet_wrap(~Age) +
  geom_bar(data=summaryM, aes(fill=Time), stat="identity", width=0.85) +
  geom_errorbar(data=summaryM, aes(ymin=foldChange-se, ymax=foldChange+se, x=Time), size=0.5, width=0.25) + 
  geom_point(aes(fill=Time), colour="black", pch=21, position=position_jitter(width = 0.1), size=1) +
  scale_fill_manual(values=c("#f7941d", "#262262")) +
  coord_cartesian(ylim = c(0, 1.68)) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig4.Bladder Bmal1_foldChange_M.pdf", width=5.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)

## Bladder Per2:

df %>% filter(TissueGene == "Bladder Per2") %>% 
  dplyr::summarise(max = max(foldChange, na.rm = TRUE))

summaryF <- summarySE(data = df %>% filter(TissueGene == "Bladder Per2", Sex == "F"),
                      measurevar = "foldChange", groupvars = c("Age", "Time"), na.rm=T)
ggplot(data=df %>% filter(TissueGene == "Bladder Per2", Sex == "F"), 
       aes(x=Time, y=foldChange, fill = Time)) +
  facet_wrap(~Age) +
  geom_bar(data=summaryF, aes(fill=Time), stat="identity", width=0.85) +
  geom_errorbar(data=summaryF, aes(ymin=foldChange-se, ymax=foldChange+se, x=Time), size=0.5, width=0.25) + 
  geom_point(aes(fill=Time), colour="black", pch=21, position=position_jitter(width = 0.1), size=1) +
  scale_fill_manual(values=c("#f7941d", "#262262")) +
  coord_cartesian(ylim = c(0, 2.49)) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig4.Bladder Per2_foldChange_F.pdf", width=5.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)

summaryM <- summarySE(data = df %>% filter(TissueGene == "Bladder Per2", Sex == "M"),
                      measurevar = "foldChange", groupvars = c("Age", "Time"), na.rm=T)
ggplot(data=df %>% filter(TissueGene == "Bladder Per2", Sex == "M"), 
       aes(x=Time, y=foldChange, fill = Time)) +
  facet_wrap(~Age) +
  geom_bar(data=summaryM, aes(fill=Time), stat="identity", width=0.85) +
  geom_errorbar(data=summaryM, aes(ymin=foldChange-se, ymax=foldChange+se, x=Time), size=0.5, width=0.25) + 
  geom_point(aes(fill=Time), colour="black", pch=21, position=position_jitter(width = 0.1), size=1) +
  scale_fill_manual(values=c("#f7941d", "#262262")) +
  coord_cartesian(ylim = c(0, 2.49)) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig4.Bladder Per2_foldChange_M.pdf", width=5.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)

## Bladder Piezo:

df %>% filter(TissueGene == "Bladder Piezo") %>% 
  dplyr::summarise(max = max(foldChange, na.rm = TRUE))

summaryF <- summarySE(data = df %>% filter(TissueGene == "Bladder Piezo", Sex == "F"),
                      measurevar = "foldChange", groupvars = c("Age", "Time"), na.rm=T)
ggplot(data=df %>% filter(TissueGene == "Bladder Piezo", Sex == "F"), 
       aes(x=Time, y=foldChange, fill = Time)) +
  facet_wrap(~Age) +
  geom_bar(data=summaryF, aes(fill=Time), stat="identity", width=0.85) +
  geom_errorbar(data=summaryF, aes(ymin=foldChange-se, ymax=foldChange+se, x=Time), size=0.5, width=0.25) + 
  geom_point(aes(fill=Time), colour="black", pch=21, position=position_jitter(width = 0.1), size=1) +
  scale_fill_manual(values=c("#f7941d", "#262262")) +
  coord_cartesian(ylim = c(0, 2.66)) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig4.Bladder Piezo_foldChange_F.pdf", width=5.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)

summaryM <- summarySE(data = df %>% filter(TissueGene == "Bladder Piezo", Sex == "M"),
                      measurevar = "foldChange", groupvars = c("Age", "Time"), na.rm=T)
ggplot(data=df %>% filter(TissueGene == "Bladder Piezo", Sex == "M"), 
       aes(x=Time, y=foldChange, fill = Time)) +
  facet_wrap(~Age) +
  geom_bar(data=summaryM, aes(fill=Time), stat="identity", width=0.85) +
  geom_errorbar(data=summaryM, aes(ymin=foldChange-se, ymax=foldChange+se, x=Time), size=0.5, width=0.25) + 
  geom_point(aes(fill=Time), colour="black", pch=21, position=position_jitter(width = 0.1), size=1) +
  scale_fill_manual(values=c("#f7941d", "#262262")) +
  coord_cartesian(ylim = c(0, 2.66)) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig4.Bladder Piezo_foldChange_M.pdf", width=5.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)

## KIDNEY PLOTS:
## kidney Bmal1:

df %>% filter(TissueGene == "kidney Bmal1") %>% 
  dplyr::summarise(max = max(foldChange, na.rm = TRUE))

summaryF <- summarySE(data = df %>% filter(TissueGene == "kidney Bmal1", Sex == "F"),
                      measurevar = "foldChange", groupvars = c("Age", "Time"), na.rm=T)
ggplot(data=df %>% filter(TissueGene == "kidney Bmal1", Sex == "F"), 
       aes(x=Time, y=foldChange, fill = Time)) +
  facet_wrap(~Age) +
  geom_bar(data=summaryF, aes(fill=Time), stat="identity", width=0.85) +
  geom_errorbar(data=summaryF, aes(ymin=foldChange-se, ymax=foldChange+se, x=Time), size=0.5, width=0.25) + 
  geom_point(aes(fill=Time), colour="black", pch=21, position=position_jitter(width = 0.1), size=1) +
  scale_fill_manual(values=c("#f7941d", "#262262")) +
  coord_cartesian(ylim = c(0, 1.47)) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig4.kidney Bmal1_foldChange_F.pdf", width=5.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)

summaryM <- summarySE(data = df %>% filter(TissueGene == "kidney Bmal1", Sex == "M"),
                      measurevar = "foldChange", groupvars = c("Age", "Time"), na.rm=T)
ggplot(data=df %>% filter(TissueGene == "kidney Bmal1", Sex == "M"), 
       aes(x=Time, y=foldChange, fill = Time)) +
  facet_wrap(~Age) +
  geom_bar(data=summaryM, aes(fill=Time), stat="identity", width=0.85) +
  geom_errorbar(data=summaryM, aes(ymin=foldChange-se, ymax=foldChange+se, x=Time), size=0.5, width=0.25) + 
  geom_point(aes(fill=Time), colour="black", pch=21, position=position_jitter(width = 0.1), size=1) +
  scale_fill_manual(values=c("#f7941d", "#262262")) +
  coord_cartesian(ylim = c(0, 1.47)) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig4.kidney Bmal1_foldChange_M.pdf", width=5.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)

## kidney Per2:

df %>% filter(TissueGene == "kidney Per2") %>% 
  dplyr::summarise(max = max(foldChange, na.rm = TRUE))

summaryF <- summarySE(data = df %>% filter(TissueGene == "kidney Per2", Sex == "F"),
                      measurevar = "foldChange", groupvars = c("Age", "Time"), na.rm=T)
ggplot(data=df %>% filter(TissueGene == "kidney Per2", Sex == "F"), 
       aes(x=Time, y=foldChange, fill = Time)) +
  facet_wrap(~Age) +
  geom_bar(data=summaryF, aes(fill=Time), stat="identity", width=0.85) +
  geom_errorbar(data=summaryF, aes(ymin=foldChange-se, ymax=foldChange+se, x=Time), size=0.5, width=0.25) + 
  geom_point(aes(fill=Time), colour="black", pch=21, position=position_jitter(width = 0.1), size=1) +
  scale_fill_manual(values=c("#f7941d", "#262262")) +
  coord_cartesian(ylim = c(0, 2.41)) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig4.kidney Per2_foldChange_F.pdf", width=5.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)

summaryM <- summarySE(data = df %>% filter(TissueGene == "kidney Per2", Sex == "M"),
                      measurevar = "foldChange", groupvars = c("Age", "Time"), na.rm=T)
ggplot(data=df %>% filter(TissueGene == "kidney Per2", Sex == "M"), 
       aes(x=Time, y=foldChange, fill = Time)) +
  facet_wrap(~Age) +
  geom_bar(data=summaryM, aes(fill=Time), stat="identity", width=0.85) +
  geom_errorbar(data=summaryM, aes(ymin=foldChange-se, ymax=foldChange+se, x=Time), size=0.5, width=0.25) + 
  geom_point(aes(fill=Time), colour="black", pch=21, position=position_jitter(width = 0.1), size=1) +
  scale_fill_manual(values=c("#f7941d", "#262262")) +
  coord_cartesian(ylim = c(0, 2.41)) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig4.kidney Per2_foldChange_M.pdf", width=5.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)

