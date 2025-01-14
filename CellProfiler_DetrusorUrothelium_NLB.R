# Libraries
library(tidyverse)
library(smplot2)
library(lme4)
library(lmerTest)
library(Rmisc)

# Set working directory and plot directory
setwd("/Volumes/bedfordlab/Nocturia/Data_Analysis/Cell Profiler/Bladder_RNAscope")
getwd()
plotdir = paste(getwd(),"/Plots/", sep = "") 

### Import Detrusor Data ###
expanded1 <- read.csv("nuclei_expanded1_DET.csv")   
expanded2 <- read.csv("nuclei_expanded2_DET.csv")  
expanded3 <- read.csv("nuclei_expanded3_DET.csv")

# Extract the mouse ID and place it into the MouseID column to fix the mouseID 
expanded3new <- expanded3 %>%
  mutate(Metadata_Mouse = ifelse(is.na(Metadata_Mouse) | Metadata_Mouse == "", 
                                 sub("^[^_]+_[^_]+_[^_]+_[^_]+_[^_]+_([^_]+)_.*$", "\\1", FileName_Arntl), 
                                 MouseID))

### Import Urothelium Data ###
expanded1_uro <- read.csv("nuclei_expanded1_URO.csv") 
expanded2_uro <- read.csv("nuclei_expanded2_URO.csv")

# Extract the mouse ID and place it into the MouseID column to fix the mouseID 
expanded2_uroNew <- expanded2_uro %>%
  mutate(Metadata_Mouse = ifelse(is.na(Metadata_Mouse) | Metadata_Mouse == "", 
                                 sub("^[^_]+_[^_]+_[^_]+_[^_]+_[^_]+_([^_]+)_.*$", "\\1", FileName_DAPI), 
                                 MouseID))

# Combine data frames:
df_det = rbind(expanded1, expanded2, expanded3new)
df_uro = rbind(expanded1_uro, expanded2_uroNew)

df_det = df_det %>% dplyr::select(FileName_DAPI, Metadata_Mouse, Metadata_Person, Metadata_Time, Metadata_Tissue, Metadata_Slide, Metadata_Section,
                          Parent_nuclei, Children_Arntl_puncta_Count, Children_Per1_puncta_Count, Children_Piezo1_puncta_Count)

df_uro = df_uro %>% dplyr::select(FileName_DAPI, Metadata_Mouse, Metadata_Person, Metadata_Time, Metadata_Tissue, Metadata_Slide, Metadata_Section,
                          Parent_nuclei, Children_Arntl_puncta_Count, Children_Per1_puncta_Count, Children_Piezo1_puncta_Count)

df = rbind(df_det, df_uro)
df = df %>% mutate_if(is.character, as.factor)
unique(df$Metadata_Mouse)

# Consolidate MouseID
df = df %>%
  mutate(Metadata_Mouse = case_when(Metadata_Mouse == "81.00" ~ "81",
                                    Metadata_Mouse == "8101" ~ "81.01",
                                    Metadata_Mouse == "8103" ~ "81.03",
                                    Metadata_Mouse == "8602" ~ "86.02",
                                    Metadata_Mouse == "8704" ~ "87.04",
                                    Metadata_Mouse == "243304" ~ "2433.04",
                                    Metadata_Mouse == "ctrl6" ~ "ctrl.6",
                                    .default = Metadata_Mouse)) %>% droplevels()

# Make new Tissue column:
df = df %>% mutate(Tissue = case_when(grepl("trus", FileName_DAPI) ~ "Detrusor",
                          grepl("thel", FileName_DAPI) ~ "Urothelium",
                          .default = NA))

# Convert puncta counts to binary (pos or neg):
df <- df %>%
  mutate(
    Arntl_pos = case_when(
      Children_Arntl_puncta_Count == 0 ~ 0,
      Children_Arntl_puncta_Count > 0 ~ 1,
      .default = NA
    ),
    Per1_pos = case_when(
      Children_Per1_puncta_Count == 0 ~ 0,
      Children_Per1_puncta_Count > 0 ~ 1,
      .default = NA
    ),
    Piezo1_pos = case_when(
      Children_Piezo1_puncta_Count == 0 ~ 0,
      Children_Piezo1_puncta_Count > 0 ~ 1,
      .default = NA
    )
  )

# Summarise data
df_bin = df %>% group_by(FileName_DAPI, Metadata_Mouse, Metadata_Person, Metadata_Time, Tissue, Metadata_Slide, Metadata_Section) %>%
  dplyr::summarise(numDAPI = n(),
                   numPos_Arntl = sum(Arntl_pos),
                   numPos_Per1 = sum(Per1_pos),
                   numPos_Piezo1 = sum(Piezo1_pos)) %>%
  ungroup()

df_bin = df_bin %>% 
  mutate(pct_Arntl = numPos_Arntl / numDAPI * 100) %>%
  mutate(pct_Per1 = numPos_Per1 / numDAPI * 100) %>%
  mutate(pct_Piezo1 = numPos_Piezo1 / numDAPI * 100)

# Read in lookup table
lut <- read.csv("/Volumes/bedfordlab/Nocturia/Data_Analysis/Nocturia Lookup Tables - Histology_LUT.csv", stringsAsFactors = T)
str(lut)

# Merge df_bin and lut
df_bin = left_join(df_bin, lut)

## Make "Quantitative" data frame:
hist(df$Children_Arntl_puncta_Count)
hist(df$Children_Per1_puncta_Count)
hist(df$Children_Piezo1_puncta_Count)

df_quant_Arntl = df %>% filter(Children_Arntl_puncta_Count > 0) %>% droplevels()
df_quant_Per1 = df %>% filter(Children_Per1_puncta_Count > 0) %>% droplevels()
df_quant_Piezo1 = df %>% filter(Children_Piezo1_puncta_Count > 0) %>% droplevels()

df_quant_Arntl = df_quant_Arntl %>% group_by(FileName_DAPI) %>%
  dplyr::summarise(avg_Arntl_puncta_perCell = mean(Children_Arntl_puncta_Count, na.rm = T))

df_quant_Per1 = df_quant_Per1 %>% group_by(FileName_DAPI) %>%
  dplyr::summarise(avg_Per1_puncta_perCell = mean(Children_Per1_puncta_Count, na.rm = T))

df_quant_Piezo1 = df_quant_Piezo1 %>% group_by(FileName_DAPI) %>%
  dplyr::summarise(avg_Piezo1_puncta_perCell = mean(Children_Piezo1_puncta_Count, na.rm = T))

# Join databases:
df2 = df_bin %>% 
  full_join(df_quant_Arntl, by = "FileName_DAPI") %>%
  full_join(df_quant_Per1, by = "FileName_DAPI") %>%
  full_join(df_quant_Piezo1, by = "FileName_DAPI")

df2<- df2 %>%
  mutate(Age = fct_relevel(Age, c("Young", "Old")))

df2 = df2 %>% 
  unite(Group, c("Sex", "Tissue"), remove = FALSE) 

Group = unique(df2$Group)

for (i in Group){
  fo <- df2 %>%
    filter(Group == i)
  
  ggplot(data = fo %>% filter(!is.na(avg_Piezo1_puncta_perCell)),
         aes(x = Time, y = avg_Piezo1_puncta_perCell)) +
    facet_wrap(~Age) +
    sm_bar() +
    ylim(0, 6.2) +
    ggtitle(i)
  ggsave(paste0(plotdir, "avg_Piezo1_puncta_perCell_", i, ".pdf"), width=10, height=8, units = "cm", useDingbats = FALSE)
  
  ggplot(data = fo %>% filter(pct_Piezo1 > 0 & pct_Piezo1 < 100),
         aes(x = Time, y = pct_Piezo1)) +
    facet_wrap(~Age) +
    sm_bar() +
    ylim(0, 100) +
    ggtitle(i)
  ggsave(paste0(plotdir, "pct_Piezo1_", i, ".pdf"), width=10, height=8, units = "cm", useDingbats = FALSE)
  
}

# Is there a strong expression difference between urothelium and detrusor?
summary(lm(pct_Arntl ~ Tissue, data = df2)) # NS
summary(lm(avg_Arntl_puncta_perCell ~ Tissue, data = df2)) # NS
summary(lm(pct_Per1 ~ Tissue, data = df2)) # NS
summary(lm(avg_Per1_puncta_perCell ~ Tissue, data = df2)) # NS
summary(lm(pct_Piezo1 ~ Tissue, data = df2)) # 0.17
summary(lm(avg_Piezo1_puncta_perCell ~ Tissue, data = df2)) # 0.0413 *

# Combine Detrusor and Urothelium data (control for Tissue in Piezo1 models only)
Sex = unique(df2$Sex)

for (i in Sex){
  foo <- df2 %>%
    filter(Sex == i)
  
  ggplot(data = foo %>% filter(!is.na(avg_Piezo1_puncta_perCell)),
         aes(x = Time, y = avg_Piezo1_puncta_perCell, color = Tissue)) +
    facet_wrap(~Age) +
    sm_bar() +
    ylim(0, 6.2) +
    ggtitle(i)
  ggsave(paste0(plotdir, "avg_Piezo1_puncta_perCell_Tissue_", i, ".pdf"), width=10, height=8, units = "cm", useDingbats = FALSE)
  
  ggplot(data = foo %>% filter(!is.na(avg_Piezo1_puncta_perCell)),
         aes(x = Time, y = avg_Piezo1_puncta_perCell, color = Metadata_Mouse)) +
    facet_wrap(~Age) +
    sm_bar() +
    ylim(0, 6.2) +
    ggtitle(i)
  ggsave(paste0(plotdir, "avg_Piezo1_puncta_perCell_Mouse_", i, ".pdf"), width=10, height=8, units = "cm", useDingbats = FALSE)
  
}

# Sample size: 
ss = df2 %>% select(Age, Sex, Mouse_ID) %>% unique()
length(unique(ss$Mouse_ID))
table(ss$Sex, ss$Age)

## Stats/Plots for all gene markers:
ggplot(data=df2, aes(x=Time, y=pct_Arntl, fill = Time)) +
  facet_wrap(Sex ~ Age) +
  geom_boxplot() +
  ggtitle("pct_Arntl")
summary(lm(pct_Arntl ~ Time*Age, data = df2 %>% filter(Sex == "F"))) # 5.90e-08 ***
summary(lm(pct_Arntl ~ Time*Age, data = df2 %>% filter(Sex == "M"))) # 1.28e-05 ***

summary(lm(pct_Arntl ~ Time, data = df2 %>% filter(Sex == "F" & Age == "Young"))) # 4.27e-16 ***
summary(lm(pct_Arntl ~ Time, data = df2 %>% filter(Sex == "F" & Age == "Old"))) # 0.000166 ***
summary(lm(pct_Arntl ~ Time, data = df2 %>% filter(Sex == "M" & Age == "Young"))) # NS
summary(lm(pct_Arntl ~ Time, data = df2 %>% filter(Sex == "M" & Age == "Old"))) # 1.73e-07 ***

summary(lmer(pct_Arntl ~ Time*Age + (1|Metadata_Mouse), data = df2 %>% filter(Sex == "F"))) # 0.017885 *
summary(lmer(pct_Arntl ~ Time*Age + (1|Metadata_Mouse), data = df2 %>% filter(Sex == "M"))) # NS

summary(lmer(pct_Arntl ~ Time + (1|Metadata_Mouse), data = df2 %>% filter(Sex == "F" & Age == "Young"))) # 0.001834 **
summary(lmer(pct_Arntl ~ Time + (1|Metadata_Mouse), data = df2 %>% filter(Sex == "F" & Age == "Old"))) # NS
summary(lmer(pct_Arntl ~ Time + (1|Metadata_Mouse), data = df2 %>% filter(Sex == "M" & Age == "Young"))) # NS
summary(lmer(pct_Arntl ~ Time + (1|Metadata_Mouse), data = df2 %>% filter(Sex == "M" & Age == "Old"))) # NS

# Percent Arntl Effect size:
summarySE(data = df2, measurevar = "pct_Arntl", groupvars = c("Age", "Sex", "Time"), na.rm = T)

ggplot(data=df2, aes(x=Time, y=avg_Arntl_puncta_perCell, fill = Time)) +
  facet_wrap(Sex ~ Age) +
  geom_boxplot() +
  ggtitle("avg_Arntl_puncta_perCell")
summary(lm(avg_Arntl_puncta_perCell ~ Time*Age, data = df2 %>% filter(Sex == "F"))) # NS
summary(lm(avg_Arntl_puncta_perCell ~ Time*Age, data = df2 %>% filter(Sex == "M"))) # NS

ggplot(data=df2, aes(x=Time, y=pct_Per1, fill = Time)) +
  facet_wrap(Sex ~ Age) +
  geom_boxplot() +
  ggtitle("pct_Per1")
summary(lm(pct_Per1 ~ Time*Age, data = df2 %>% filter(Sex == "F"))) # 2.85e-08 ***
summary(lm(pct_Per1 ~ Time*Age, data = df2 %>% filter(Sex == "M"))) # 0.007962 **

summary(lm(pct_Per1 ~ Time, data = df2 %>% filter(Sex == "F" & Age == "Young"))) # 0.00457 **
summary(lm(pct_Per1 ~ Time, data = df2 %>% filter(Sex == "F" & Age == "Old"))) # 6.97e-07 ***
summary(lm(pct_Per1 ~ Time, data = df2 %>% filter(Sex == "M" & Age == "Young"))) # NS
summary(lm(pct_Per1 ~ Time, data = df2 %>% filter(Sex == "M" & Age == "Old"))) # 0.02390 *

summary(lmer(pct_Per1 ~ Time*Age + (1|Metadata_Mouse), data = df2 %>% filter(Sex == "F"))) # 0.0608 .
summary(lmer(pct_Per1 ~ Time*Age + (1|Metadata_Mouse), data = df2 %>% filter(Sex == "M"))) # NS

summary(lmer(pct_Per1 ~ Time + (1|Metadata_Mouse), data = df2 %>% filter(Sex == "F" & Age == "Young"))) # 0.1333
summary(lmer(pct_Per1 ~ Time + (1|Metadata_Mouse), data = df2 %>% filter(Sex == "F" & Age == "Old"))) # NS
summary(lmer(pct_Per1 ~ Time + (1|Metadata_Mouse), data = df2 %>% filter(Sex == "M" & Age == "Young"))) # NS
summary(lmer(pct_Per1 ~ Time + (1|Metadata_Mouse), data = df2 %>% filter(Sex == "M" & Age == "Old"))) # NS

# Percent Per1 Effect size:
summarySE(data = df2, measurevar = "pct_Per1", groupvars = c("Age", "Sex", "Time"), na.rm = T)

ggplot(data=df2, aes(x=Time, y=avg_Per1_puncta_perCell, fill = Time)) +
  facet_wrap(Sex ~ Age) +
  geom_boxplot() +
  ggtitle("avg_Per1_puncta_perCell")
summary(lm(avg_Per1_puncta_perCell ~ Time*Age, data = df2 %>% filter(Sex == "F"))) # NS
summary(lm(avg_Per1_puncta_perCell ~ Time*Age, data = df2 %>% filter(Sex == "M"))) # NS

ggplot(data=df2, aes(x=Time, y=pct_Piezo1, fill = Time)) +
  facet_wrap(Sex ~ Age) +
  geom_boxplot() +
  ggtitle("pct_Piezo1")
summary(lm(pct_Piezo1 ~ Time*Age, data = df2 %>% filter(Sex == "F"))) # 8.52e-07 ***
summary(lm(pct_Piezo1 ~ Time*Age, data = df2 %>% filter(Sex == "M"))) # NS

summary(lmer(pct_Piezo1 ~ Time + Tissue + (1|Metadata_Mouse), data = df2 %>% filter(Sex == "F"))) # NS
summary(lmer(pct_Piezo1 ~ Time + Tissue + (1|Metadata_Mouse), data = df2 %>% filter(Sex == "M"))) # NS

summary(lmer(pct_Piezo1 ~ Time + (1|Metadata_Mouse), data = df2 %>% filter(Sex == "F" & Age == "Young"))) # NS
summary(lmer(pct_Piezo1 ~ Time + (1|Metadata_Mouse), data = df2 %>% filter(Sex == "F" & Age == "Old"))) # NS
summary(lmer(pct_Piezo1 ~ Time + (1|Metadata_Mouse), data = df2 %>% filter(Sex == "M" & Age == "Young"))) # NS
summary(lmer(pct_Piezo1 ~ Time + (1|Metadata_Mouse), data = df2 %>% filter(Sex == "M" & Age == "Old"))) # NS

summary(lmer(pct_Piezo1 ~ Time + Tissue + Sex + Age + (1|Metadata_Mouse), data = df2))
summary(lmer(avg_Piezo1_puncta_perCell ~ Time + Tissue + Sex + Age + (1|Metadata_Mouse), data = df2))

ggplot(data=df2, aes(x=Time, y=avg_Piezo1_puncta_perCell, fill = Time)) +
  facet_wrap(Age ~ Sex) +
  geom_boxplot() +
  ggtitle("avg_Piezo1_puncta_perCell")
summary(lm(avg_Piezo1_puncta_perCell ~ Time*Age, data = df2 %>% filter(Sex == "F"))) # 0.034 *
summary(lm(avg_Piezo1_puncta_perCell ~ Time*Age, data = df2 %>% filter(Sex == "M"))) # 0.00104 **

summary(lmer(avg_Piezo1_puncta_perCell ~ Time*Age + Tissue + (1|Metadata_Mouse), data = df2 %>% filter(Sex == "F"))) # NS
summary(lmer(avg_Piezo1_puncta_perCell ~ Time*Age + Tissue + (1|Metadata_Mouse), data = df2 %>% filter(Sex == "M"))) # 0.030507 *

summary(lmer(avg_Piezo1_puncta_perCell ~ Time + Tissue + (1|Metadata_Mouse), data = df2 %>% filter(Sex == "F" & Age == "Young"))) # NS
summary(lmer(avg_Piezo1_puncta_perCell ~ Time + Tissue + (1|Metadata_Mouse), data = df2 %>% filter(Sex == "F" & Age == "Old"))) # NS
summary(lmer(avg_Piezo1_puncta_perCell ~ Time + Tissue + (1|Metadata_Mouse), data = df2 %>% filter(Sex == "M" & Age == "Young"))) # 0.000192 ***
summary(lmer(avg_Piezo1_puncta_perCell ~ Time + Tissue + (1|Metadata_Mouse), data = df2 %>% filter(Sex == "M" & Age == "Old"))) # NS

# Piezo1 puncta Effect size:
sum = summarySE(data = df2, measurevar = "avg_Piezo1_puncta_perCell", groupvars = c("Age", "Sex", "Time"), na.rm = T)
round((sum[4,5] / sum[3,5]), digits = 3) # adult male
round((sum[8,5] / sum[7,5]), digits = 3) # aged male

# Make GOOD plots:
library(Rmisc)

## ARNTL/BMAL Plots:

summaryF <- summarySE(data = df2 %>% filter(Sex == "F"),
                      measurevar = "pct_Arntl", groupvars = c("Age", "Time"), na.rm=T)
ggplot(data=df2 %>% filter(Sex == "F"), 
       aes(x=Time, y=pct_Arntl, fill = Time)) +
  facet_wrap(~Age) +
  geom_bar(data=summaryF, aes(fill=Time), stat="identity", width=0.85) +
  geom_errorbar(data=summaryF, aes(ymin=pct_Arntl-se, ymax=pct_Arntl+se, x=Time), size=0.5, width=0.25) + 
  geom_point(aes(fill=Time), colour="black", pch=21, position=position_jitter(width = 0.1), size=1) +
  scale_fill_manual(values=c("#f7941d", "#262262")) +
  coord_cartesian(ylim = c(0,85)) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig5.pct_Arntl_F.pdf", width=5.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)

summaryM <- summarySE(data = df2 %>% filter(Sex == "M"),
                      measurevar = "pct_Arntl", groupvars = c("Age", "Time"), na.rm=T)
ggplot(data=df2 %>% filter(Sex == "M"), 
       aes(x=Time, y=pct_Arntl, fill = Time)) +
  facet_wrap(~Age) +
  geom_bar(data=summaryM, aes(fill=Time), stat="identity", width=0.85) +
  geom_errorbar(data=summaryM, aes(ymin=pct_Arntl-se, ymax=pct_Arntl+se, x=Time), size=0.5, width=0.25) + 
  geom_point(aes(fill=Time), colour="black", pch=21, position=position_jitter(width = 0.1), size=1) +
  scale_fill_manual(values=c("#f7941d", "#262262")) +
  coord_cartesian(ylim = c(0,85)) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig5.pct_Arntl_M.pdf", width=5.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)

## PER1 Plots:

summaryF <- summarySE(data = df2 %>% filter(Sex == "F"),
                      measurevar = "pct_Per1", groupvars = c("Age", "Time"), na.rm=T)
ggplot(data=df2 %>% filter(Sex == "F"), 
       aes(x=Time, y=pct_Per1, fill = Time)) +
  facet_wrap(~Age) +
  geom_bar(data=summaryF, aes(fill=Time), stat="identity", width=0.85) +
  geom_errorbar(data=summaryF, aes(ymin=pct_Per1-se, ymax=pct_Per1+se, x=Time), size=0.5, width=0.25) + 
  geom_point(aes(fill=Time), colour="black", pch=21, position=position_jitter(width = 0.1), size=1) +
  scale_fill_manual(values=c("#f7941d", "#262262")) +
  coord_cartesian(ylim = c(0,86)) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig5.pct_Per1_F.pdf", width=5.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)

summaryM <- summarySE(data = df2 %>% filter(Sex == "M"),
                      measurevar = "pct_Per1", groupvars = c("Age", "Time"), na.rm=T)
ggplot(data=df2 %>% filter(Sex == "M"), 
       aes(x=Time, y=pct_Per1, fill = Time)) +
  facet_wrap(~Age) +
  geom_bar(data=summaryM, aes(fill=Time), stat="identity", width=0.85) +
  geom_errorbar(data=summaryM, aes(ymin=pct_Per1-se, ymax=pct_Per1+se, x=Time), size=0.5, width=0.25) + 
  geom_point(aes(fill=Time), colour="black", pch=21, position=position_jitter(width = 0.1), size=1) +
  scale_fill_manual(values=c("#f7941d", "#262262")) +
  coord_cartesian(ylim = c(0,86)) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig5.pct_Per1_M.pdf", width=5.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)

## PIEZO1 Plots:

summaryF <- summarySE(data = df2 %>% filter(Sex == "F"),
                      measurevar = "avg_Piezo1_puncta_perCell", groupvars = c("Age", "Time"), na.rm=T)
ggplot(data=df2 %>% filter(Sex == "F"), 
       aes(x=Time, y=avg_Piezo1_puncta_perCell, fill = Time)) +
  facet_wrap(~Age) +
  geom_bar(data=summaryF, aes(fill=Time), stat="identity", width=0.85) +
  geom_errorbar(data=summaryF, aes(ymin=avg_Piezo1_puncta_perCell-se, ymax=avg_Piezo1_puncta_perCell+se, x=Time), size=0.5, width=0.25) + 
  geom_point(aes(fill=Time), colour="black", pch=21, position=position_jitter(width = 0.1), size=1) +
  scale_fill_manual(values=c("#f7941d", "#262262")) +
  coord_cartesian(ylim = c(0,6.2)) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig5.avg_Piezo1_puncta_perCell_F.pdf", width=5.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)

summaryM <- summarySE(data = df2 %>% filter(Sex == "M"),
                      measurevar = "avg_Piezo1_puncta_perCell", groupvars = c("Age", "Time"), na.rm=T)
ggplot(data=df2 %>% filter(Sex == "M"), 
       aes(x=Time, y=avg_Piezo1_puncta_perCell, fill = Time)) +
  facet_wrap(~Age) +
  geom_bar(data=summaryM, aes(fill=Time), stat="identity", width=0.85) +
  geom_errorbar(data=summaryM, aes(ymin=avg_Piezo1_puncta_perCell-se, ymax=avg_Piezo1_puncta_perCell+se, x=Time), size=0.5, width=0.25) + 
  geom_point(aes(fill=Time), colour="black", pch=21, position=position_jitter(width = 0.1), size=1) +
  scale_fill_manual(values=c("#f7941d", "#262262")) +
  coord_cartesian(ylim = c(0,6.2)) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig5.avg_Piezo1_puncta_perCell_M.pdf", width=5.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)


