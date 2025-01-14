# Libraries
library(tidyverse)
library(lme4)
library(lmerTest)

# Set working directory and plot directory
setwd("/Volumes/bedfordlab/Nocturia/Data_Analysis")
plotdir = paste(getwd(), "/Plots/", sep = "") 

## TEMP LOGGER DATA

# Read in Tb table
Tb = read.csv("Nocturia Lookup Tables - Tb.csv", header = T, sep = ",")
Tb$Age = factor(Tb$Age, levels = c("Young", "Old"))
range(Tb$Amplitude)
TbF = Tb %>% filter(Sex == "F")
TbM = Tb %>% filter(Sex == "M")

# Make Plots:
library(Rmisc)

summaryF <- summarySE(TbF, measurevar="Amplitude", groupvars="Age", na.rm=T)
ggplot(data=TbF, aes(x=Age, y=Amplitude)) +
  geom_bar(data=summaryF, aes(fill=Age), stat="identity", width=0.85) +
  geom_errorbar(data=summaryF, aes(ymin=Amplitude-se, ymax=Amplitude+se, x=Age), size=0.5, width=0.25) + 
  geom_point(aes(fill=Age), colour="black", pch=21, position=position_jitter(0.1), size=1) +
  scale_fill_manual(values=c("#8a1e8e", "#642c66")) +
  xlab(NULL) + ylab(NULL) + 
  theme_classic() +
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank())
ggsave("Plots/Fig1.ClockLab_Amp_F.pdf", width=3, height=6, units="cm", dpi=1500, useDingbats=FALSE)

summaryM <- summarySE(TbM, measurevar="Amplitude", groupvars="Age", na.rm=T)
ggplot(data=TbM, aes(x=Age, y=Amplitude)) +
  geom_bar(data=summaryM, aes(fill=Age), stat="identity", width=0.85) +
  geom_errorbar(data=summaryM, aes(ymin=Amplitude-se, ymax=Amplitude+se, x=Age), size=0.5, width=0.25) + 
  geom_point(aes(fill=Age), colour="black", pch=21, position=position_jitter(0.1), size=1) +
  scale_fill_manual(values=c("#6c8c24", "#4d5c2c")) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() +
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank())
ggsave("Plots/Fig1.ClockLab_Amp_M.pdf", width=3, height=6, units="cm", dpi=1500, useDingbats=FALSE)

## Amplitude Stats:
summary(lm(Amplitude ~ Age*Sex, data = Tb)) # age: 0.01121 *
summary(lm(Amplitude ~ Age, data = TbF)) # 0.000759 ***
summary(lm(Amplitude ~ Age, data = TbM)) # 0.135

## Period Stats:
summary(lm(Period ~ Age*Sex, data = Tb)) # age: 0.700
summary(lm(Period ~ Age, data = TbF)) # 0.0185 *
summary(lm(Period ~ Age, data = TbM)) # 0.545

## WATER CONSUMPTION DATA

# Read in Lookup Table data (need weight info)
lut = read.csv("Nocturia Lookup Tables - LUT.csv", header = T, sep = ",")

weight = lut %>% group_by(Cage, Age, Sex) %>%
  dplyr::summarise(totWeight = sum(Weight),
            avgWeight = mean(Weight),
            numMice = n())

write.csv(weight, "weight.csv")

# Read in Water table
wc = read.csv("Nocturia Lookup Tables - Water.csv", header = T, sep = ",")
wc = left_join(wc, weight)
wc$Age = factor(wc$Age, levels = c("Young", "Old"))

# Standardize by mouse weight
wc = wc %>% mutate(Intake_Std_tot = Intake_24h / totWeight,
                   Intake_Std_avg = Intake_24h / avgWeight)

hist(wc$Intake_24h)
hist(wc$Intake_Std_tot)
hist(wc$Intake_Std_avg)

# Sample size
ss = wc %>% group_by(Cage, Age, Sex) %>% dplyr::summarise(n=n())
mean(ss$n)

## Stats:
summary(lmer(Intake_24h ~ Age*Sex + (1|Cage), data = wc)) # NS
summary(lmer(Intake_Std_tot ~ Age*Sex + (1|Cage), data = wc)) # NS
summary(lmer(Intake_Std_avg ~ Age*Sex + (1|Cage), data = wc)) # NS

## Variance Statistics:
var = summarySE(wc, measurevar="Intake_Std_tot", groupvars=c("Cage", "Age", "Sex"), na.rm=T)
summary(lm(sd ~ N, data = var)) # 0.0564 .
summary(lm(sd ~ Sex + Age, data = var)) # NS

# Make Plots:

summary = summarySE(wc, measurevar="Intake_Std_tot", groupvars=c("Age", "Sex"), na.rm=T)
ggplot(data = wc, aes(x=Age, y=Intake_Std_tot)) +
  facet_grid(.~ Sex) +
  geom_bar(data=summary, aes(fill=Age), stat="identity", width=0.85) +
  geom_point(aes(fill=Age), colour="black", pch=21, position=position_jitter(0.1), size=1) +
  geom_errorbar(data=summary, aes(ymin=Intake_Std_tot-se, ymax=Intake_Std_tot+se, x=Age), size=0.5, width=0.25) +
  ylab("24h H2O intake (mL/g mouse)") +
  theme_classic() + theme(legend.position = "none")
#ggsave("Plots/H2O_intake_BothSexes_allData.pdf", width=10, height=10, units="cm", dpi=1500, useDingbats=FALSE)

range(wc$Intake_Std_tot)

summaryF = summarySE(wc %>% filter(Sex == "F"), measurevar="Intake_Std_tot", groupvars="Age", na.rm=T)
ggplot(data=wc %>% filter(Sex == "F"), aes(x=Age, y=Intake_Std_tot)) +
  geom_bar(data=summaryF, aes(fill=Age), stat="identity", width=0.85) +
  geom_point(aes(fill=Age), colour="black", pch=21, position=position_jitter(0.1), size=1) +
  geom_errorbar(data=summaryF, aes(ymin=Intake_Std_tot-se, ymax=Intake_Std_tot+se, x=Age), size=0.5, width=0.25) + 
  scale_fill_manual(values=c("#8a1e8e", "#642c66")) +
  ylim(0, 0.47) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), panel.spacing = unit(0, "lines"),  
                          strip.background = element_blank(), strip.text.x = element_blank())
#ggsave("Plots/H2O_intake_F_allData.pdf", width=3, height=6, units="cm", dpi=1500, useDingbats=FALSE)

summaryM = summarySE(wc %>% filter(Sex == "M"), measurevar="Intake_Std_tot", groupvars="Age", na.rm=T)
ggplot(data=wc %>% filter(Sex == "M"), aes(x=Age, y=Intake_Std_tot)) +
  geom_bar(data=summaryM, aes(fill=Age), stat="identity", width=0.85) +
  geom_point(aes(fill=Age), colour="black", pch=21, position=position_jitter(0.1), size=1) +
  geom_errorbar(data=summaryM, aes(ymin=Intake_Std_tot-se, ymax=Intake_Std_tot+se, x=Age), size=0.5, width=0.25) + 
  scale_fill_manual(values=c("#6c8c24", "#4d5c2c")) +
  ylim(0, 0.47) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), panel.spacing = unit(0, "lines"),  
                          strip.background = element_blank(), strip.text.x = element_blank())
#ggsave("Plots/H2O_intake_M_allData.pdf", width=3, height=6, units="cm", dpi=1500, useDingbats=FALSE)

# Take average per cage
wc_avg = wc %>% group_by(Cage, Age, Sex) %>%
  dplyr::summarise(Avg_Intake = mean(Intake_24h), n = n())

wcF = wc_avg %>% filter(Sex == "F")
wcM = wc_avg %>% filter(Sex == "M")

summaryF <- summarySE(wcF, measurevar="Avg_Intake", groupvars="Age", na.rm=T)
ggplot(data=wcF, aes(x=Age, y=Avg_Intake)) +
  geom_bar(data=summaryF, aes(fill=Age), stat="identity", width=0.85) +
  geom_errorbar(data=summaryF, aes(ymin=Avg_Intake-se, ymax=Avg_Intake+se, x=Age), size=0.5, width=0.25) + 
  geom_point(aes(fill=Age), colour="black", pch=21, position=position_jitter(0.1), size=1) +
  scale_fill_manual(values=c("#8a1e8e", "#642c66")) +
  ylim(0, 30) +
  xlab(NULL) + ylab(NULL) + 
  theme_classic() +
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank())
ggsave("Plots/Fig2.H2O_intake_F.pdf", width=3, height=6, units="cm", dpi=1500, useDingbats=FALSE)

summaryM <- summarySE(wcM, measurevar="Avg_Intake", groupvars="Age", na.rm=T)
ggplot(data=wcM, aes(x=Age, y=Avg_Intake)) +
  geom_bar(data=summaryM, aes(fill=Age), stat="identity", width=0.85) +
  geom_errorbar(data=summaryM, aes(ymin=Avg_Intake-se, ymax=Avg_Intake+se, x=Age), size=0.5, width=0.25) + 
  geom_point(aes(fill=Age), colour="black", pch=21, position=position_jitter(0.1), size=1) +
  scale_fill_manual(values=c("#6c8c24", "#4d5c2c")) +
  ylim(0, 30) +
  xlab(NULL) + ylab(NULL) + 
  theme_classic() +
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank())
ggsave("Plots/Fig2.H2O_intake_M.pdf", width=3, height=6, units="cm", dpi=1500, useDingbats=FALSE)

# Join water intake data with filter paper data: (From Nocturia_FP_v4_NLB.R)
hc_new = read.csv("HCLC_filterPapers.csv", header = T, sep = ",")
hc_new = hc_new %>% dplyr::mutate(across(where(is_character), as_factor))
hc_new = hc_new %>% dplyr::mutate(across(where(is_integer), as_factor))
hc_new = hc_new[,-1]
str(hc_new)

wc_long = wc %>% pivot_longer(
  cols = starts_with("FP_"),
  names_to = "Time_cageType",
  names_prefix = "FP_",
  values_to = "file_name")

wc_long = wc_long %>% drop_na(file_name)
wc_long = wc_long %>% dplyr::mutate(across(where(is_character), as_factor))
wc_long = wc_long %>% dplyr::mutate(across(where(is_integer), as_factor))
wc_long = left_join(wc_long, hc_new)

wc_long_trim = wc_long %>% select(file_name, Cage, Sex, Age, trialType, startDate, numMice, totWeight, avgWeight, Intake_24h, Intake_Std_tot, Intake_Std_avg,
                                  Time_cageType, Time, cageType, pct_cover)

wc_long_trim = wc_long_trim %>%
  mutate_at(vars(pct_cover), ~replace_na(., 0))

# What is the relationship between 24h water intake and urine cover (latrine cage only)?
wc_long_trim_sum = wc_long_trim %>% filter(trialType %in% c("Accl", "Record")) %>%
  group_by(Cage, Sex, Age, trialType, startDate, numMice, totWeight, Intake_24h) %>%
  dplyr::summarise(tot_pct_cover = sum(pct_cover),
                   n = n()) %>%
  ungroup()

# One Model per Sex:
summary(lmer(tot_pct_cover ~ Intake_24h + Age + (1|Cage), data = wc_long_trim_sum %>% filter(Sex == "F"))) # 0.0451 *
summary(lmer(tot_pct_cover ~ Intake_24h + Age + totWeight + (1|Cage), data = wc_long_trim_sum %>% filter(Sex == "F"))) # NS

summary(lmer(tot_pct_cover ~ Intake_24h + Age + (1|Cage), data = wc_long_trim_sum %>% filter(Sex == "M"))) # 0.0911 .
summary(lmer(tot_pct_cover ~ Intake_24h + Age + totWeight + (1|Cage), data = wc_long_trim_sum %>% filter(Sex == "M"))) # NS

# Is there an absolute difference in water intake between age classes? NO
summary(lmer(Intake_24h ~ Age + (1|Cage), data = wc_long_trim_sum %>% filter(Sex == "F"))) # 0.115
summary(lmer(Intake_24h ~ Age + totWeight + (1|Cage), data = wc_long_trim_sum %>% filter(Sex == "F"))) # NS

summary(lmer(Intake_24h ~ Age + (1|Cage), data = wc_long_trim_sum %>% filter(Sex == "M"))) # NS
summary(lmer(Intake_24h ~ Age + totWeight + (1|Cage), data = wc_long_trim_sum %>% filter(Sex == "M"))) # NS
