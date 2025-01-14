# Libraries
library(tidyverse)
library(ggbeeswarm)
library(plotly)
library(lme4)
library(lmerTest)

########## 90 minute Void Spot Assay ##########

# Set working directory and plot directory 
setwd("/Volumes/bedfordlab/Nocturia/FilterPaper/VSA")
plotdir = "/Volumes/bedfordlab/Nocturia/Data_Analysis/Plots"

# Read in Data:
list_of_files <- list.files(pattern = "*.csv", full.names = TRUE, recursive = TRUE) # Use recursive to access files in all subfolders
vsa <- list_of_files %>%
  setNames(nm = .) %>% 
  map_df(~read_csv(.x, col_types = cols(), col_names = FALSE), .id = "file_name") 
colnames(vsa) = c("file_name", "spot", "Label", "Area", "X", "Y")
vsa = vsa[!grepl("Area", vsa$Area), ]
vsa = vsa[!grepl("Drawing of duplicate", vsa$Label), ]
vsa = vsa %>% dplyr::select(-Label)
vsa$file_name = gsub(".*/", "", vsa$file_name)
vsa = vsa %>% mutate_at(vars(c(3:5)), function(x) as.numeric(as.character(x))) 
length(unique(vsa$file_name))

# Convert Area to Volume (uL): y = 12.7722x (see /Social Hierarchy/VSA_StandardCurve/OUTPUT/VSA_standardCurve.R)
vsa$Vol = 12.7722*(vsa$Area)

# Void spot histogram:
hist(vsa$Vol)

# Calculate summary values for each filter paper:
vs_sum = vsa %>% 
  filter(Vol >= 1) %>% # only consider voids larger than 1uL
  group_by(file_name) %>%
  dplyr::summarise(numVoids = n(),
            minVol = min(Vol),
            maxVol = max(Vol),
            meanVol = mean(Vol),
            medVol = median(Vol),
            totVol = sum(Vol)) %>% ungroup()

# Read in Lookup Tables:
mouse_lut = read.csv("/Volumes/bedfordlab/Nocturia/Data_Analysis/Nocturia Lookup Tables - LUT.csv")

meanWeight = as.data.frame(mouse_lut %>% group_by(Age, Sex) %>%
  dplyr::summarise(meanWeight = mean(Weight, na.rm = TRUE)))

# For mice with missing weights, assign mean value for Age/Sex class
mouse_lut = mouse_lut %>% mutate(Weight = case_when(is.na(Weight) == T & Age == "Old" & Sex == "F" ~ meanWeight[1,3],
                                                   is.na(Weight) == T & Age == "Old" & Sex == "M" ~ meanWeight[2,3],
                                                   is.na(Weight) == T & Age == "Young" & Sex == "M" ~ meanWeight[4,3],
                                                   .default = Weight))

vsa_lut = read.csv("/Volumes/bedfordlab/Nocturia/Data_Analysis/Nocturia Lookup Tables - VSA.csv")
vsa_lut = vsa_lut %>% mutate(Cohort = paste("C", Cohort, sep = ""))
vsa_lut = full_join(vsa_lut %>% select(-Weight), 
                mouse_lut %>% select(-DOB, -Sac_date, -Sac_time), 
                by = c("Round", "Cohort", "Age", "Cage", "Mouse_ID", "Sex"))
vsa_lut = vsa_lut %>% dplyr::rename(file_name = CSV)

# Justification for 90 minute assay:
vsa_lut_filt = vsa_lut %>% filter(!(is.na(X60m_voids))) %>% droplevels()
table(vsa_lut_filt$Sex, vsa_lut_filt$Time)
vsa_lut_filt = vsa_lut_filt %>% 
  mutate(across(c(X60m_voids, X75m_voids, X90m_voids), factor))

vsa_lut_filt_long = vsa_lut_filt %>% 
  select(Age, Cage, Mouse_ID, Sex, Time, Trial, X60m_voids, X75m_voids, X90m_voids) %>%
  pivot_longer(X60m_voids:X90m_voids, names_to = "TimePoint", values_to = "num_voids")

ggplot(vsa_lut_filt_long, aes(TimePoint, fill = num_voids)) + 
  geom_bar(position = position_fill(reverse = FALSE)) +
  scale_fill_manual(values=c("lightgrey", "#FCDBB4", "#FAB868", "#F7941D")) +
  xlab(NULL) + ylab(NULL) + 
  theme_classic() +
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank())
ggsave("FigS3.VSA_TrialDuration.pdf", width=4.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)

vsa_lut_filt_long$num_voids_bin = ifelse(vsa_lut_filt_long$num_voids == 0, 0, 1)

vsa_lut_filt_long %>% filter(TimePoint == "X60m_voids") %>% count("num_voids_bin") 
11/(11+9) # 55%

vsa_lut_filt_long %>% filter(TimePoint == "X75m_voids") %>% count("num_voids_bin") 
14/(14+6) # 70%

vsa_lut_filt_long %>% filter(TimePoint == "X90m_voids") %>% count("num_voids_bin") 
16/(16+4) # 80%
  
# Join dataframes:
vs = full_join(vsa_lut, vs_sum)

# Re-level factors:
vs = vs %>% mutate(across(where(is_character), as_factor))
vs$Round = as.factor(vs$Round)
vs = vs %>%
  mutate(Age = fct_relevel(Age, c("Young", "Old"))) %>%
  mutate(Time = fct_relevel(Time, c("light", "dark"))) %>%
  mutate(Cage = as.factor(Cage)) %>%
  mutate(trialDur = as.POSIXct(trialDur, format = "%H:%M")) %>%
  mutate(across(Start:End, ~ as.POSIXct(., format = "%H:%M:%S")))
str(vs)

vs_new = vs %>% 
  dplyr::select(Round, Cohort, Age, Sex, Cage, Mouse_ID, Weight, Trial, Date, Time, Start, End, trialDur, cFOS, VSA_control, arenaSize, arenaColor, Fecal, fileName, file_name, TrueZero, Chewed,
         numVoids, minVol, maxVol, meanVol, medVol, totVol) %>%
  dplyr::mutate(across(.cols = c(numVoids, minVol, maxVol, meanVol, medVol, totVol),
                ~if_else(is.na(numVoids) == TRUE, 0, .)))

p1 = ggplot(vs_new, aes(Time, maxVol, color = Round)) +
  facet_grid(Sex ~ Age) +
  geom_beeswarm(size = 1.25, cex = 3)
ggplotly(p1)

p2 = ggplot(vs_new, aes(Time, totVol, color = VSA_control)) +
  facet_grid(Sex ~ Age) +
  geom_beeswarm(size = 1.25, cex = 3)
ggplotly(p2)

# Calculate mean values per Mouse:
vs_avg = vs_new %>%
  group_by(Round, Cohort, Age, Sex, Cage, Mouse_ID, Time) %>%
  dplyr::summarise(meanVol = mean(meanVol),
            medVol = mean(medVol),
            numVoids = mean(numVoids),
            totVol = mean(totVol),
            maxVol = max(maxVol),
            n = n())

p3 = ggplot(vs_avg, aes(Time, totVol, color = Cohort)) +
  facet_grid(Sex ~ Age) +
  geom_point(size = 1.25) +
  ylab("total urine volume (uL)") +
  geom_line(aes(group=interaction(Mouse_ID)))
ggplotly(p3)

table(vs_avg$Sex, vs_avg$Age, vs_avg$Round)

## First pass statistical models
## Make new filtered data set
## Possible filters: VSA control cages, cFos trials
vs_filt = vs_new %>% 
  filter(VSA_control != "Y") %>% 
  filter(cFOS != "Y")

table(vs_filt$Sex, vs_filt$Age, vs_filt$Round)

## Male Models:
# Interaction model:
summary(lmer(totVol ~ Time*Age + (1|Cage/Mouse_ID), data = vs_filt %>% filter(Sex == "M"))) # interaction NS
summary(lmer(maxVol ~ Time*Age + (1|Cage/Mouse_ID), data = vs_filt %>% filter(Sex == "M"))) # interaction NS

# Young Males
m1_Y = lmer(totVol ~ Time + (1|Cage/Mouse_ID), 
             data = vs_filt %>% filter(Sex == "M" & Age == "Young"))
summary(m1_Y) # Effect size = 137.154 (nested RE)

m2_Y = lmer(maxVol ~ Time + (1|Cage/Mouse_ID), 
            data = vs_filt %>% filter(Sex == "M" & Age == "Young"))
summary(m2_Y) # Effect size = 115.381 (nested RE)

# Old Males
m1_O = lmer(totVol ~ Time + (1|Cage/Mouse_ID), 
           data = vs_filt %>% filter(Sex == "M" & Age == "Old"))
summary(m1_O) # Effect size = 102.048 (nested RE)

m2_O = lmer(maxVol ~ Time + (1|Cage/Mouse_ID), 
            data = vs_filt %>% filter(Sex == "M" & Age == "Old"))
summary(m2_O) # Effect size = 93.05 (nested RE)

## Female Models:
# Interaction model:
summary(lmer(totVol ~ Time*Age + (1|Cage/Mouse_ID), data = vs_filt %>% filter(Sex == "F"))) # interaction NS
summary(lmer(maxVol ~ Time*Age + (1|Cage/Mouse_ID), data = vs_filt %>% filter(Sex == "F"))) # interaction NS

# Young Females
f1_Y = lmer(totVol ~ Time + (1|Cage/Mouse_ID), 
            data = vs_filt %>% filter(Sex == "F" & Age == "Young"))
summary(f1_Y) # NS

f2_Y = lmer(maxVol ~ Time + (1|Cage/Mouse_ID), 
            data = vs_filt %>% filter(Sex == "F" & Age == "Young"))
summary(f2_Y) # Effect size = 73.26

# Old Females
f1_O = lmer(totVol ~ Time + (1|Cage/Mouse_ID), 
            data = vs_filt %>% filter(Sex == "F" & Age == "Old"))
summary(f1_O) # NS

f2_O = lmer(maxVol ~ Time + (1|Cage/Mouse_ID), 
            data = vs_filt %>% filter(Sex == "F" & Age == "Old"))
summary(f2_O) # Effect size = 108.06

# Calculate mean values per Cage:
vs_filt_avg = vs_filt %>%
  group_by(Round, Cohort, Age, Sex, Cage, Mouse_ID, Time) %>%
  dplyr::summarise(meanVol = mean(meanVol),
            medVol = mean(medVol),
            numVoids = mean(numVoids),
            totVol = mean(totVol),
            maxVol = max(maxVol),
            n = n()) %>%
  ungroup() %>%
  droplevels()

# Sample size:
range(vs_filt_avg$n) # 2-3
round(mean(vs_filt_avg$n), digits = 1) # 2.1

# How many mice don't have matching light and dark trials? 
pairs = as.data.frame(table(vs_filt_avg$Mouse_ID))
colnames(pairs)[1] = "Mouse_ID"
goodPairs = pairs %>% filter(Freq == 2) %>% droplevels() %>% select(Mouse_ID)
goodPairs = goodPairs[["Mouse_ID"]]
goodPairs # 74

# How many mice showed increased toVol or maxVol from light to dark?
vs_filt_avg_wide = vs_filt_avg %>% filter(Mouse_ID %in% goodPairs) %>%
  pivot_wider(names_from = Time, values_from = meanVol:n) %>%
  mutate(totVol_diff = totVol_dark - totVol_light) %>% # total difference
  mutate(maxVol_diff = maxVol_dark - maxVol_light) %>%
  mutate(totVol_pctDiff = (totVol_dark - totVol_light)/totVol_dark) %>% # percent difference
  mutate(maxVol_pctDiff = (maxVol_dark - maxVol_light)/maxVol_dark) %>%
  mutate(totVol_diffBin = case_when(totVol_diff > 0 ~ "increase",
                                    totVol_diff < 0 ~ "decrease")) %>%
  mutate(maxVol_diffBin = case_when(maxVol_diff > 0 ~ "increase",
                                    maxVol_diff < 0 ~ "decrease"))

vs_filt_avg_wide[mapply(is.infinite, vs_filt_avg_wide)] = NA

# Fisher's Exact Test
table(vs_filt_avg_wide$Age, vs_filt_avg_wide$totVol_diffBin, vs_filt_avg_wide$Sex)
table(vs_filt_avg_wide$Age, vs_filt_avg_wide$maxVol_diffBin, vs_filt_avg_wide$Sex)

# Make Plots:
setwd(plotdir)
library(Rmisc)

# Calculate 'effect size' for totVol and maxVol:
ES_totVol = summarySE(data = vs_filt_avg_wide, 
                      measurevar = "totVol_diff", na.rm = T, 
                      groupvars = c("Sex", "Age"))

round(ES_totVol$totVol_diff, digits = 1)

ES_maxVol = summarySE(data = vs_filt_avg_wide, 
                      measurevar = "maxVol_diff", na.rm = T, 
                      groupvars = c("Sex", "Age"))

round(ES_maxVol$maxVol_diff, digits = 0)

# Intra- and inter-individual variance:
var.test(totVol ~ Age, data = vs_filt %>% filter (Sex == "F")) # p-value = 0.02924
var.test(totVol ~ Age, data = vs_filt %>% filter (Sex == "M")) # p-value = 0.07353

# Inter-individual variance only:
summarySE(vs_filt_avg, measurevar = "totVol", groupvars = c("Sex", "Age"))
var.test(totVol ~ Age, data = vs_filt_avg %>% filter (Sex == "F")) # p-value = 0.01091
var.test(totVol ~ Age, data = vs_filt_avg %>% filter (Sex == "M")) # p-value = 0.7706

# Intra- and inter-individual variance:
var.test(maxVol ~ Age, data = vs_filt %>% filter (Sex == "F")) # p-value = 8.967e-07
var.test(maxVol ~ Age, data = vs_filt %>% filter (Sex == "M")) # p-value = 0.09954

# Inter-individual variance only:
summarySE(vs_filt_avg, measurevar = "maxVol", groupvars = c("Sex", "Age"))
var.test(maxVol ~ Age, data = vs_filt_avg %>% filter (Sex == "F")) # p-value = 0.0001497
var.test(maxVol ~ Age, data = vs_filt_avg %>% filter (Sex == "M")) # p-value = 0.8295

# Does total urine volume decrease overall with age?

summaryM = summarySE(vs_filt %>% filter(Sex == "M"), measurevar="totVol", groupvars="Age", na.rm=T)
ggplot(data=vs_filt %>% filter(Sex == "M"), aes(x=Age, y=totVol)) +
  geom_bar(data=summaryM, aes(fill=Age), stat="identity", width=0.85) +
  geom_point(aes(fill=Age), colour="black", pch=21, position=position_jitter(0.1), size=1) +
  geom_errorbar(data=summaryM, aes(ymin=totVol-se, ymax=totVol+se, x=Age), linewidth=0.5, width=0.25) + 
  scale_fill_manual(values=c("#6c8c24", "#4d5c2c")) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
#ggsave("VSA_totVol_M_allData.pdf", width=3, height=6, units="cm", dpi=1500, useDingbats=FALSE)

summary(lmer(totVol ~ Age + (1|Cage/Mouse_ID), data = vs_filt %>% filter(Sex == "M"))) # 0.182

summaryF = summarySE(vs_filt %>% filter(Sex == "F"), measurevar="totVol", groupvars="Age", na.rm=T)
ggplot(data=vs_filt %>% filter(Sex == "F"), aes(x=Age, y=totVol)) +
  geom_bar(data=summaryF, aes(fill=Age), stat="identity", width=0.85) +
  geom_point(aes(fill=Age), colour="black", pch=21, position=position_jitter(0.1), size=1) +
  geom_errorbar(data=summaryF, aes(ymin=totVol-se, ymax=totVol+se, x=Age), linewidth=0.5, width=0.25) + 
  scale_fill_manual(values=c("#8a1e8e", "#642c66")) +
  theme_classic() +  
  theme(legend.position = "none", axis.text.y=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
#ggsave("VSA_totVol_F_allData.pdf", width=3, height=6, units="cm", dpi=1500, useDingbats=FALSE)

summary(lmer(totVol ~ Age + (1|Cage/Mouse_ID), data = vs_filt %>% filter(Sex == "F"))) # 0.0471 *

# Does max void volume decrease overall with age?

summaryM = summarySE(vs_filt %>% filter(Sex == "M"), measurevar="maxVol", groupvars="Age", na.rm=T)
ggplot(data=vs_filt %>% filter(Sex == "M"), aes(x=Age, y=maxVol)) +
  geom_bar(data=summaryM, aes(fill=Age), stat="identity", width=0.85) +
  geom_point(aes(fill=Age), colour="black", pch=21, position=position_jitter(0.1), size=1) +
  geom_errorbar(data=summaryM, aes(ymin=maxVol-se, ymax=maxVol+se, x=Age), linewidth=0.5, width=0.25) + 
  scale_fill_manual(values=c("#6c8c24", "#4d5c2c")) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
#ggsave("VSA_maxVol_M_allData.pdf", width=3, height=6, units="cm", dpi=1500, useDingbats=FALSE)

summary(lmer(maxVol ~ Age + (1|Cage/Mouse_ID), data = vs_filt %>% filter(Sex == "M"))) # 0.368

summaryF = summarySE(vs_filt %>% filter(Sex == "F"), measurevar="maxVol", groupvars="Age", na.rm=T)
ggplot(data=vs_filt %>% filter(Sex == "F"), aes(x=Age, y=maxVol)) +
  geom_bar(data=summaryF, aes(fill=Age), stat="identity", width=0.85) +
  geom_point(aes(fill=Age), colour="black", pch=21, position=position_jitter(0.1), size=1) +
  geom_errorbar(data=summaryF, aes(ymin=maxVol-se, ymax=maxVol+se, x=Age), linewidth=0.5, width=0.25) + 
  scale_fill_manual(values=c("#8a1e8e", "#642c66")) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
#ggsave("VSA_maxVol_F_allData.pdf", width=3, height=6, units="cm", dpi=1500, useDingbats=FALSE)

summary(lmer(maxVol ~ Age + (1|Cage/Mouse_ID), data = vs_filt %>% filter(Sex == "F"))) # 0.407

# Does total urine volume during light and dark phase change with age? (see models above)

ggplot(data=vs_filt_avg %>% filter(Sex == "M" & Mouse_ID %in% goodPairs), aes(Time, totVol, color = Age)) +
  facet_wrap(~Age) +
  geom_line(aes(group=interaction(Mouse_ID))) +
  geom_point(size = 1.25) +
  scale_color_manual(values=c("#6c8c24", "#4d5c2c")) +
  ylim(0, 1308) + 
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Fig3.VSA_totVol_M.pdf", width=6, height=6, units="cm", dpi=1500, useDingbats=FALSE)

ggplot(data=vs_filt_avg %>% filter(Sex == "F" & Mouse_ID %in% goodPairs), aes(Time, totVol, color = Age)) +
  facet_wrap(~Age) +
  geom_line(aes(group=interaction(Mouse_ID))) +
  geom_point(size = 1.25) +
  scale_color_manual(values=c("#8a1e8e", "#642c66")) +
  ylim(0, 1308) + 
  xlab(NULL) + ylab(NULL) +
  theme_classic() +
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Fig3.VSA_totVol_F.pdf", width=6, height=6, units="cm", dpi=1500, useDingbats=FALSE)

# Does max void volume during light and dark phase change with age? (see models above)

ggplot(data=vs_filt_avg %>% filter(Sex == "M" & Mouse_ID %in% goodPairs), aes(Time, maxVol, color = Age)) +
  facet_wrap(~Age) +
  geom_line(aes(group=interaction(Mouse_ID))) +
  geom_point(size = 1.25) +
  scale_color_manual(values=c("#6c8c24", "#4d5c2c")) +
  ylim(0, 886) + 
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Fig3.VSA_maxVol_M.pdf", width=6, height=6, units="cm", dpi=1500, useDingbats=FALSE)

ggplot(data=vs_filt_avg %>% filter(Sex == "F" & Mouse_ID %in% goodPairs), aes(Time, maxVol, color = Age)) +
  facet_wrap(~Age) +
  geom_line(aes(group=interaction(Mouse_ID))) +
  geom_point(size = 1.25) +
  scale_color_manual(values=c("#8a1e8e", "#642c66")) +
  ylim(0, 886) + 
  xlab(NULL) + ylab(NULL) +
  theme_classic() +
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Fig3.VSA_maxVol_F.pdf", width=6, height=6, units="cm", dpi=1500, useDingbats=FALSE)

########## 12 hour Home Cage Assay ##########

# Set working directory and plot directory 
setwd("/Volumes/bedfordlab/Nocturia/FilterPaper/HomeCage")
plotdir = "/Volumes/bedfordlab/Nocturia/Data_Analysis/Plots" 

# Read in Data:
list_of_files <- list.files(pattern = "*.csv", full.names = TRUE, recursive = TRUE) # Use recursive to access files in all subfolders
hc <- list_of_files %>%
  setNames(nm = .) %>% 
  map_df(~read_csv(.x, col_types = cols(), col_names = FALSE), .id = "file_name") 
colnames(hc) = hc[1, ]
colnames(hc)[1] = "file_name"
colnames(hc)[2] = "spot"
hc = hc[!grepl("Area", hc$Area), ]
hc = hc[!grepl("Drawing of duplicate", hc$Label), ]
hc = hc %>% filter(MinThr != 255)
hc = hc %>% filter(Circ. > 0.1)
hc$file_name = gsub(".*/", "", hc$file_name)
hc = hc %>% mutate_at(vars(c(4:10)), function(x) as.numeric(as.character(x))) 
hc = hc %>% mutate(Label = case_when(spot == 1 ~ "FP_Area", .default = NA))
hc = hc[ ,-c(5:12)]
length(unique(hc$file_name)) # 338

# Calculate summary values for each filter paper:
hc_area = hc %>% filter(Label == "FP_Area")
colnames(hc_area)[4] = "FP_Area"
hc_sum = hc %>% filter(!spot == 1) %>%
  group_by(file_name) %>%
  dplyr::summarise(numVoids = n(),
            totUrine = sum(Area)) %>% ungroup()

hc_sum = full_join(hc_sum, hc_area) %>% 
  dplyr::select(-spot, -Label)

# Read in Lookup Table:
hc_lut = read.csv("/Volumes/bedfordlab/Nocturia/Data_Analysis/Nocturia Lookup Tables - HomeCage.csv")
str(hc_lut)
length(unique(hc_lut$fileName)) # 375
length(unique(hc_lut$CSV)) # 339
hc_lut = hc_lut %>% dplyr::rename(file_name = CSV)

# Join dataframes:
hc_new = full_join(hc_lut, hc_sum)
str(hc_new)

meanArea = hc_new %>% filter(Chewed != "Y") %>% select(FP_Area) %>% unlist()
meanArea = mean(meanArea, na.rm = TRUE)

hc_new = hc_new %>%
  dplyr::select(Round, Age, Sex, Cage, numMice, Trial, segment, Time, endDate, endTime, trialDur, cageType, trialType, TrueZero, Chewed, file_name,
                numVoids, totUrine, FP_Area) %>%
  dplyr::mutate(across(.cols = c(numVoids, totUrine),
                       ~if_else(TrueZero == "Y", 0, .))) %>%
  dplyr::mutate(FP_Area = case_when(TrueZero == "Y" & Chewed != "Y" ~ meanArea,
                                 .default = FP_Area))

hc_new = hc_new %>% mutate_if(is.character, as.factor)
hc_new = as_tibble(hc_new)

# Calculate % Cover
str(hc_new)
hc_new = hc_new %>%
  dplyr::mutate(pct_cover = case_when(TrueZero != "Y" & Chewed != "Y" ~ ((totUrine / FP_Area)*100),
                                      TrueZero != "Y" & Chewed == "Y" & FP_Area >= meanArea ~ ((totUrine / FP_Area)*100), # smallChew
                                      TrueZero != "Y" & Chewed == "Y" & FP_Area < meanArea ~ ((totUrine + (meanArea - FP_Area)) / meanArea)*100, # bigChew
                                      TrueZero == "Y" & Chewed == "Y" ~ ((meanArea - FP_Area) / meanArea)*100,
                                      .default = 0))

hc_new = hc_new %>% 
  dplyr::mutate(ChewSize = case_when(TrueZero != "Y" & Chewed == "Y" & FP_Area < meanArea ~ "bigChew",
                                     TrueZero != "Y" & Chewed == "Y" & FP_Area >= meanArea ~ "smallChew"))

# Proportion of Chewed Papers:
chew = hc_new %>% group_by(Age, Sex) %>%
  dplyr::summarise(num_bigChew = sum(ChewSize == "bigChew", na.rm = T), n = n())
chew$pct_chew = (chew$num_bigChew / chew$n)*100

# Re-level factors:
hc_new = hc_new %>%
  mutate(Age = fct_relevel(Age, c("Young", "Old"))) %>%
  mutate(Time = fct_relevel(Time, c("light", "dark")))

# There are two FPs for 08182022_792 L2. Unclear which is correct --> arbitrarily remove 20220818_18.06_C2_792_M_latrine.tiff
fo = hc_new %>% filter(Trial == "08182022_792" & segment == "L2") # Values are very similar
hc_new = hc_new %>% filter(file_name != "20220818_18.06_C2_792_M_latrine_ps.csv")

# Save output file (This file is used for further analysis with Ethovision data)
setwd("/Volumes/bedfordlab/Nocturia/Data_Analysis")
write.csv(hc_new, "HCLC_filterPapers.csv")

# Remove 48h trials:
hc_new = hc_new %>% filter(trialDur != 48)

# Plot the raw data:
setwd(plotdir)

ggplot(hc_new %>% filter(trialType %in% c("Accl", "Record") & Sex == "M"), 
       aes(Time, pct_cover, color = Age)) +
  facet_grid(~Age) +
  geom_beeswarm(alpha = 0.5, size = 1, cex = 2.5, method = "hex") +
  ggtitle("Male LC Papers")
#ggsave("LatrineCageFP_M_allData.pdf", width=15, height=10, units="cm", dpi=1500, useDingbats=FALSE)

ggplot(hc_new %>% filter(trialType %in% c("Accl", "Record") & Sex == "F"), 
       aes(Time, pct_cover, color = Age)) +
  facet_grid(~Age) +
  geom_beeswarm(alpha = 0.5, size = 1, cex = 2.5, method = "hex") +
  ggtitle("Female LC Papers")
#ggsave("LatrineCageFP_F_allData.pdf", width=15, height=10, units="cm", dpi=1500, useDingbats=FALSE)

## First Pass Statistical Models:

# Male Models:
m1 = lmer(pct_cover ~ Time*Age + (1|Cage), hc_new %>% filter(trialType %in% c("Accl", "Record") & Sex == "M"))
summary(m1)

m1_Y = lmer(pct_cover ~ Time + (1|Cage), 
            data = hc_new %>% filter(trialType %in% c("Accl", "Record") & Sex == "M" & Age == "Young"))
summary(m1_Y) # p = 0.00202 **; Effect size = 7.288

m1_O = lmer(pct_cover ~ Time + (1|Cage), 
            data = hc_new %>% filter(trialType %in% c("Accl", "Record") & Sex == "M" & Age == "Old"))
summary(m1_O) # p = 0.000138 ***; Effect size = 10.276

## Male Control Experiment Models:

# Interaction:
summary(lmer(pct_cover ~ Time*cageType + (1|Cage),
             data = hc_new %>% filter(trialType %in% c("Ctrl") & Sex == "M" & Age == "Young"))) # cageType: 0.0130 *
summary(lmer(pct_cover ~ Time*cageType + (1|Cage),
             data = hc_new %>% filter(trialType %in% c("Ctrl") & Sex == "M" & Age == "Old"))) # cageType: NS

summary(lmer(pct_cover ~ Time + (1|Cage),
             data = hc_new %>% filter(trialType %in% c("Ctrl") & Sex == "M" & Age == "Young" & cageType == "home"))) # p = 0.042 *
summary(lmer(pct_cover ~ Time + (1|Cage),
             data = hc_new %>% filter(trialType %in% c("Ctrl") & Sex == "M" & Age == "Young" & cageType == "latrine"))) # p = 0.00287 **
summary(lmer(pct_cover ~ Time + (1|Cage),
             data = hc_new %>% filter(trialType %in% c("Ctrl") & Sex == "M" & Age == "Old" & cageType == "home"))) # 0.04883 *
summary(lmer(pct_cover ~ Time + (1|Cage),
             data = hc_new %>% filter(trialType %in% c("Ctrl") & Sex == "M" & Age == "Old" & cageType == "latrine"))) # NS

# Female Models:
f1 = lmer(pct_cover ~ Time*Age + (1|Cage), hc_new %>% filter(trialType %in% c("Accl", "Record") & Sex == "F"))
summary(f1)

f1_Y = lmer(pct_cover ~ Time + (1|Cage), 
            data = hc_new %>% filter(trialType %in% c("Accl", "Record") & Sex == "F" & Age == "Young"))
summary(f1_Y) # p = 0.00242 **; Effect size = 3.421

f1_O = lmer(pct_cover ~ Time + (1|Cage), 
            data = hc_new %>% filter(trialType %in% c("Accl", "Record") & Sex == "F" & Age == "Old"))
summary(f1_O) # p = 0.0029 **; Effect size = 6.402

## Female Control Experiment Models:

# Interaction:
summary(lmer(pct_cover ~ Time*cageType + (1|Cage),
             data = hc_new %>% filter(trialType %in% c("Ctrl") & Sex == "F" & Age == "Young"))) # cageType: 0.00211 **
summary(lmer(pct_cover ~ Time*cageType + (1|Cage),
             data = hc_new %>% filter(trialType %in% c("Ctrl") & Sex == "F" & Age == "Old"))) # cageType: NS

summary(lmer(pct_cover ~ Time + (1|Cage),
             data = hc_new %>% filter(trialType %in% c("Ctrl") & Sex == "F" & Age == "Young" & cageType == "home"))) # p = 0.00465 **
summary(lmer(pct_cover ~ Time + (1|Cage),
             data = hc_new %>% filter(trialType %in% c("Ctrl") & Sex == "F" & Age == "Young" & cageType == "latrine"))) # p = 0.0223 *
summary(lmer(pct_cover ~ Time + (1|Cage),
             data = hc_new %>% filter(trialType %in% c("Ctrl") & Sex == "F" & Age == "Old" & cageType == "home"))) # NS
summary(lmer(pct_cover ~ Time + (1|Cage),
             data = hc_new %>% filter(trialType %in% c("Ctrl") & Sex == "F" & Age == "Old" & cageType == "latrine"))) # p = 0.024165 *

# Calculate mean values per Cage, per trialType:
hc_avg = hc_new %>% filter(trialType %in% c("Accl", "Record")) %>%
  group_by(Age, Sex, Cage, Time, cageType) %>%
  dplyr::summarise(mean_pct_cover = mean(pct_cover),
            n = n())

# Sample size:
hc_avg_tot = hc_avg %>% group_by(Cage) %>%
  group_by(Age, Sex, Cage, cageType) %>%
  dplyr::summarise(mean_pct_cover = mean(mean_pct_cover),
                   n = sum(n))

range(hc_avg_tot$n) # 4-15
round(mean(hc_avg_tot$n), digits = 1) # 8.7

# Make Plots:
ggplot(data=hc_avg %>% filter(Sex == "M"), 
       aes(Time, mean_pct_cover, color = Age)) +
  facet_wrap(~Age) +
  geom_line(aes(group=interaction(Cage))) +
  geom_point(size = 1.25) +
  scale_color_manual(values=c("#6c8c24", "#4d5c2c")) +
  ylim(0, 60) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("FigS2.pct_cover_Accl_Record_M.pdf", width=6, height=6, units="cm", dpi=1500, useDingbats=FALSE)

ggplot(data=hc_avg %>% filter(Sex == "F"), 
       aes(Time, mean_pct_cover, color = Age)) +
  facet_wrap(~Age) +
  geom_line(aes(group=interaction(Cage))) +
  geom_point(size = 1.25) +
  scale_color_manual(values=c("#8a1e8e", "#642c66")) +
  ylim(0, 60) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("FigS2.pct_cover_Accl_Record_F.pdf", width=6, height=6, units="cm", dpi=1500, useDingbats=FALSE)

# Plot the control trial data:
# Calculate mean values per Cage:
hc_avg2 = hc_new %>% filter(trialType == "Ctrl") %>%
  group_by(Age, Sex, Cage, cageType, Time) %>%
  dplyr::summarise(mean_pct_cover = mean(pct_cover),
            n = n())

ggplot(hc_avg2 %>% filter(Sex == "F"), 
       aes(Time, mean_pct_cover, color = Age)) +
  facet_grid(Age ~ cageType) +
  geom_point(size = 1.25) +
  geom_line(aes(group=interaction(Cage))) +
  ggtitle("Female HC LC Papers")
#ggsave("CtrlTrials_FP_F_allData.pdf", width=12, height=10, units="cm", dpi=1500, useDingbats=FALSE)

ggplot(data=hc_avg2 %>% filter(Sex == "F" & Age == "Young"), 
       aes(Time, mean_pct_cover, color = cageType)) +
  facet_wrap(~cageType) +
  geom_line(aes(group=interaction(Cage))) +
  geom_point(size = 1.25) +
  scale_color_manual(values=c("#00a08a", "#c6340a")) +
  ylim(0, 60) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("FigS2.pct_cover_Ctrl_F_Young.pdf", width=6, height=6, units="cm", dpi=1500, useDingbats=FALSE)

ggplot(data=hc_avg2 %>% filter(Sex == "F" & Age == "Old"), 
       aes(Time, mean_pct_cover, color = cageType)) +
  facet_wrap(~cageType) +
  geom_line(aes(group=interaction(Cage))) +
  geom_point(size = 1.25) +
  scale_color_manual(values=c("#00a08a", "#c6340a")) +
  ylim(0, 60) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("FigS2.pct_cover_Ctrl_F_Old.pdf", width=6, height=6, units="cm", dpi=1500, useDingbats=FALSE)

ggplot(hc_avg2 %>% filter(Sex == "M"), 
       aes(Time, mean_pct_cover, color = Age)) +
  facet_grid(Age ~ cageType) +
  geom_point(size = 1.25) +
  geom_line(aes(group=interaction(Cage))) +
  ggtitle("Male HC LC Papers")
#ggsave("CtrlTrials_FP_M_allData.pdf", width=12, height=10, units="cm", dpi=1500, useDingbats=FALSE)

ggplot(data=hc_avg2 %>% filter(Sex == "M" & Age == "Young"), 
       aes(Time, mean_pct_cover, color = cageType)) +
  facet_wrap(~cageType) +
  geom_line(aes(group=interaction(Cage))) +
  geom_point(size = 1.25) +
  scale_color_manual(values=c("#00a08a", "#c6340a")) +
  ylim(0, 60) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("FigS2.pct_cover_Ctrl_M_Young.pdf", width=6, height=6, units="cm", dpi=1500, useDingbats=FALSE)

ggplot(data=hc_avg2 %>% filter(Sex == "M" & Age == "Old"), 
       aes(Time, mean_pct_cover, color = cageType)) +
  facet_wrap(~cageType) +
  geom_line(aes(group=interaction(Cage))) +
  geom_point(size = 1.25) +
  scale_color_manual(values=c("#00a08a", "#c6340a")) +
  ylim(0, 60) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("FigS2.pct_cover_Ctrl_M_Old.pdf", width=6, height=6, units="cm", dpi=1500, useDingbats=FALSE)