require(tidyverse)
require(RSQLite)
require(DBI)

# Set working directory and plot directory (Skip to Line #144 to load R Workspace with fully processed dataframe)
setwd("/Users/nicole/Desktop/NocEtho_Outputs")
plotdir <- paste(getwd(), "/Plots/", sep = "") 

# set up a SQLite database to contain the processed files. Maybe this will make it easier to share and retrieve data on demand?
dbConn<- dbConnect(SQLite(), "nocturia_db.sqlite") # might be better to make this be local on your machine

# Read in lookup table: https://docs.google.com/spreadsheets/d/1f7IaATHpR9399gbgsBQUwJI2blYPv9oDj9bNo60KNEo/edit?gid=266372834#gid=266372834
lut <- read_csv("/Volumes/bedfordlab/Nocturia/Data_Analysis/Nocturia Lookup Tables - LUT.csv")
str(lut)
lut <- lut %>% mutate_at(c("Round", "Cohort", "Age", "Cage", "Mouse_ID", "Sex"), as.factor)
lut = lut %>% mutate_at(vars(DOB, Sac_date), as.Date, format = "%m/%d/%Y")
lut$Sac_time <- as.POSIXct(lut$Sac_time, format = "%I:%M:%S %p")
str(lut)

# PART 1: DATA MASSAGING
# List file names and sort them in order, so they're in the right sequence

getwd()
etho.files <- dir(pattern = ".txt$") 

# set the column names given that they are actually in 2 rows instead of only 1
cnames<- c("Trialtime_s", "Recordingtime_s", "Xcenter_cm", "Ycenter_cm", "Area_cm2", "Areachange_cm2", "Elongation", "Distancemoved_cm", "Velocity_cm_s", "Inzone", "Activity_pct", "Result1")

# loop through the files, write each output to db
for(i in 1:length(etho.files)){
  # columns appear to not always be in the same order, but we have a two row column name, try and work around that
  test<- read_csv(etho.files[i], skip = 34, n_max = 1)
  # get the start time
  st<- str_split(read_lines(etho.files[i], skip = 12, n_max = 1), ",", simplify = TRUE)[1,2] %>% 
    str_replace_all(pattern = "\"", replacement = "") %>% 
    as.POSIXct(format = "%m/%d/%Y %H:%M:%OS", tz = "America/Denver")
  
  # normalized start time
  nt<- as.POSIXct(paste("2000-01-01", format(st, format = "%H:%M:%OS")), tz = "America/Denver")
  
  # read in the data
  dat<- read_csv(etho.files[i], skip = 35, na = c("", "-"))
  names(dat)<- names(test)
  
  # do some changes/calculations to the data.frame
  dat<- dat %>% 
    rename(Trialtime_s = `Trial time`,
           Recordingtime_s = `Recording time`,
           Xcenter_cm = `X center`,
           Ycenter_cm = `Y center`,
           Area_cm2 = Area,
           Areachange_cm2 = Areachange,
           Distancemoved_cm = `Distance moved`,
           Velocity_cm_s = Velocity,
           Inzone = `In zone`,
           Activity_pct = Activity,
           Result1 = `Result 1`) %>% 
    select(Trialtime_s, Recordingtime_s, Xcenter_cm, Ycenter_cm, Area_cm2, Areachange_cm2, Elongation, Distancemoved_cm, Velocity_cm_s, Inzone, Activity_pct, Result1) %>% 
    mutate(Timestamp = Recordingtime_s + st,
           TimestampNormal = Recordingtime_s + nt,
           FilePath = str_replace(etho.files[i], pattern = ".txt", replacement = ""),
           Inzone = as.integer(Inzone),
           Result1 = as.integer(Result1)) %>% 
    mutate(Datestamp = date(Timestamp),
           LD = case_when(hour(Timestamp) >= 6 & hour(Timestamp) < 18 ~ "light",
                          TRUE ~ "dark"),
           TimestampRound = round_date(Timestamp, "minute"),
           TimestampNormalRound = round_date(TimestampNormal, "minute"))  %>% 
    select(Activity_pct, FilePath, Timestamp, Datestamp, LD, TimestampNormal, TimestampNormalRound) %>% 
    separate_wider_delim(FilePath, delim = "_", names = c("startDate", "startTime", "Cam", "Cohort", "Cage", "CageType"), cols_remove = FALSE) %>% 
    mutate(startDate = case_when(startDate == "08172022" ~ "08182022",
                                 TRUE ~ startDate),
           Timestamp = as.character(Timestamp), # need to coerce datetime to character since SQLite doesn't have a datetime data type
           Datestamp = as.character(Datestamp),
           TimestampNormal = as.character(TimestampNormal),
           TimestampNormalRound = as.character(TimestampNormalRound))
  
  # save the data into central database
  dbWriteTable(dbConn, name = "tblRawProccessed", value = dat, append = TRUE)
  
  rm(dat)
  gc()
  
  print(paste("Done with file:", etho.files[i]))
  
}

# release the database connection
dbDisconnect(dbConn)

# Try query the db on local disk, see if performance increases...The answer is absolutely have the db local, the IO is much improved
# This is the same database, I just copied from Alcova
dbConn<- dbConnect(SQLite(), "/Users/nicole/Desktop/NocEtho_Outputs/nocturia_db.sqlite")

# get some metadata to help describe all of the data and store into a table
fls<- dbGetQuery(dbConn, "SELECT FilePath, startDate, startTime, Cam, Cohort, Cage, CageType, 
                 MIN(Timestamp) AS MinTimestamp, 
                 MAX(Timestamp) AS MaxTimestamp FROM tblRawProccessed GROUP BY FilePath, 
                 startDate, startTime, Cam, Cohort, Cage, CageType")

# Fix FilePath typos 
length(unique(fls$FilePath))

# Save this to the database
dbWriteTable(dbConn, name = "tblExperiments", value = fls, overwrite = TRUE)

dbDisconnect(dbConn)

# PART 2: DATA ANALYSIS
dbConn<- dbConnect(SQLite(), "/Users/nicole/Desktop/NocEtho_Outputs/nocturia_db.sqlite")

expmnts<- dbReadTable(dbConn, "tblExperiments")

dbListTables(dbConn)
dbListFields(dbConn, "tblRawProccessed")

dbGetQuery(dbConn, "SELECT COUNT(*) AS cnt FROM tblRawProccessed")

## Step 1: Remove activity "outliers" (convert activity values > 30 to NA)
## Step 2: Calculate "mean" activity per 1 minute bin
## Step 3: Normalize "per minute mean" activity (per trial)
## Step 4: Calculate cumulative "normalized" activity per 12 hour bin (per trial)

# Use this query to avoid reading entirety of dataset into RAM prior to summarizing/plotting
allTrials<- dbGetQuery(dbConn, "SELECT FilePath, startDate, Cam, Cohort, Cage, CageType, TimestampNormalRound, AVG(Activity_pct) AS min_avg FROM tblRawProccessed WHERE Activity_pct < 30 GROUP BY FilePath, startDate, Cam, Cohort, Cage, CageType, TimestampNormalRound") %>%
  mutate(TimestampNormalRound = as.POSIXct(TimestampNormalRound, format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")) %>% 
  group_by(FilePath) %>%
  mutate(z_score = scale(min_avg, center = TRUE, scale = TRUE)) %>% 
  mutate(mm_norm = (min_avg - min(min_avg)) / (max(min_avg) - min(min_avg))) %>%
  unite("Trial", c("startDate", "Cage"), remove = FALSE) %>%
  ungroup() %>%  
  unite("Trial.CageType", c("startDate", "Cage", "CageType"), remove = FALSE) %>%
  ungroup() %>% 
  left_join(lut %>% 
              select(Round:numMice, Sex) %>% 
              mutate(Cage = as.character(Cage)) %>% 
              distinct(), by = c("Cohort", "Cage"))

# ####### Save the Workspace #######
# setwd("/Volumes/bedfordlab/Nocturia/Data_Analysis")
# save.image(file = "latrineNoldus_v10_NLB.RData")

####### Load the Workspace #######
setwd("/Volumes/bedfordlab/Nocturia/Data_Analysis")
load("latrineNoldus_v10_NLB.RData")

min(allTrials$TimestampNormalRound, na.rm = T) # "2000-01-01 06:02:00 America/Denver"
max(allTrials$TimestampNormalRound, na.rm = T) # "2000-01-03 18:30:00 America/Denver"
length(unique(allTrials$FilePath))
table(allTrials$FilePath)

# Adjust start date for Cage 792 latrine cage recording to align with home cage recording
allTrials2<- allTrials %>% 
  mutate(TimestampNormalRound = case_when(FilePath == "08182022_08.21.40_Pi10_C2_792_latrine" ~ TimestampNormalRound + lubridate::days(1),
                                          .default = TimestampNormalRound))

# Remove NA Timestamps
allTrials2<- allTrials2 %>% drop_na(TimestampNormalRound)

allTrials2<- allTrials2 %>% group_by(FilePath) %>%
  mutate(segment = case_when(
    TimestampNormalRound >= as.POSIXct("2000-01-01 06:00:00 America/Denver") & TimestampNormalRound < as.POSIXct("2000-01-01 18:00:00 America/Denver") ~ "L1",
    TimestampNormalRound >= as.POSIXct("2000-01-01 18:00:00 America/Denver") & TimestampNormalRound < as.POSIXct("2000-01-02 06:00:00 America/Denver") ~ "D1",
    TimestampNormalRound >= as.POSIXct("2000-01-02 06:00:00 America/Denver") & TimestampNormalRound < as.POSIXct("2000-01-02 18:00:00 America/Denver") ~ "L2",
    TimestampNormalRound >= as.POSIXct("2000-01-02 18:00:00 America/Denver") & TimestampNormalRound < as.POSIXct("2000-01-03 06:00:00 America/Denver") ~ "D2",
    TimestampNormalRound >= as.POSIXct("2000-01-03 06:00:00 America/Denver") & TimestampNormalRound < as.POSIXct("2000-01-03 18:00:00 America/Denver") ~ "L3",
    TimestampNormalRound >= as.POSIXct("2000-01-03 18:00:00 America/Denver") & TimestampNormalRound < as.POSIXct("2000-01-04 06:00:00 America/Denver") ~ "D3",
    TimestampNormalRound >= as.POSIXct("2000-01-04 06:00:00 America/Denver") & TimestampNormalRound < as.POSIXct("2000-01-05 18:00:00 America/Denver") ~ "L4",
    .default = NA)) %>% ungroup() 

table(allTrials2$segment)
sum(is.na(allTrials2$segment))

allTrials2<- allTrials2 %>% group_by(FilePath) %>%
  mutate(Time = case_when(segment %in% c("L1", "L2", "L3", "L4") ~ "light",
                          segment %in% c("D1", "D2", "D3") ~ "dark",
                          .default = NA)) %>% ungroup()

table(allTrials2$Time)
sum(is.na(allTrials2$Time))

# Re-order factor levels
str(allTrials2)
allTrials2<- allTrials2 %>%
  mutate(Age = fct_relevel(Age, c("Young", "Old"))) %>%
  mutate(Time = fct_relevel(Time, c("light", "dark"))) %>%
  mutate(CageType = as.factor(CageType))

# Remove incomplete (<12h) or unpaired segments (based on visual inspection of plots in /Users/nicole/Desktop/NocEtho_Outputs/Plots/cumSum_total)
allTrials3<- {allTrials2 %>% 
  mutate(keep = case_when(Trial == "04232024_88" & segment %in% c("L2", "D2", "L3") ~ "remove",
                          Trial == "04252024_86" & segment %in% c("L1", "D3") ~ "remove",
                          Trial == "04252024_87" & segment %in% c("D3") ~ "remove",
                          Trial == "04252024_88" & segment %in% c("D3") ~ "remove",
                          Trial == "04272024_89" & segment %in% c("D3") ~ "remove",
                          Trial == "04272024_90" & segment %in% c("D3") ~ "remove",
                          Trial == "04272024_1106" & segment %in% c("D3") ~ "remove",
                          Trial == "05022024_89" & segment %in% c("D3") ~ "remove",
                          Trial == "05022024_90" & segment %in% c("L3", "D3") ~ "remove",
                          Trial == "05022024_1106" & segment %in% c("D3") ~ "remove",
                          Trial == "05082024_91" & segment %in% c("L1") ~ "remove",
                          Trial == "05082024_92" & segment %in% c("L1") ~ "remove",
                          Trial == "05082024_93" & segment %in% c("L1") ~ "remove",
                          Trial == "05142024_94" & segment %in% c("L1", "L3") ~ "remove",
                          Trial == "05142024_95" & segment %in% c("L1", "L3") ~ "remove",
                          Trial == "05142024_1253" & segment %in% c("L1") ~ "remove",
                          Trial == "05162024_1253" & segment %in% c("D3") ~ "remove",
                          Trial == "05192024_2158" & segment %in% c("D3") ~ "remove",
                          Trial == "05192024_2159" & segment %in% c("D3") ~ "remove",
                          Trial == "05192024_2160" & segment %in% c("D3") ~ "remove",
                          Trial == "05212024_2158" & segment %in% c("D3") ~ "remove",
                          Trial == "05212024_2159" & segment %in% c("D3") ~ "remove",
                          Trial == "05212024_2160" & segment %in% c("D3") ~ "remove",
                          Trial == "05212024_2161" & segment %in% c("D3") ~ "remove",
                          Trial == "08032022_746" & segment %in% c("L3") ~ "remove",
                          Trial == "08032022_747" & segment %in% c("L3") ~ "remove",
                          Trial == "08182022_791" & segment %in% c("L3") ~ "remove",
                          Trial == "08182022_792" & segment %in% c("D1", "L2", "D3", "L4") ~ "remove",
                          Trial == "08202022_791" & segment %in% c("L3") ~ "remove",
                          Trial == "08202022_792" & segment %in% c("L3") ~ "remove",
                          Trial == "09232023_81" & segment %in% c("L3") ~ "remove",
                          Trial == "10042023_82" & segment %in% c("L3") ~ "remove",
                          Trial == "10042023_83" & segment %in% c("D2", "L3") ~ "remove",
                          Trial == "10252023_84" & segment %in% c("L3") ~ "remove",
                          Trial == "10252023_85" & segment %in% c("L3") ~ "remove",
                          Trial == "10272022_875" & segment %in% c("L3") ~ "remove",
                          Trial == "10272022_876" & segment %in% c("L3") ~ "remove",
                          Trial == "10272023_84" & segment %in% c("L2") ~ "remove",
                          Trial == "10272023_85" & segment %in% c("L2") ~ "remove",
                          Trial == "10292022_875" & segment %in% c("L3") ~ "remove",
                          Trial == "10292022_876" & segment %in% c("L3") ~ "remove",
                          .default = "keep")) %>%
  filter(keep == "keep")}

# Make hour and minute time stamps (drop date)
allTrials3<- allTrials3 %>%
  mutate(TimestampNormalRound_Date_hr = as.POSIXct(strftime(TimestampNormalRound, format="%m/%d/%Y %H:%M:%S"), format="%m/%d/%Y %H"), # 48hr, standard date
         TimestampNormalRound_hr = as.POSIXct(strftime(TimestampNormalRound, format="%H"), format="%H"), # 24hr, system date
         TimestampNormalRound_min = as.POSIXct(strftime(TimestampNormalRound, format="%H:%M"), format="%H:%M")) # 24hr, system date

# Convert character columns to factors
allTrials3<- allTrials3 %>% mutate(across(where(is_character), as_factor))
str(allTrials3)

# Make new shifted Timestamp column so plot starts at ZT0
library(lubridate)
allTrials3<- allTrials3 %>%
  mutate(TimestampNormalRound_min_shift = TimestampNormalRound_min - hours(6)) %>%
  mutate(TimestampNormalRound_min_shift = gsub(as.Date(Sys.time()) - 1, as.Date(Sys.time()), TimestampNormalRound_min_shift))

allTrials3$TimestampNormalRound_min_shift = as.POSIXct(allTrials3$TimestampNormalRound_min_shift,
                                                            format = "%Y-%m-%d %H:%M:%S", tz = "America/Denver")

# Add in filter paper data (see Notcuria_FP_v4_NLB.R)
HCLC_filterPapers<- read.csv("/Volumes/bedfordlab/Nocturia/Data_Analysis/HCLC_filterPapers.csv")
HCLC_filterPapers<- HCLC_filterPapers[,-1]
HCLC_filterPapers<- HCLC_filterPapers %>% mutate(across(where(is_character), as_factor))
HCLC_filterPapers<- HCLC_filterPapers %>% mutate(across(where(is_integer), as_factor))
colnames(HCLC_filterPapers)[12]<- "CageType"
str(HCLC_filterPapers)
HCLC_filterPapers<- HCLC_filterPapers %>%
  mutate(Age = fct_relevel(Age, c("Young", "Old")))

FP_info<- HCLC_filterPapers %>%
  select(Round, Age, Sex, Cage, numMice, Trial, segment, Time, trialType) %>%
  filter(!segment == "") %>%
  unique()

allTrials3<- left_join(allTrials3, FP_info)
table(allTrials3$trialType)
sum(is.na(allTrials3$trialType)) # all Cage 746 and 747 Record trials had 48h filter paper
allTrials3<- allTrials3 %>% mutate(trialType = replace_na(trialType, "Record"))
sum(is.na(allTrials3$trialType)) # 0

# Calculate per minute relative latrine cage activity
allTrials3_wide<- allTrials3 %>%
  select(!c(FilePath, Trial.CageType, Cam, min_avg, z_score, keep)) %>%
  pivot_wider(names_from = CageType, values_from = mm_norm) %>%
  drop_na() %>%
  mutate(totAct_norm = home + latrine,
         rel_latAct = latrine / (home + latrine))

# Plot Activity Data (using Record trials only)
ggplot(data = allTrials3_wide %>% filter(trialType == "Record" & Sex == "M" & Cage != "876"),
       aes(x = TimestampNormalRound_min_shift, y = rel_latAct, group = Age, colour = Age, fill = Age)) +
  geom_smooth(method = "gam", span = 1, linewidth = 0.4) +
  scale_colour_manual(values=c("#6c8c24", "#4d5c2c")) +
  scale_fill_manual(values=c("#6c8c24", "#4d5c2c")) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig2.Relative_LatCageActivity_M_24h.pdf", width=8, height=6, units="cm", dpi=1500, useDingbats=FALSE)

ggplot(data = allTrials3_wide %>% filter(trialType == "Record" & Sex == "F"),
       aes(x = TimestampNormalRound_min_shift, y = rel_latAct, group = Age, colour = Age, fill = Age)) +
  geom_smooth(method = "gam", span = 1, linewidth = 0.4) +
  scale_colour_manual(values=c("#8a1e8e", "#642c66")) +
  scale_fill_manual(values=c("#8a1e8e", "#642c66")) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
          strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig2.Relative_LatCageActivity_F_24h.pdf", width=8, height=6, units="cm", dpi=1500, useDingbats=FALSE)

# Plot Home and Latrine Activity for each Age/Sex class and each trial type:
plotdir = "/Volumes/bedfordlab/Nocturia/Data_Analysis/Plots/"
allTrials3<- allTrials3 %>% 
  unite(class_Trial, c("Sex", "Age", "trialType"), remove = FALSE) 

class_Trial<- unique(allTrials3$class_Trial)

for (i in class_Trial){
  df <- allTrials3 %>%
    filter(class_Trial == i)
  
  ggplot(data = df %>% filter(Cage != "876"), 
         aes(x = TimestampNormalRound_min_shift, y = mm_norm, group = CageType, colour = CageType, fill = CageType)) +
           geom_smooth(method = "gam", span = 1, linewidth = 0.4) +
           scale_colour_manual(values=c("#00a08a", "#c6340a")) +
           scale_fill_manual(values=c("#00a08a", "#c6340a")) +
           coord_cartesian(ylim = c(0, 0.15)) +  
           xlab(NULL) + ylab(NULL) +
           theme_classic() +  
           theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
                 strip.background = element_blank(), strip.text.x = element_blank())
  ggsave(paste0(plotdir, "FigS2.HomeLat_Activity_24h_", i, ".pdf"), width=8, height=6, units = "cm", useDingbats = FALSE)
  
}

# Make Raster Plot
library(viridis)
perHour<- allTrials3 %>% group_by(across(c(1:8, 13:17, 21, 25))) %>%
  dplyr::summarise(mm_norm_Hour_AVG = mean(mm_norm, na.rm = T)) %>%
  ungroup()

trials<- unique(perHour$Trial)

for (i in trials){
  df <- perHour %>%
    filter(Trial == i)
  
  ggplot(data = df, (aes(x = TimestampNormalRound_Date_hr, y = CageType))) +
    geom_raster(aes(fill=mm_norm_Hour_AVG)) + scale_fill_viridis(option = "viridis") +
    ggtitle(i) +
    theme_classic()
  ggsave(paste0(plotdir, "ActivityRasters/", i, ".pdf"), width = 20, height = 10, units = "cm", useDingbats = FALSE)
  
}

ggplot(data=perHour %>% filter(Trial == "04272024_89"), aes(x=TimestampNormalRound_Date_hr, y=CageType)) +
  geom_raster(aes(fill=mm_norm_Hour_AVG)) + scale_fill_viridis(option = "viridis") +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank())
ggsave("Plots/Fig2.ExampleRaster_04272024_89.pdf", width=10, height=6, units="cm", dpi=1500, useDingbats=FALSE)

# Model Activity Data
library(nlme)
library(lmerTest)
sum(is.na(allTrials3_wide$rel_latAct)) # These are time stamps where both home and latrine cage activity is zero

# # Use time series analysis: https://rdrr.io/cran/nlme/man/corAR1.html
# # Cage 876 is an outlier (mice sniffing fish-eye camera lens according to DST) --> Remove this cage from analysis
# 
# modelLME_M = lme(rel_latAct ~ Time*Age + TimestampNormalRound_min,
#                random = ~1|Cage,
#                correlation = corAR1(form = ~1|Cage),
#                data = allTrials3_wide %>% filter(Sex == "M" & trialType == "Record" & !is.na(rel_latAct) & Cage != "876"),
#                method = "REML")
# summary(modelLME_M) # Timedark:AgeOld p = 0.0000 ; Age: p = 0.4216
# ACF(modelLME_M)
# 
# modelLME_F = lme(rel_latAct ~ Time*Age + TimestampNormalRound_min,
#                  random = ~1|Cage,
#                  correlation = corAR1(form = ~1|Cage),
#                  data = allTrials3_wide %>% filter(Sex == "F" & trialType == "Record" & !is.na(rel_latAct)),
#                  method = "REML")
# summary(modelLME_F) # Timedark:AgeOld p = 0.0000 ; Age: p = 0.0268
# ACF(modelLME_F)

# Use simpler statistics on averaged data (over each light and dark segment)
allTrials3_wide_sum = allTrials3_wide %>%
  group_by(Cohort, Cage, Round, Age, numMice, Sex, segment, Time, trialType) %>%
  dplyr::summarise(avg_totAct_norm = mean(totAct_norm, na.rm = T),
                   avg_rel_latAct = mean(rel_latAct, na.rm = T),
                   n = n()) %>% ungroup() %>% droplevels()

# How many cages don't have matching light and dark trials? 
pairs = as.data.frame(table(allTrials3_wide_sum$Cage, allTrials3_wide_sum$Time, allTrials3_wide_sum$trialType))
colnames(pairs)[1:3] = c("Cage", "Time", "trialType")
goodPairs = pairs %>% filter(Freq >= 1 & trialType == "Record" & Cage != "876") %>% droplevels() 
goodPairs = as.data.frame(table(goodPairs$Cage)) %>% filter(Freq == 2) %>% droplevels() 
colnames(goodPairs)[1] = "Cage"
goodPairs = goodPairs[["Cage"]]
goodPairs # 26

# Now combine all light segments and all dark segments for plotting
allTrials3_wide_sum2 = allTrials3_wide %>%
  filter(Cage %in% goodPairs) %>%
  group_by(Cohort, Cage, Round, Age, numMice, Sex, Time, trialType) %>%
  dplyr::summarise(avg_totAct_norm = mean(totAct_norm, na.rm = T),
                   avg_rel_latAct = mean(rel_latAct, na.rm = T),
                   n = n()) %>% ungroup() %>% droplevels()

ggplot(data=allTrials3_wide_sum2 %>% filter(Sex == "M" & trialType == "Record" & Cage != "876"), 
       aes(Time, avg_rel_latAct, color = Age)) +
  facet_wrap(~Age) +
  geom_line(aes(group=interaction(Cage))) +
  geom_point(size = 1.25) +
  scale_color_manual(values=c("#6c8c24", "#4d5c2c")) +
  ylim(0, 0.6) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig2.Relative_LatCageActivity_M.pdf", width=6, height=6, units="cm", dpi=1500, useDingbats=FALSE)

# Interaction Model:
summary(lmer(avg_rel_latAct ~ Age*Time + as.integer(numMice) + (1|Cage), data = allTrials3_wide_sum %>% 
               filter(Sex == "M" & trialType == "Record" & Cage != "876"))) # AgeOld:Timedark p = 0.0478 *

# One Model per Age class:
summary(lmer(avg_rel_latAct ~ Time + as.integer(numMice) + (1|Cage), data = allTrials3_wide_sum %>% 
               filter(Sex == "M" & Age == "Young" & trialType == "Record" & Cage != "876"))) # Time: Estimate = 0.21349 ***
summary(lmer(avg_rel_latAct ~ Time + as.integer(numMice) + (1|Cage), data = allTrials3_wide_sum %>% 
               filter(Sex == "M" & Age == "Old" & trialType == "Record" & Cage != "876"))) # Time: Estimate = 0.13118 ***

ggplot(data=allTrials3_wide_sum2 %>% filter(Sex == "F" & trialType == "Record"), 
       aes(Time, avg_rel_latAct, color = Age)) +
  facet_wrap(~Age) +
  geom_line(aes(group=interaction(Cage))) +
  geom_point(size = 1.25) +
  scale_color_manual(values=c("#8a1e8e", "#642c66")) +
  ylim(0, 0.6) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig2.Relative_LatCageActivity_F.pdf", width=6, height=6, units="cm", dpi=1500, useDingbats=FALSE)

# Interaction Model:
summary(lmer(avg_rel_latAct ~ Age*Time + as.integer(numMice) + (1|Cage), data = allTrials3_wide_sum %>% 
               filter(Sex == "F" & trialType == "Record"))) # AgeOld:Timedark p = 0.0470 *

# One Model per Age class:
summary(lmer(avg_rel_latAct ~ Time + as.integer(numMice) + (1|Cage), data = allTrials3_wide_sum %>% 
               filter(Sex == "F" & Age == "Young" & trialType == "Record"))) # Time: Estimate = 0.22362 ***
summary(lmer(avg_rel_latAct ~ Time + as.integer(numMice) + (1|Cage), data = allTrials3_wide_sum %>% 
               filter(Sex == "F" & Age == "Old" & trialType == "Record"))) # Time: Estimate = 0.146283 ***

# What are the effect size differences?
library(Rmisc)

allTrials3_wide_sum3<- allTrials3_wide_sum2 %>% 
  filter(trialType == "Record") %>%
  select(-avg_totAct_norm, -n) %>%
  pivot_wider(names_from = Time, values_from = avg_rel_latAct)

allTrials3_wide_sum3$Diff<- allTrials3_wide_sum3$dark - allTrials3_wide_sum3$light
ES_rel_latAct<- summarySE(data=allTrials3_wide_sum3, measurevar = "Diff", groupvars = c("Sex", "Age"))
round(ES_rel_latAct$Diff, digits = 3)

####### Calculate cumulative activity statistics #######

# CumSum of normalized activity (over the whole 48h trial)
cumulative<- allTrials3 %>% 
  group_by(Trial.CageType) %>%
  dplyr::mutate(cumsumA = cumsum(mm_norm)) %>% ungroup()

trials<- unique(cumulative$Trial)
plotdir = "/Volumes/bedfordlab/Nocturia/Data_Analysis/Plots/"

for (i in trials){
  df <- cumulative %>%
    filter(Trial == i)

  ggplot(data = df, (aes(x = TimestampNormalRound, y = cumsumA, group = CageType, color = segment))) +
    geom_line() +
    labs(y = "cumulative activity (%)") +
    geom_vline(xintercept = as.numeric(as.POSIXct("2000-01-01 06:00:00 America/Denver")), linetype = "dotted") +
    geom_vline(xintercept = as.numeric(as.POSIXct("2000-01-01 18:00:00 America/Denver")), linetype = "dotted") +
    geom_vline(xintercept = as.numeric(as.POSIXct("2000-01-02 06:00:00 America/Denver")), linetype = "dotted") +
    geom_vline(xintercept = as.numeric(as.POSIXct("2000-01-02 18:00:00 America/Denver")), linetype = "dotted") +
    geom_vline(xintercept = as.numeric(as.POSIXct("2000-01-03 06:00:00 America/Denver")), linetype = "dotted") +
    geom_vline(xintercept = as.numeric(as.POSIXct("2000-01-03 18:00:00 America/Denver")), linetype = "dotted") +
    ggtitle(i) +
    theme_classic()
  ggsave(paste0(plotdir, "cumSum_total/", i, ".pdf"), width = 20, height = 15, units = "cm", useDingbats = FALSE)

}

# CumSum of normalized activity (over each 12h light/dark segment)
cumulative.segment = allTrials3 %>% 
  group_by(Trial.CageType, startDate, segment) %>%
  dplyr::mutate(cumsumA.segment = cumsum(mm_norm),
         segment_start = min(TimestampNormalRound, na.rm = T),
         segment_end = max(TimestampNormalRound, na.rm = T)) %>% ungroup() 

for (i in trials){
  df <- cumulative.segment %>%
    filter(Trial == i)

  ggplot(data = df, aes(x = TimestampNormalRound, y = cumsumA.segment, group = CageType, color = CageType)) +
    geom_line() +
    labs(y="cumulative activity (%)") +
    ggtitle(i) +
    theme_classic()
  ggsave(paste0(plotdir, "cumSum_segment/", i, ".pdf"), width = 20, height = 15, units = "cm", useDingbats = FALSE)

}

# Add in Ethovision Trial data
Etho_info<- read.csv("/Volumes/bedfordlab/Nocturia/Data_Analysis/Nocturia Lookup Tables - Ethovision.csv", stringsAsFactors = T)
Etho_info<- Etho_info %>% dplyr::mutate(across(where(is_integer), as_factor))
Etho_info<- Etho_info %>% 
  dplyr::mutate(FilePath = as.factor(str_replace_all(fileName_ANSI, pattern = ".txt", replacement = "")))
str(Etho_info)
cumulative.segment = full_join(cumulative.segment, Etho_info %>% select(-startDate))

# Make a time stamp index so every plot starts at 0
library(data.table)
cumulative.segment = cumulative.segment %>% 
  group_by(Trial.CageType, segment) %>%
  dplyr::mutate(TimestampNormalRound_Index = data.table::rleid(TimestampNormalRound_min)) %>%
  dplyr::mutate(maxIndex = max(TimestampNormalRound_Index)) %>% ungroup()

# Re-order factor levels
str(cumulative.segment)
cumulative.segment<- cumulative.segment %>%
  dplyr::mutate(Age = fct_relevel(Age, c("Young", "Old"))) %>%
  dplyr::mutate(Time = fct_relevel(Time, c("light", "dark")))

ggplot(data = cumulative.segment %>% filter(trialType == "Record" & maxIndex >= 690 & Sex == "M"),
       aes(x = TimestampNormalRound_Index, y = cumsumA.segment, group = cageType, colour = cageType, fill = cageType)) +
  facet_wrap(Age ~ Time) +
  geom_smooth(method = "gam", span = 1, linewidth = 0.4) +
  ggtitle("Male Activity Record Trials")
#ggsave("Plots/HCLC_Activity_RecordTrials_M.pdf", width=14, height=10, units="cm", dpi=1500, useDingbats=FALSE)

ggplot(data = cumulative.segment %>% filter(trialType == "Ctrl" & maxIndex >= 690 & Sex == "M"),
       aes(x = TimestampNormalRound_Index, y = cumsumA.segment, group = cageType, colour = cageType, fill = cageType)) +
  facet_wrap(Age ~ Time) +
  geom_smooth(method = "gam", span = 1, linewidth = 0.4) +
  ggtitle("Male Activity Record Trials")
#ggsave("Plots/HCLC_Activity_ControlTrials_M.pdf", width=14, height=10, units="cm", dpi=1500, useDingbats=FALSE)

ggplot(data = cumulative.segment %>% filter(trialType == "Record" & maxIndex >= 690 & Sex == "F"),
       aes(x = TimestampNormalRound_Index, y = cumsumA.segment, group = cageType, colour = cageType, fill = cageType)) +
  facet_wrap(Age ~ Time) +
  geom_smooth(method = "gam", span = 1, linewidth = 0.4) +
  ggtitle("Female Activity Record Trials")
#ggsave("Plots/HCLC_Activity_RecordTrials_F.pdf", width=14, height=10, units="cm", dpi=1500, useDingbats=FALSE)

ggplot(data = cumulative.segment %>% filter(trialType == "Ctrl" & maxIndex >= 690 & Sex == "F"),
       aes(x = TimestampNormalRound_Index, y = cumsumA.segment, group = cageType, colour = cageType, fill = cageType)) +
  facet_grid(Age ~ Time) +
  geom_smooth(method = "gam", span = 1, linewidth = 0.4) +
  ggtitle("Female Activity Control Trials")
#ggsave("Plots/HCLC_Activity_ControlTrials_F.pdf", width=14, height=10, units="cm", dpi=1500, useDingbats=FALSE)

# Calculate the maximum cumulative activity per light/dark segment
cumulative.segment_max<- cumulative.segment %>% 
  group_by(FilePath, Trial, Trial.CageType, trialType, startDate, Cam, Cohort, Cage, CageType, Age, numMice, Sex, segment, Time, segment_start, segment_end) %>%
  dplyr::summarise(maxcum = max(cumsumA.segment)) %>% ungroup()

cumulative.segment_max_wide<- cumulative.segment_max %>% 
  select(-FilePath, -Trial.CageType, -Cam, -segment_start, -segment_end) %>%
  group_by(Trial, trialType, startDate, Cohort, Cage, Age, numMice, Sex, segment) %>%
  pivot_wider(names_from = CageType, values_from = maxcum) %>%
  ungroup() %>%
  dplyr::mutate(totAct_norm = home + latrine,
         rel_latAct = latrine / (home + latrine))

cumulative.segment_max_wide<- cumulative.segment_max_wide %>% mutate_if(is.character, as.factor)
str(cumulative.segment_max_wide)

# Calculate average cumulative normalized total activity
cumulative.segment_max_wide_sum<- cumulative.segment_max_wide %>%
  group_by(Cage, Age, numMice, Sex, trialType, Time) %>%
  dplyr::summarise(avg_totAct_norm = mean(totAct_norm),
                   avg_rel_latAct = mean(rel_latAct),
                   n = n())

# Sample size:
fo = cumulative.segment_max_wide_sum %>% 
  filter(trialType == "Record") %>%
  filter(Cage %in% goodPairs)

range(fo$n) # 1-4
round(mean(fo$n), digits = 1) # 2.1

ggplot(data=cumulative.segment_max_wide_sum %>% filter(Sex == "M" & trialType == "Record" & Cage %in% goodPairs & Cage != "876"), 
       aes(Time, avg_totAct_norm, color = Age)) +
  facet_wrap(~Age) +
  geom_line(aes(group=interaction(Cage))) +
  geom_point(size = 1.25) +
  scale_color_manual(values=c("#6c8c24", "#4d5c2c")) +
  ylim(0, 170) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig2.TotalActivity_Etho_M.pdf", width=6, height=6, units="cm", dpi=1500, useDingbats=FALSE)

# Interaction Model:
summary(lmer(totAct_norm ~ Age*Time + as.integer(numMice) + (1|Cage), data = cumulative.segment_max_wide %>% filter(Sex == "M" & trialType == "Record" & Cage != "876"))) # AgeOld:Timedark p = 0.00204 **
summary(lmer(rel_latAct ~ Age*Time + as.integer(numMice) + (1|Cage), data = cumulative.segment_max_wide %>% filter(Sex == "M" & trialType == "Record" & Cage != "876"))) # AgeOld:Timedark p = 0.64521

# One Model per Age class:
summary(lmer(totAct_norm ~ Time + as.integer(numMice) + (1|Cage), data = cumulative.segment_max_wide %>% filter(Sex == "M" & Age == "Young" & trialType == "Record" & Cage != "876"))) # Time: Estimate = 77.729 ***
summary(lmer(totAct_norm ~ Time + as.integer(numMice) + (1|Cage), data = cumulative.segment_max_wide %>% filter(Sex == "M" & Age == "Old" & trialType == "Record" & Cage != "876"))) # Time: Estimate = 52.851 ***

ggplot(data=cumulative.segment_max_wide_sum %>% filter(Sex == "F" & trialType == "Record" & Cage %in% goodPairs), 
       aes(Time, avg_totAct_norm, color = Age)) +
  facet_wrap(~Age) +
  geom_line(aes(group=interaction(Cage))) +
  geom_point(size = 1.25) +
  scale_color_manual(values=c("#8a1e8e", "#642c66")) +
  ylim(0, 170) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig2.TotalActivity_Etho_F.pdf", width=6, height=6, units="cm", dpi=1500, useDingbats=FALSE)

# Interaction Model:
summary(lmer(totAct_norm ~ Age*Time + as.integer(numMice) + (1|Cage), data = cumulative.segment_max_wide %>% filter(Sex == "F" & trialType == "Record"))) # AgeOld:Timedark p = 0.0862 .
summary(lmer(rel_latAct ~ Age*Time + as.integer(numMice) + (1|Cage), data = cumulative.segment_max_wide %>% filter(Sex == "F" & trialType == "Record"))) # AgeOld:Timedark p = 0.09293 .

# One Model per Age class:
summary(lmer(totAct_norm ~ Time + as.integer(numMice) + (1|Cage), data = cumulative.segment_max_wide %>% filter(Sex == "F" & Age == "Young" & trialType == "Record"))) # Time: Estimate = 104.083 ***
summary(lmer(totAct_norm ~ Time + as.integer(numMice) + (1|Cage), data = cumulative.segment_max_wide %>% filter(Sex == "F" & Age == "Old" & trialType == "Record"))) # Time: Estimate = 77.965 ***

# Effect size?
cumulative.segment_max_wide_sum2<- cumulative.segment_max_wide_sum %>% 
  filter(trialType == "Record") %>%
  select(-avg_rel_latAct, -n) %>%
  pivot_wider(names_from = Time, values_from = avg_totAct_norm) %>%
  na.omit()
  
# Absolute Difference
cumulative.segment_max_wide_sum2$Diff<- cumulative.segment_max_wide_sum2$dark - cumulative.segment_max_wide_sum2$light
ES_totAct<- summarySE(data=cumulative.segment_max_wide_sum2, measurevar = "Diff", groupvars = c("Sex", "Age"))
round(ES_totAct$Diff, digits = 1)

# Percent Difference
cumulative.segment_max_wide_sum2$Diff_pct<- (cumulative.segment_max_wide_sum2$dark - cumulative.segment_max_wide_sum2$light) / cumulative.segment_max_wide_sum2$light
ES_totAct2<- summarySE(data=cumulative.segment_max_wide_sum2, measurevar = "Diff_pct", groupvars = c("Sex", "Age"))
round(ES_totAct2$Diff_pct, digits = 1)

# Add filter paper data
table(HCLC_filterPapers$trialType) # Only 245 FPs will have corresponding video (others are Acclimation trials)
cumulative.segment_max<- left_join(cumulative.segment_max, HCLC_filterPapers) # Home Cage Activity during Record trials will not have an associated FP

# Add weight data
weight<- read.csv("/Volumes/bedfordlab/Nocturia/Data_Analysis/weight.csv", stringsAsFactors = T)
weight<- weight[,-1]
weight<- weight %>% dplyr::mutate(across(where(is_integer), as_factor))
cumulative.segment_max<- left_join(cumulative.segment_max, weight)

ggplot(cumulative.segment_max %>% filter(trialType == "Record"), 
       aes(x = maxcum, y = pct_cover, colour = Cage)) +
  facet_grid(Age ~ Sex) +
  geom_point(size = 1.25, alpha = 0.5) +
  ggtitle("Latrine Cage Activity and FP") 
#ggsave("LatrineCage_Activity_pctCover_allData.pdf", width=14, height=10, units="cm", dpi=1500, useDingbats=FALSE)

# Does Urine Cover correlate with numMice or weight?
summary(lm(pct_cover ~ as.integer(numMice), data = cumulative.segment_max %>% filter(trialType == "Record"))) # NS
summary(lm(pct_cover ~ totWeight, data = cumulative.segment_max %>% filter(trialType == "Record"))) # p = 6.14e-06 ***
summary(lm(pct_cover ~ avgWeight, data = cumulative.segment_max %>% filter(trialType == "Record"))) # p = 2.45e-08 ***

# Does cumulative activity correlate with numMice? 
summary(lm(maxcum ~ as.integer(numMice), data = cumulative.segment_max %>% filter(trialType == "Record"))) # p = 0.00930 **

# Does Urine Cover correlate with max cumulative activity in the latrine cage?
summary(lm(pct_cover ~ maxcum, data = cumulative.segment_max %>% filter(trialType == "Record" & Sex == "M" & Cage != "876"))) # NS

M1 = lmer(pct_cover ~ maxcum + as.numeric(numMice) + Age + totWeight + (1|Cage), data = cumulative.segment_max %>% filter(trialType == "Record" & Sex == "M")) # maxcum p = 0.0107 *; Age p = 0.1506
summary(M1)
# M2 = lmer(pct_cover ~ maxcum + Age + (1|Cage), data = cumulative.segment_max %>% filter(trialType == "Record" & Sex == "M" & Cage != "876")) # maxcum p = 0.0146 *; Age p = 0.0812 .
# summary(M2)

summary(lm(pct_cover ~ maxcum, data = cumulative.segment_max %>% filter(trialType == "Record" & Sex == "F"))) # NS

F1 = lmer(pct_cover ~ maxcum + as.numeric(numMice) + Age + totWeight + (1|Cage), data = cumulative.segment_max %>% filter(trialType == "Record" & Sex == "F")) # maxcum p = 0.00286 **; Age p = 0.78802
summary(F1)
# F2 = lmer(pct_cover ~ maxcum + Age + (1|Cage), data = cumulative.segment_max %>% filter(trialType == "Record" & Sex == "F")) # maxcum p = 0.00325 **; Age p = 0.02129 *
# summary(F2)

# Extract R-squared value from model
library(MuMIn)
r.squaredGLMM(M1)
# r.squaredGLMM(M2)
r.squaredGLMM(F1)
# r.squaredGLMM(F2)

# Is the Age difference explained by body weight?
summary(lmer(pct_cover ~ maxcum*Age + totWeight + as.numeric(numMice) + (1|Cage), data = cumulative.segment_max %>% filter(trialType == "Record" & Sex == "M" & Cage != "876")))
summary(lmer(pct_cover ~ maxcum*Age + totWeight + as.numeric(numMice) + (1|Cage), data = cumulative.segment_max %>% filter(trialType == "Record" & Sex == "F")))

# Plot latrine cage activity and filter paper
library(broom)

# Apply the linear model to each group
# Male Data
MalePred<- cumulative.segment_max %>%
  filter(trialType == "Record" & Sex == "M" & !is.na(pct_cover)) %>%
  group_by(Age) %>%
  group_modify(~{ 
    model <- lm(pct_cover ~ maxcum, data = .x)      # Fit the model within the group
    predictions <- predict(model, se.fit = TRUE)    # Get predictions and their standard errors within the group
    augmented_data <- augment(model, .x)            # Augment the data with fitted values and standard errors
    augmented_data$se_fit <- predictions$se.fit     # Add the standard error of the fitted values to the augmented data
    augmented_data                                  # Return the augmented data for this group
  }) %>%
  ungroup()

ggplot(MalePred, aes(x = maxcum, y = pct_cover, group = Age, colour = Age)) +
  geom_point(size = 1.25) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_ribbon(aes(ymin = .fitted - se_fit, ymax = .fitted + se_fit), alpha = 0.3) +
  scale_color_manual(values=c("#6c8c24", "#4d5c2c")) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank())
ggsave("Plots/Fig2.cumLatAct_pctCover_M.pdf", width=6, height=6, units="cm", dpi=1500, useDingbats=FALSE)

ggplot(cumulative.segment_max %>% filter(trialType == "Record" & Sex == "M"),
       aes(x = maxcum, y = pct_cover)) +
  geom_point(size = 1.25, aes(colour = Age)) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values=c("#6c8c24", "#4d5c2c")) +
  theme_classic() 

# Apply the linear model to each group
# Female Data
FemalePred<- cumulative.segment_max %>%
  filter(trialType == "Record" & Sex == "F" & !is.na(pct_cover)) %>%
  group_by(Age) %>%
  group_modify(~{ 
    model <- lm(pct_cover ~ maxcum, data = .x)      # Fit the model within the group
    predictions <- predict(model, se.fit = TRUE)    # Get predictions and their standard errors within the group
    augmented_data <- augment(model, .x)            # Augment the data with fitted values and standard errors
    augmented_data$se_fit <- predictions$se.fit     # Add the standard error of the fitted values to the augmented data
    augmented_data                                  # Return the augmented data for this group
  }) %>%
  ungroup()

ggplot(FemalePred, aes(x = maxcum, y = pct_cover, group = Age, colour = Age)) +
  geom_point(size = 1.25) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_ribbon(aes(ymin = .fitted - se_fit, ymax = .fitted + se_fit), alpha = 0.3) +
  scale_color_manual(values=c("#8a1e8e", "#642c66")) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank())
ggsave("Plots/Fig2.cumLatAct_pctCover_F.pdf", width=6, height=6, units="cm", dpi=1500, useDingbats=FALSE)

ggplot(cumulative.segment_max %>% filter(trialType == "Record" & Sex == "F"),
       aes(x = maxcum, y = pct_cover)) +
  geom_point(size = 1.25, aes(colour = Age)) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values=c("#8a1e8e", "#642c66")) +
  theme_classic()

####### Plot Filter Paper Data (Accl + Record Trials only) #######
HCLC_filterPapers<- HCLC_filterPapers %>%
  dplyr::mutate(Age = fct_relevel(Age, c("Young", "Old")))

HCLC_filterPapers_sum<- HCLC_filterPapers %>% 
  filter(trialType %in% c("Accl", "Record") & trialDur == 12) %>%
  group_by(Cage, Age, numMice, Sex, Time) %>%
  dplyr::summarise(avg_pct_cover = mean(pct_cover),
                   n = n()) %>%
  ungroup()

# Sample size:
foo = HCLC_filterPapers_sum %>% 
  group_by(Cage, Age, numMice, Sex) %>%
  dplyr::summarise(n = sum(n))

range(foo$n) # 4-15
round(mean(foo$n), digits = 1) # 8.7

table(HCLC_filterPapers_sum$Cage) # All cages have matching light and dark trials 

ggplot(data=HCLC_filterPapers_sum %>% filter(Sex == "M"), 
       aes(Time, avg_pct_cover, color = Age)) +
  facet_wrap(~Age) +
  geom_line(aes(group=interaction(Cage))) +
  geom_point(size = 1.25) +
  scale_color_manual(values=c("#6c8c24", "#4d5c2c")) +
  ylim(0, 45) +
  theme_classic()
#ggsave("Plots/AVG_FilterPaper_pctCover_Accl_Record_M.pdf", width=6, height=6, units="cm", dpi=1500, useDingbats=FALSE)

ggplot(data=HCLC_filterPapers_sum %>% filter(Sex == "F"), 
       aes(Time, avg_pct_cover, color = Age)) +
  facet_wrap(~Age) +
  geom_line(aes(group=interaction(Cage))) +
  geom_point(size = 1.25) +
  scale_color_manual(values=c("#8a1e8e", "#642c66")) +
  ylim(0, 45) +
  theme_classic()
#ggsave("Plots/AVG_FilterPaper_pctCover_Accl_Record_F.pdf", width=6, height=6, units="cm", dpi=1500, useDingbats=FALSE)

# Replace missing value with average totWeight per class:
weight_sum<- summarySE(weight, measurevar="totWeight", groupvars=c("Age", "Sex"), na.rm=T)

# Merge the 'weight_sum' values into 'weight' using Age and Sex
weight<- weight %>%
  left_join(weight_sum %>%
              select(Age, Sex, totWeight), by = c("Age", "Sex")) %>%  # Fill missing totWeight based on the joined totWeight from weight_sum
  mutate(totWeight = coalesce(totWeight.x, totWeight.y)) %>% # coalesce will use values from weight_sum if totWeight is NA
  select(-totWeight.x, -totWeight.y)  # Remove the temporary columns after merge

# Plot total Urine Cover per treatment group (standardized by totWeight)
HCLC_filterPapers<- left_join(HCLC_filterPapers, weight)

# Does total urine cover differ between adult and aged mice?
summary(lm(pct_cover ~ Age, data = HCLC_filterPapers %>% filter(trialType %in% c("Accl", "Record") & trialDur == 12 & Sex == "F"))) # p = 0.000514 ***
summary(lm(pct_cover ~ Age, data = HCLC_filterPapers %>% filter(trialType %in% c("Accl", "Record") & trialDur == 12 & Sex == "M"))) # p = 0.0156 *

# Need to "control" for numMice and totWeight?
summary(lm(pct_cover ~ as.numeric(numMice), data = HCLC_filterPapers %>% filter(trialType %in% c("Accl", "Record") & trialDur == 12))) # p = 0.0179 *
summary(lm(pct_cover ~ totWeight, data = HCLC_filterPapers %>% filter(trialType %in% c("Accl", "Record") & trialDur == 12))) # p = 4.62e-12 ***

# Is there a linear relationship?
ggplot(HCLC_filterPapers %>% filter(trialType %in% c("Accl", "Record") & trialDur == 12),
       aes(x = totWeight, y = pct_cover)) +
  geom_point(size = 1.25, aes(colour = Age)) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic()

summary(lm(pct_cover ~ totWeight, data = HCLC_filterPapers %>% filter(trialType %in% c("Accl", "Record") & trialDur == 12))) # Adjusted R-squared:  0.1818 ***

# No effect of Age:
summary(lmer(pct_cover ~ Age + Time + totWeight + (1|Cage), data = HCLC_filterPapers %>% filter(trialType %in% c("Accl", "Record") & trialDur == 12 & Sex == "F"))) # p = 0.173
summary(lmer(pct_cover ~ Age + (1|Cage), data = HCLC_filterPapers %>% filter(trialType %in% c("Accl", "Record") & trialDur == 12 & Sex == "F"))) # p = 0.0731

summary(lmer(pct_cover ~ Age + Time + totWeight + (1|Cage), data = HCLC_filterPapers %>% filter(trialType %in% c("Accl", "Record") & trialDur == 12 & Sex == "M"))) # p = 0.9549
summary(lmer(pct_cover ~ Age + (1|Cage), data = HCLC_filterPapers %>% filter(trialType %in% c("Accl", "Record") & trialDur == 12 & Sex == "M"))) # p = 0.34173

HCLC_filterPapers_sum2<- left_join(HCLC_filterPapers_sum, weight)
HCLC_filterPapers_sum2<- HCLC_filterPapers_sum2 %>% group_by(Cage, Age, numMice, Sex, totWeight) %>%
  dplyr:: summarise(pct_cover_LD_avg = mean(avg_pct_cover),
                    n = n())

summaryF <- summarySE(data = HCLC_filterPapers_sum2 %>% filter(Sex == "F"),
                      measurevar = "pct_cover_LD_avg", groupvars = "Age", na.rm=T)

ggplot(data=HCLC_filterPapers_sum2 %>% filter(Sex == "F"), aes(x=Age, y=pct_cover_LD_avg)) +
  geom_bar(data=summaryF, aes(fill=Age), stat="identity", width=0.85) +
  geom_errorbar(data=summaryF, aes(ymin=pct_cover_LD_avg-se, ymax=pct_cover_LD_avg+se, x=Age), size=0.5, width=0.25) + 
  geom_point(aes(fill=Age), colour="black", pch=21, position=position_jitter(0.1), size=1) +
  scale_fill_manual(values=c("#8a1e8e", "#642c66")) +
  xlab(NULL) + ylab(NULL) + 
  theme_classic() +
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank())
ggsave("Plots/Fig2.pct_cover_LD_avg_Record_Accl_F.pdf", width=3, height=6, units="cm", dpi=1500, useDingbats=FALSE)

summaryM <- summarySE(data = HCLC_filterPapers_sum2 %>% filter(Sex == "M"),
                      measurevar = "pct_cover_LD_avg", groupvars = "Age", na.rm=T)

ggplot(data=HCLC_filterPapers_sum2 %>% filter(Sex == "M"), aes(x=Age, y=pct_cover_LD_avg)) +
  geom_bar(data=summaryM, aes(fill=Age), stat="identity", width=0.85) +
  geom_errorbar(data=summaryM, aes(ymin=pct_cover_LD_avg-se, ymax=pct_cover_LD_avg+se, x=Age), size=0.5, width=0.25) + 
  geom_point(aes(fill=Age), colour="black", pch=21, position=position_jitter(0.1), size=1) +
  scale_fill_manual(values=c("#6c8c24", "#4d5c2c")) +
  xlab(NULL) + ylab(NULL) + 
  theme_classic() +
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank())
ggsave("Plots/Fig2.pct_cover_LD_avg_Record_Accl_M.pdf", width=3, height=6, units="cm", dpi=1500, useDingbats=FALSE)

## Light Phase ONLY plots:
summaryF <- summarySE(data = HCLC_filterPapers_sum %>% filter(Sex == "F" & Time == "light"),
                      measurevar = "avg_pct_cover", groupvars = "Age", na.rm=T)

ggplot(data=HCLC_filterPapers_sum %>% filter(Sex == "F" & Time == "light"), aes(x=Age, y=avg_pct_cover)) +
  geom_bar(data=summaryF, aes(fill=Age), stat="identity", width=0.85) +
  geom_errorbar(data=summaryF, aes(ymin=avg_pct_cover-se, ymax=avg_pct_cover+se, x=Age), size=0.5, width=0.25) + 
  geom_point(aes(fill=Age), colour="black", pch=21, position=position_jitter(0.1), size=1) +
  scale_fill_manual(values=c("#8a1e8e", "#642c66")) +
  xlab(NULL) + ylab(NULL) + 
  theme_classic() +
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank())
ggsave("Plots/Fig2.avg_pct_cover_LIGHT_ONLY_Record_Accl_F.pdf", width=3, height=6, units="cm", dpi=1500, useDingbats=FALSE)

summary(lmer(pct_cover ~ Age + (1|Cage), data = HCLC_filterPapers %>% filter(trialType %in% c("Accl", "Record") & trialDur == 12 & Sex == "F" & Time == "light"))) # p = 0.173

summaryM <- summarySE(data = HCLC_filterPapers_sum %>% filter(Sex == "M" & Time == "light"),
                      measurevar = "avg_pct_cover", groupvars = "Age", na.rm=T)

ggplot(data=HCLC_filterPapers_sum %>% filter(Sex == "M" & Time == "light"), aes(x=Age, y=avg_pct_cover)) +
  geom_bar(data=summaryM, aes(fill=Age), stat="identity", width=0.85) +
  geom_errorbar(data=summaryM, aes(ymin=avg_pct_cover-se, ymax=avg_pct_cover+se, x=Age), size=0.5, width=0.25) + 
  geom_point(aes(fill=Age), colour="black", pch=21, position=position_jitter(0.1), size=1) +
  scale_fill_manual(values=c("#6c8c24", "#4d5c2c")) +
  xlab(NULL) + ylab(NULL) + 
  theme_classic() +
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank())
ggsave("Plots/Fig2.avg_pct_cover_LIGHT_ONLY_Record_Accl_M.pdf", width=3, height=6, units="cm", dpi=1500, useDingbats=FALSE)

summary(lmer(pct_cover ~ Age + (1|Cage), data = HCLC_filterPapers %>% filter(trialType %in% c("Accl", "Record") & trialDur == 12 & Sex == "M" & Time == "light"))) # p = 0.4479

## Weight Plots:
lut<- read.csv("/Volumes/bedfordlab/Nocturia/Data_Analysis/Nocturia Lookup Tables - LUT.csv", stringsAsFactors = T)
summary(lm(Weight ~ Age, data = lut %>% filter(Sex == "F"))) # p = 5.87e-12 ***
summary(lm(Weight ~ Age, data = lut %>% filter(Sex == "M"))) # p = 2.28e-06 ***

weight<- weight %>%
  dplyr::mutate(Age = fct_relevel(Age, c("Young", "Old")))

summaryF <- summarySE(data = weight %>% filter(Sex == "F"),
                      measurevar = "totWeight", groupvars = "Age", na.rm=T)

ggplot(data=weight %>% filter(Sex == "F"), aes(x=Age, y=totWeight)) +
  geom_bar(data=summaryF, aes(fill=Age), stat="identity", width=0.85) +
  geom_errorbar(data=summaryF, aes(ymin=totWeight-se, ymax=totWeight+se, x=Age), size=0.5, width=0.25) + 
  geom_point(aes(fill=Age), colour="black", pch=21, position=position_jitter(0.1), size=1) +
  scale_fill_manual(values=c("#8a1e8e", "#642c66")) +
  theme_classic() +
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank())
#ggsave("Plots/totWeight_perCage_F.pdf", width=3, height=6, units="cm", dpi=1500, useDingbats=FALSE)

summary(lm(totWeight ~ Age, data = weight %>% filter(Sex == "F"))) # NS

summaryM <- summarySE(data = weight %>% filter(Sex == "M"),
                      measurevar = "totWeight", groupvars = "Age", na.rm=T)

ggplot(data=weight %>% filter(Sex == "M"), aes(x=Age, y=totWeight)) +
  geom_bar(data=summaryM, aes(fill=Age), stat="identity", width=0.85) +
  geom_errorbar(data=summaryM, aes(ymin=totWeight-se, ymax=totWeight+se, x=Age), size=0.5, width=0.25) + 
  geom_point(aes(fill=Age), colour="black", pch=21, position=position_jitter(0.1), size=1) +
  scale_fill_manual(values=c("#6c8c24", "#4d5c2c")) +
  theme_classic() +
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank())
#ggsave("Plots/totWeight_perCage_M.pdf", width=3, height=6, units="cm", dpi=1500, useDingbats=FALSE)

summary(lm(totWeight ~ Age, data = weight %>% filter(Sex == "M"))) # NS
