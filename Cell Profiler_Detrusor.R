# Libraries
library(Rmisc)
library(dplyr)
library(data.table)
library(tidyverse)
library(stringr)
library(magrittr)
library(ggplot2)
library(ggdark)
library(ggbeeswarm)
library(plotly)
library(gridExtra)
library(lme4)
library(lmerTest)
library(wesanderson)
#rbind
library(plyr)
library(Matrix)
library(data.table)
library(ggpubr)
library(lme4)
library(emmeans)



#### DETRUSOR ####

# Set working directory and plot directory
setwd("/Volumes/bedfordlab/Nocturia/Bladder Images/CellProfiler/OUTPUT/Detrusor")
getwd()
plotdir = paste(getwd(),"/Plots/", sep = "") 


### Import Data ###
expanded1 <- read.csv("nuclei_expanded1.csv")   
expanded2 <- read.csv("nuclei_expanded2.csv")  
expanded3 <- read.csv("nuclei_expanded3.csv")

#removed unwanted columns (should have 23 variables after this)
expanded1 <- expanded1[ -c(1:4,6:7,13:21,23:24) ]
expanded2 <- expanded2[ -c(1:4,6:7,13:21,23:24) ]
expanded3 <- expanded3[ -c(1:4,6:7,13:21,23:24) ]


# Extract the mouse ID and place it into the MouseID column to fix the mouseID 
expanded3new <- expanded3 %>%
  mutate(Metadata_Mouse = ifelse(is.na(Metadata_Mouse) | Metadata_Mouse == "", 
                                 sub("^[^_]+_[^_]+_[^_]+_[^_]+_[^_]+_([^_]+)_.*$", "\\1", FileName_Arntl), 
                                 MouseID))

#bind the two data sheets together
Expanded_det <- rbind(expanded1, expanded2, expanded3new)


### Plot count distributions:

##Piezo
#with zeros
mean(Expanded_det$Children_Piezo1_puncta_Count) #0.6432272
range(Expanded_det$Children_Piezo1_puncta_Count) #0 21

#without zeros
Expanded_det_noZero_Piezo = Expanded_det %>% filter(Children_Piezo1_puncta_Count > 0)
hist(Expanded_det_noZero_Piezo$Children_Piezo1_puncta_Count, breaks = 100)
mean(Expanded_det_noZero_Piezo$Children_Piezo1_puncta_Count) #2.142602
range(Expanded_det_noZero_Piezo$Children_Piezo1_puncta_Count) #1 21


##Arntl
#with zeros
mean(Expanded_det$Children_Arntl_puncta_Count) #0.4031557
range(Expanded_det$Children_Arntl_puncta_Count) #0 33

#without zeros
Expanded_det_noZero_Arntl = Expanded_det %>% filter(Children_Arntl_puncta_Count > 0)
hist(Expanded_det_noZero_Arntl$Children_Arntl_puncta_Count, breaks = 100)
mean(Expanded_det_noZero_Arntl$Children_Arntl_puncta_Count) #2.336439
range(Expanded_det_noZero_Arntl$Children_Arntl_puncta_Count) #1 33

##Per1
#with zeros
mean(Expanded_det$Children_Per1_puncta_Count) #0.5095564
range(Expanded_det$Children_Per1_puncta_Count) #0 34

#without zeros
Expanded_det_noZero_Per1 = Expanded_det %>% filter(Children_Per1_puncta_Count > 0)
hist(Expanded_det_noZero_Per1$Children_Per1_puncta_Count, breaks = 100)
mean(Expanded_det_noZero_Per1$Children_Per1_puncta_Count) #2.647896
range(Expanded_det_noZero_Per1$Children_Per1_puncta_Count) #1 34



####Data Wrangleing####
# Create the first data frame with the positive values for each count
# Turns the puncta counts into binary. If punct was dected within a cell the value is 1, if not detected the value is zero
All_Data <- Expanded_det %>%
  mutate(
    Per1_pos = case_when(
      Children_Per1_puncta_Count == 0 ~ 0,
      Children_Per1_puncta_Count > 0 ~ 1,
      .default = NA
    ),
    Arntl_pos = case_when(
      Children_Arntl_puncta_Count == 0 ~ 0,
      Children_Arntl_puncta_Count > 0 ~ 1,
      .default = NA
    ),
    Piezo1_pos = case_when(
      Children_Piezo1_puncta_Count == 0 ~ 0,
      Children_Piezo1_puncta_Count > 0 ~ 1,
      .default = NA
    )
  )


# Summarize the data by FileName_DAPI to show the amount of positive cells per image
All_Data2 <- All_Data %>%
  group_by(FileName_DAPI, Metadata_Date, Metadata_Mouse, Metadata_Time) %>%
  summarise(
    numPos_Piezo1 = sum(Piezo1_pos, na.rm = TRUE),
    numPos_Arntl = sum(Arntl_pos, na.rm = TRUE),
    numPos_Per1 = sum(Per1_pos, na.rm = TRUE),
    numDAPI = n()
  )


#make the puncta data into percentages
All_Data3 = All_Data2 %>%
  mutate(pct_Peizo1 = numPos_Piezo1 / numDAPI * 100) %>%
  mutate(pct_Arntl = numPos_Arntl / numDAPI * 100) %>%
  mutate(pct_Per1 = numPos_Per1 / numDAPI * 100)


####Plotting####

#assign sex to mice
All_Data3$Sex <- ifelse(All_Data3$Metadata_Mouse %in% c("81.02", "81", "84.04", "85.06", "85.08", "84.01", "84.03", "85.07","81.00", "8103", "8101"), "Female",
                        ifelse(All_Data3$Metadata_Mouse %in% c("791.384", "Ctrl.1", "Ctrl.2", "83.02", "Ctrl.3", "83.01", "8602", "243304", "ctrl6", "8704"), "Male", NA))
#assign age to mice
All_Data3$Age <- ifelse(All_Data3$Metadata_Mouse %in% c("791.384", "Ctrl.1", "Ctrl.2", "83.02", "Ctrl.3", "83.01", "81", "81.02", "81.00", "8103", "8101"), "Old",
                        ifelse(All_Data3$Metadata_Mouse %in% c("84.04", "85.06", "85.08", "84.01", "84.03", "85.07", "8602", "243304", "ctrl6", "8704"), "Young", NA))

#fix the mouse names to all have the same format
All_Data3$Metadata_Mouse[All_Data3$Metadata_Mouse == 8101] <- 81.01
All_Data3$Metadata_Mouse[All_Data3$Metadata_Mouse == 81] <- 81.01
All_Data3$Metadata_Mouse[All_Data3$Metadata_Mouse == 8103] <- 81.03
All_Data3$Metadata_Mouse[All_Data3$Metadata_Mouse == 243304] <- 2433.04
All_Data3$Metadata_Mouse[All_Data3$Metadata_Mouse == 8602] <- 86.02
All_Data3$Metadata_Mouse[All_Data3$Metadata_Mouse == 8704] <- 87.04
All_Data3$Metadata_Mouse[All_Data3$Metadata_Mouse == "ctrl6"] <- "Ctrl.6"

#Creating a "Sex" column 
female_data <- All_Data3[All_Data3$Sex == "Female", ]
male_data <- All_Data3[All_Data3$Sex == "Male", ]
female_data$Sex <- "Female"
male_data$Sex <- "Male"

#create new data frame
combined_data <- rbind(female_data, male_data)

# Create a new column that combines Age and Sex for grouping
combined_data$Group <- paste(combined_data$Age, combined_data$Sex)


### Piezo ###
#ALL
ggplot(combined_data, aes(x = Group, y = pct_Peizo1, fill = Sex)) +
  geom_boxplot() +
  labs(title = "Percentage of Postive Piezo Puncta by Age and Sex in Mice",
       x = "Group",
       y = "Piezo Percentage") +
  theme_minimal() +
  scale_fill_manual(values = c("Female" = "salmon", "Male" = "skyblue"))+
  theme(legend.position = "none")
ggsave("/Volumes/bedfordlab/Nocturia/Data_Analysis/Cell Profiler/Piezo_boxplot.pdf", width=12, height=6, units="cm", dpi=1500, useDingbats=FALSE)


# Convert Group to a factor
combined_data$Group <- factor(combined_data$Group, levels = c("Young Female", "Old Female", "Young Male", "Old Male"))

# Run ANOVA
anova_result_peizo <- aov(pct_Peizo1 ~ Group, data = combined_data)
summary(anova_result_peizo)

# Tukey's HSD for pairwise comparisons
tukey_result_peizo <- TukeyHSD(anova_result_peizo)
tukey_result_peizo

#sex differeneces in young are significant but as the mice age the sex differences become insignifcant. There is also signifcance betwen O and Y Male and O and Y females


#Female
ggplot(female_data, aes(x = Metadata_Time, y = pct_Peizo1, fill = Metadata_Time)) +   
  geom_bar(stat = "identity", position = "dodge") +   
  geom_point(aes(color = Metadata_Mouse), size = 3) +  # Add black points with size adjustment  
  labs(
    title = "Female",
    x = "Time (AM/PM)",
    y = "Piezo Percentage"
  ) +   
  theme_minimal() +   
  scale_fill_manual(values = c("AM" = "darkgoldenrod1", "PM" = "lightblue")) +   
  scale_x_discrete(labels = c("AM", "PM")) +   
  theme(
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +   
  facet_grid(. ~ Age, scales = "free", space = "free")
ggsave("/Volumes/bedfordlab/Nocturia/Data_Analysis/Cell Profiler/Piezo_F_barplot.pdf", width=12, height=6, units="cm", dpi=1500, useDingbats=FALSE)

# Run ANOVA
anova_result_peizo_F <- aov(pct_Peizo1 ~ Metadata_Time * Age, data = female_data)
summary(anova_result_peizo_F)

# Tukey's HSD for pairwise comparisons
tukey_result_peizo_F <- TukeyHSD(anova_result_peizo_F)
tukey_result_peizo_F

#pm old to am old sing, young am to old am NOT significant, pm young to pm old sign, 

#Male
ggplot(male_data, aes(x = Metadata_Time, y = pct_Peizo1, fill = Metadata_Time)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_point(aes(color = Metadata_Mouse), size = 3) +  # Add black points
  labs(
    title = "Male",
    x = "Time (AM/PM)",
    y = "Piezo Percentage"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("AM" = "darkgoldenrod1", "PM" = "lightblue")) +
  scale_x_discrete(labels = c("AM", "PM")) +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +   
  facet_grid(. ~ Age, scales = "free", space = "free")

ggsave("/Volumes/bedfordlab/Nocturia/Data_Analysis/Cell Profiler/Piezo_M_barplot.pdf", width=12, height=6, units="cm", dpi=1500, useDingbats=FALSE)


# Run ANOVA
anova_result_peizo_M <- aov(pct_Peizo1 ~ Metadata_Time * Age, data = male_data)
summary(anova_result_peizo_M)

# Tukey's HSD for pairwise comparisons
tukey_result_peizo_M <- TukeyHSD(anova_result_peizo_M)
tukey_result_peizo_M




### ARNTL ###
#ALL
ggplot(combined_data, aes(x = Group, y = pct_Arntl, fill = Sex)) +
  geom_boxplot() +
  labs(title = "Percentage of Positive Arntl Puncta by Age and Sex in Mice",
       x = "Group",
       y = "Arntl Percentage") +
  theme_minimal() +
  scale_fill_manual(values = c("Female" = "salmon", "Male" = "skyblue"))+
  theme(legend.position = "none")
ggsave("/Volumes/bedfordlab/Nocturia/Data_Analysis/Cell Profiler/Arntl_boxplot.pdf", width=12, height=6, units="cm", dpi=1500, useDingbats=FALSE)


# Run ANOVA
anova_result_arntl <- aov(pct_Arntl ~ Group, data = combined_data)
summary(anova_result_arntl)

# Tukey's HSD for pairwise comparisons
tukey_result_arntl <- TukeyHSD(anova_result_arntl)
tukey_result_arntl

#sex differeneces in young are significant but as the mice age the sex differences become insignifcant. There is also signifcance betwen O and Y Male and O and Y females

#Female
ggplot(female_data, aes(x = Metadata_Time, y = pct_Arntl, fill = Metadata_Time)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_point(aes(color = Metadata_Mouse), size = 3) +  # Add black points
  labs(
    title = "Female",
    x = "Time (AM/PM)",
    y = "Arntl Percentage"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("AM" = "darkgoldenrod1", "PM" = "lightblue")) +
  scale_x_discrete(labels = c("AM", "PM")) +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +   
  facet_grid(. ~ Age, scales = "free", space = "free")
ggsave("/Volumes/bedfordlab/Nocturia/Data_Analysis/Cell Profiler/Arntl_F_barplot.pdf", width=12, height=6, units="cm", dpi=1500, useDingbats=FALSE)



# Run ANOVA
anova_result_arntl_F <- aov(pct_Arntl ~ Metadata_Time * Age, data = female_data)
summary(anova_result_arntl_F)

# Tukey's HSD for pairwise comparisons
tukey_result_arntl_F <- TukeyHSD(anova_result_arntl_F)
tukey_result_arntl_F


#Male
ggplot(male_data, aes(x = Metadata_Time, y = pct_Arntl, fill = Metadata_Time)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_point(aes(color = Metadata_Mouse), size = 3) +  # Add black points
  labs(
    title = "Male",
    x = "Time (AM/PM)",
    y = "Arntl Percentage"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("AM" = "darkgoldenrod1", "PM" = "lightblue")) +
  scale_x_discrete(labels = c("AM", "PM")) +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +   
  facet_grid(. ~ Age, scales = "free", space = "free")
ggsave("/Volumes/bedfordlab/Nocturia/Data_Analysis/Cell Profiler/Arntl_M_barplot.pdf", width=12, height=6, units="cm", dpi=1500, useDingbats=FALSE)


# Run ANOVA
anova_result_arntl_M <- aov(pct_Arntl ~ Metadata_Time * Age, data = male_data)
summary(anova_result_arntl_M)

# Tukey's HSD for pairwise comparisons
tukey_result_arntl_M <- TukeyHSD(anova_result_arntl_M)
tukey_result_arntl_M




### Per1 ###
ggplot(combined_data, aes(x = Group, y = pct_Per1, fill = Sex)) +
  geom_boxplot() +
  labs(title = "Percentage of Positive Per1 Puncta by Age and Sex in Mice",
       x = "Group",
       y = "Arntl Percentage") +
  theme_minimal() +
  scale_fill_manual(values = c("Female" = "salmon", "Male" = "skyblue"))+
  theme(legend.position = "none") 
ggsave("/Volumes/bedfordlab/Nocturia/Data_Analysis/Cell Profiler/Per1_boxplot.pdf", width=18, height=10, units="cm", dpi=1500, useDingbats=FALSE)


# Run ANOVA
anova_result_per1 <- aov(pct_Per1 ~ Group, data = combined_data)
summary(anova_result_per1)

# Tukey's HSD for pairwise comparisons
tukey_result_per1 <- TukeyHSD(anova_result_per1)
tukey_result_per1

# OF YF is significant and so is young male young female



#Female
ggplot(female_data, aes(x = Metadata_Time, y = pct_Per1, fill = Metadata_Time)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_point(aes(color = Metadata_Mouse), size = 3) +  # Add black points
  labs(
    title = "Female",
    x = "Time (AM/PM)",
    y = "Per1 Percentage"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("AM" = "darkgoldenrod1", "PM" = "lightblue")) +
  scale_x_discrete(labels = c("AM", "PM")) +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +   
  facet_grid(. ~ Age, scales = "free", space = "free")
ggsave("/Volumes/bedfordlab/Nocturia/Data_Analysis/Cell Profiler/Per1_F_barplot.pdf", width=12, height=6, units="cm", dpi=1500, useDingbats=FALSE)


# Run ANOVA
anova_result_per1_F <- aov(pct_Per1 ~ Metadata_Time * Age, data = female_data)
summary(anova_result_per1_F)

# Tukey's HSD for pairwise comparisons
tukey_result_per1_F <- TukeyHSD(anova_result_per1_F)
tukey_result_per1_F



#Male
ggplot(male_data, aes(x = Metadata_Time, y = pct_Per1, fill = Metadata_Time)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_point(aes(color = Metadata_Mouse), size = 3) +  # Add black points
  labs(
    title = "Male",
    x = "Time (AM/PM)",
    y = "Per1 Percentage"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("AM" = "darkgoldenrod1", "PM" = "lightblue")) +
  scale_x_discrete(labels = c("AM", "PM")) +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +   
  facet_grid(. ~ Age, scales = "free", space = "free")
ggsave("/Volumes/bedfordlab/Nocturia/Data_Analysis/Cell Profiler/Per1_M_barplot.pdf", width=12, height=6, units="cm", dpi=1500, useDingbats=FALSE)


# Run ANOVA
anova_result_per1_M <- aov(pct_Per1 ~ Metadata_Time * Age, data = male_data)
summary(anova_result_per1_M)

# Tukey's HSD for pairwise comparisons
tukey_result_per1_M <- TukeyHSD(anova_result_per1_M)
tukey_result_per1_M