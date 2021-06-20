#LOAD LIBRARIES ----
library(tidyr)
library(dplyr)
library(tidyverse)
library(readr)
library(AICcmodavg)
library(ggplot2)
library(ggrepel)
library(polynom)
library(RColorBrewer)
library(ggthemes)
library(geomorph)
library(broom)
library(gginnards)
library(ggfortify)
library(rphylopic)
library(png)
library(gridExtra)


#DATA IMPORT ----
#Import data
ossification_seq <- read_csv("Data/earbones_ossification_events.csv")
ossification_seq$state <-  as.character(ossification_seq$state)

measuraments <- read_csv("Data/measuraments.csv")

growth_curve <- read_csv("Data/growth_data.csv")

#Make natural log relevant columns
measuraments <- measuraments %>% mutate(BZW_log = log(BZW), bullaL_log  = log(bullaL), bullaW_log  = log(bullaW), 
                                        perioticL_log = log(perioticL), perioticW_log = log(perioticW))

ossification_seq <- ossification_seq %>% mutate(TL_mm_log = log(TL_mm))

growth_curve <- growth_curve %>% mutate(TL_mm_log = log(TL_mm))

#Split groups growth data
growth_curve_groups <- growth_curve %>% group_by(group) %>% group_split()
View(growth_curve_groups)

growth_curve_mysticeti <- growth_curve_groups[[1]]
growth_curve_odontoceti <- growth_curve_groups[[2]]

#Split taxa growth data
growth_curve_taxa <- growth_curve %>% group_by(taxon) %>% group_split()
View(growth_curve_taxa)

growth_curve_B.bonaerensis <- growth_curve_taxa[[1]]
growth_curve_B.physalus <- growth_curve_taxa[[2]]
growth_curve_Ph.phocoena <- growth_curve_taxa[[3]]
growth_curve_St.attenuata <- growth_curve_taxa[[4]]

#Split taxa growth sequence stages
ossification_seq_taxa <- ossification_seq %>% group_by(taxon) %>% group_split()
View(ossification_seq_taxa)

ossification_seq_B.bonaerensis <- ossification_seq_taxa[[1]]
ossification_seq_B.physalus <- ossification_seq_taxa[[2]]
ossification_seq_Ph.phocoena <- ossification_seq_taxa[[3]]
ossification_seq_St.attenuata <- ossification_seq_taxa[[4]]

#Split measurements dataset
#Create 2 datasets, 1 for periotic and 1 for bulla
bulla_meas <- measuraments %>% select(2:8, 11:13)
periotic_meas  <- measuraments %>% drop_na() %>% select(2:6, 9:11, 14:15) #getting rid of specimens that do not have the periotic

#Divide measurement data by group
bulla_meas_groups <- bulla_meas %>% group_by(group) %>% group_split()

bulla_meas_mysticeti <- bulla_meas_groups[[1]]
bulla_meas_odontoceti <- bulla_meas_groups[[2]]

periotic_meas_groups <- periotic_meas %>% group_by(group) %>% group_split()

periotic_meas_mysticeti <- periotic_meas_groups[[1]]
periotic_meas_odontoceti <- periotic_meas_groups[[2]]

##Divide measurement data by species
#Bulla
bulla_meas_species <- bulla_meas %>% group_by(taxon) %>% group_split()
View(bulla_meas_species)

bulla_meas_B.bonaerensis <- bulla_meas_species[[1]]
bulla_meas_B.acutorostrata <- bulla_meas_species[[2]]
bulla_meas_B.physalus <- bulla_meas_species[[5]]
bulla_meas_Ph.phocoena <- bulla_meas_species[[17]]
bulla_meas_St.attenuata <- bulla_meas_species[[20]]
#Combine 2 minke species
bulla_meas_minke <- bind_rows(bulla_meas_B.bonaerensis, bulla_meas_B.acutorostrata)
#Replace species to match growth data (all B.bonaerensis)
bulla_meas_minke[bulla_meas_minke == "B.acutorostrata"] <- "B.bonaerensis"
bulla_meas_minke

#Make new dataframe with selected better sampled species for taxon-level analysis
bulla_meas_taxa <- bind_rows(bulla_meas_minke, bulla_meas_B.physalus, bulla_meas_Ph.phocoena, bulla_meas_St.attenuata)
bulla_meas_taxa 

#Periotic
periotic_meas_species <- periotic_meas %>% group_by(taxon) %>% group_split()
View(periotic_meas_species)

periotic_meas_B.bonaerensis <- periotic_meas_species[[1]]
periotic_meas_B.acutorostrata <- periotic_meas_species[[2]]
periotic_meas_Ph.phocoena <- periotic_meas_species[[16]]
periotic_meas_St.attenuata <- periotic_meas_species[[17]]
#Combine 2 minke species
periotic_meas_minke <- bind_rows(periotic_meas_B.bonaerensis, periotic_meas_B.acutorostrata)
#Replace species to match growth data (all B.bonaerensis)
periotic_meas_minke[periotic_meas_minke == "B.acutorostrata"] <- "B.bonaerensis"
periotic_meas_minke

#Make new dataframe with selected better sampled species for taxon-level analysis
periotic_meas_taxa <- bind_rows(periotic_meas_minke, periotic_meas_Ph.phocoena, periotic_meas_St.attenuata)
periotic_meas_taxa 

##Calculate % of growth
#Get age at birth of all taxa
growth_max_age <- growth_curve %>% group_by(taxon) %>% summarize(max_age = max(Age_months))
growth_max_age

#Make new column with % growth for each specimen in each species
#Growth data
growth_curve_B.bonaerensis <- growth_curve_B.bonaerensis %>% mutate(growth_max_age[1,2]) %>% mutate(Age_100 = (Age_months*100/max_age))
growth_curve_B.physalus <- growth_curve_B.physalus %>% mutate(growth_max_age[2,2]) %>% mutate(Age_100 = (Age_months*100/max_age))

growth_curve_Ph.phocoena <- growth_curve_Ph.phocoena %>% mutate(growth_max_age[3,2]) %>% mutate(Age_100 = (Age_months*100/max_age))
growth_curve_St.attenuata <- growth_curve_St.attenuata  %>% mutate(growth_max_age[4,2]) %>% mutate(Age_100 = (Age_months*100/max_age))

#Ossification sequence
ossification_seq_B.bonaerensis <- ossification_seq_B.bonaerensis %>% mutate(growth_max_age[1,2]) %>% mutate(Age_100 = (Age_months*100/max_age))
ossification_seq_B.physalus <- ossification_seq_B.physalus %>% mutate(growth_max_age[2,2]) %>% mutate(Age_100 = (Age_months*100/max_age))

ossification_seq_Ph.phocoena <- ossification_seq_Ph.phocoena %>% mutate(growth_max_age[3,2]) %>% mutate(Age_100 = (Age_months*100/max_age))
ossification_seq_St.attenuata <- ossification_seq_St.attenuata  %>% mutate(growth_max_age[4,2]) %>% mutate(Age_100 = (Age_months*100/max_age))

#Reform original dataset with added column Age_100 for plots
ossification_seq <- bind_rows(ossification_seq_B.bonaerensis, ossification_seq_B.physalus,ossification_seq_Ph.phocoena,ossification_seq_St.attenuata)

#Palettes
mypalette_tableau20 <- as.matrix(ggthemes_data[["tableau"]][["color-palettes"]][["regular"]][["Tableau 20"]][["value"]])
image(1:20, 1, as.matrix(1:20), col = mypalette_tableau20, xlab = "Tableau20",
      ylab = "", xaxt = "n", yaxt = "n", bty = "n")

mypalette_Mysticeti <- c(mypalette_tableau20[2,], mypalette_tableau20[4,], mypalette_tableau20[6,], mypalette_tableau20[8,],
                         mypalette_tableau20[10,], mypalette_tableau20[12,], mypalette_tableau20[14,],mypalette_tableau20[16,], 
                         mypalette_tableau20[18,], mypalette_tableau20[20,])
image(1:10, 1, as.matrix(1:10), col = mypalette_Mysticeti, xlab = "Mysticeti",
      ylab = "", xaxt = "n", yaxt = "n", bty = "n")

mypalette_Odontoceti <- c(mypalette_tableau20[1,], mypalette_tableau20[3,], mypalette_tableau20[5,], mypalette_tableau20[7,],
                          mypalette_tableau20[9,], mypalette_tableau20[11,], mypalette_tableau20[13,], mypalette_tableau20[15,], 
                          mypalette_tableau20[17,], mypalette_tableau20[19,])
image(1:10, 1, as.matrix(1:10), col = mypalette_Odontoceti, xlab = "Odontoceti",
      ylab = "", xaxt = "n", yaxt = "n", bty = "n")

mypalette_earbones <- c(mypalette_Mysticeti[2], mypalette_Mysticeti[3], mypalette_Odontoceti[1], mypalette_Odontoceti[6],mypalette_tableau20[13:14],"#ffffff", mypalette_Odontoceti[10])
mypalette_earbones_image <- image(1:8, 1, as.matrix(1:8), col = mypalette_earbones, xlab = "ear bones plots",
      ylab = "", xaxt = "n", yaxt = "n", bty = "n")

#Images for plots
B.bonaerensis <- readPNG("Data/b.bona.png")
B.physalus <- readPNG("Data/b.physalus.png")
Ph.phocoena <- readPNG("Data/phocoena.png")
St.attenuata <- readPNG("Data/stenella.png")

#PLOTS FOR OSSIFICATION EVENTS ON GROWTH CURVE ----

#Stages data frame for plot lines
stages <- data.frame(Age_months = c(2, 12/2), Age_100 = c(20,50)) #embryo fixed to about 2 months

#Palette for fill bone type
#Check order species and bones in dataset first
as.factor(ossification_seq$bone)
#Match with order species as selected from mypalette_earbones
as.factor(growth_curve$taxon)

mypalette_earbones_fill <- c(mypalette_earbones[2],mypalette_earbones[1], mypalette_earbones[3:4], 
                             mypalette_earbones[7],mypalette_earbones[7],mypalette_earbones[7],mypalette_earbones[7])

#Plot data species with selected model, specimens and events - % gestation to make plot better
ossification_earbones <- ggplot(ossification_seq, aes(y = TL_mm, x = Age_100, xend = 1, shape = state,  fill = bone, color = taxon)) + #xend useful to make sure graphs goes close to 0
  geom_point(size = 4)+ 
  stat_smooth(data = growth_curve_B.bonaerensis, aes(y = TL_mm, x = Age_100), method = lm, formula = y ~ poly(x, 2, raw = TRUE)-1, se = F, 
              inherit.aes = F, fullrange = T,
              linetype = 1, colour = mypalette_earbones[1], show.legend = F)+
  stat_smooth(data = growth_curve_B.physalus, aes(y = TL_mm, x = Age_100), method = lm, formula = y ~ poly(x, 2, raw = TRUE)-1, se = F, 
              inherit.aes = F, fullrange = T,
              linetype = 1, colour = mypalette_earbones[2], show.legend = F)+
  stat_smooth(data = growth_curve_Ph.phocoena, aes(y = TL_mm, x = Age_100), method = lm, formula = y ~ x, se = F, 
              inherit.aes = F, fullrange = T,
              linetype = 2, colour = mypalette_earbones[3], show.legend = F)+ 
  stat_smooth(data = growth_curve_St.attenuata, aes(y = TL_mm, x = Age_100), method = lm, formula = y ~ x, se = F, 
              inherit.aes = F, fullrange = T,
              linetype = 2, colour = mypalette_earbones[4], show.legend = F)+ 
  scale_color_manual(values = mypalette_earbones[1:4])+
  scale_shape_manual(name = "Events", labels  = c("bulla/periotic visible", "processes of bulla/periotic ossified", "bulla/periotic completely ossified"), 
                                                  values = c(21, 22, 24))+
  scale_fill_manual(values = mypalette_earbones_fill,  name = "Element", labels = c("Bulla", "Periotic"))+
  geom_vline(data = stages, aes(xintercept = Age_100), color = mypalette_earbones[5], linetype = 4)+
  theme_bw()+
  xlab("Age (% gestation)")+
  ylab("Total Length (mm)")+
  ggtitle("Ossification stages ear bones")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom", legend.direction = "vertical")+
  guides(fill = guide_legend(override.aes = list(fill = c("black","white"), shape = 21)), 
         colour = guide_legend(label = F, title = NULL, override.aes = list(shape = NA)))
ossification_earbones <- move_layers(ossification_earbones, "GeomPoint", position = "top")
ossification_earbones

#Add annotations growth stages
ossification_earbones <- ossification_earbones + 
  annotate("text", x = 10, y = 5400, label = "embryo", fontface = "italic", size = 4)+
  annotate("text", x = 35, y = 5400, label = "early fetus", fontface = "italic", size = 4)+
  annotate("text", x = 65, y = 5400, label = "late fetus", fontface = "italic", size = 4)
ossification_earbones

#Add silhouettes groups
ossification_earbones <- ossification_earbones + 
  add_phylopic(B.bonaerensis, alpha = 1, x = 80, y = 2800, ysize = 600, color = mypalette_earbones[1])+
  add_phylopic(B.physalus, alpha = 1, x = 65, y = 4000, ysize = 600, color = mypalette_earbones[2])+
  add_phylopic(Ph.phocoena, alpha = 1, x = 90, y = 200, ysize = 700, color = mypalette_earbones[3])+
  add_phylopic(St.attenuata, alpha = 1, x = 90, y = 1300, ysize = 600, color = mypalette_earbones[4])
ossification_earbones

#Plot log data both groups with selected model, specimens and events
ossification_earbones_log <- ggplot(ossification_seq, aes(y = TL_mm_log, x = Age_100, xend = 1, shape = state, fill = bone, color = taxon)) + #xend useful to make sure graphs goes close to 0
  geom_point(size = 4)+ 
  stat_smooth(data = growth_curve_B.bonaerensis, aes(y = TL_mm_log, x = Age_100), method = lm, formula = y ~ x, se = F, 
               inherit.aes = F, fullrange = T,
              linetype = 1, colour = mypalette_earbones[1], show.legend = F)+
  stat_smooth(data = growth_curve_B.physalus, aes(y = TL_mm_log, x = Age_100), method = lm, formula = y ~ x, se = F, 
               inherit.aes = F, fullrange = F,
              linetype = 1, colour = mypalette_earbones[2], show.legend = F)+ 
  stat_smooth(data = growth_curve_Ph.phocoena, aes(y = TL_mm_log, x = Age_100), method = lm, formula = y ~ x, se = F, 
              inherit.aes = F, fullrange = T,
              linetype = 2, colour = mypalette_earbones[3], show.legend = F)+
  stat_smooth(data = growth_curve_St.attenuata, aes(y = TL_mm_log, x = Age_100), method = lm, formula = y ~ x, se = F, 
              inherit.aes = F, fullrange = T,
              linetype = 2, colour = mypalette_earbones[4], show.legend = F)+
  scale_color_manual(values = mypalette_earbones[1:4])+
  scale_shape_manual(name = "Events", labels  = c("bulla/periotic visible", "processes of bulla/periotic ossified", "bulla/periotic completely ossified"), 
                     values = c(21, 22, 24))+
  scale_fill_manual(values = mypalette_earbones_fill,  name = "Element", labels = c("Bulla", "Periotic"))+
  geom_vline(data = stages, aes(xintercept = Age_100), color = mypalette_earbones[5], linetype = 4)+
  theme_bw()+
  xlab("Age (% gestation)")+
  ylab("Log(Total Length)")+
  ggtitle("Ossification stages ear bones")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom", legend.direction = "vertical")+
  guides(fill = guide_legend(override.aes = list(fill = c("black","white"), shape = 21)), colour = guide_legend(label = F, title = NULL, override.aes = list(shape = NA)))
ossification_earbones_log <- move_layers(ossification_earbones_log, "GeomPoint", position = "top")
ossification_earbones_log

#Add annotations growth stages
ossification_earbones_log <- ossification_earbones_log + 
  annotate("text", x = 10, y = 8.6, label = "embryo", fontface = "italic", size = 4)+
  annotate("text", x = 35, y = 8.6, label = "early fetus", fontface = "italic", size = 4)+
  annotate("text", x = 65, y = 8.6, label = "late fetus", fontface = "italic", size = 4)
ossification_earbones_log

#Add silhouettes groups
ossification_earbones_log <- ossification_earbones_log + 
  add_phylopic(B.bonaerensis, alpha = 1, x = 85, y = 8, ysize = 4.2, color = mypalette_earbones[1])+
  add_phylopic(B.physalus, alpha = 1, x = 15, y = 6.5, ysize = 4.8, color = mypalette_earbones[2])+
  add_phylopic(Ph.phocoena, alpha = 1, x = 85, y = 6, ysize = 6, color = mypalette_earbones[3])+
  add_phylopic(St.attenuata, alpha = 1, x = 25, y = 4, ysize = 4.5, color = mypalette_earbones[4])
ossification_earbones_log

#ALLOMETRY OF BULLA AND PERIOTIC MEASURAMENTS FOR EACH GROUP (2-way ANOVA) ----
##Create basic models with no interactions with groups for comparison
allometry_bullaL_log <- lm(bullaL_log ~ BZW_log, data  = bulla_meas)
allometry_bullaW_log <- lm(bullaW_log ~ BZW_log, data  = bulla_meas)
allometry_perioticL_log <- lm(perioticL_log ~ BZW_log, data  = periotic_meas)
allometry_perioticW_log <- lm(perioticW_log ~ BZW_log, data  = periotic_meas)

#Save regressions to file
sink("allometry/allometry_log.txt", append = F)
print("Bulla length")
print(summary(allometry_bullaL_log))
print(anova(allometry_bullaL_log))

print("Bulla width")
print(summary(allometry_bullaW_log))
print(anova(allometry_bullaW_log))

print("Periotic length")
print(summary(allometry_perioticL_log))
print(anova(allometry_perioticL_log))

print("Periotic width")
print(summary(allometry_perioticW_log))
print(anova(allometry_perioticW_log))
sink()

##Regression of log values with group effect - check differences between groups
allometry_bullaL_log_int <- lm(bullaL_log ~ BZW_log * group, data  = bulla_meas)
allometry_bullaW_log_int <- lm(bullaW_log ~ BZW_log * group, data  = bulla_meas)
allometry_perioticL_log_int <- lm(perioticL_log ~ BZW_log * group, data  = periotic_meas)
allometry_perioticW_log_int <- lm(perioticW_log ~ BZW_log * group, data  = periotic_meas)

#Check results
summary(allometry_bullaL_log_int)
anova(allometry_bullaL_log_int)
summary(allometry_bullaW_log_int)
anova(allometry_bullaW_log_int)
summary(allometry_perioticL_log_int)
anova(allometry_perioticL_log_int)
summary(allometry_perioticW_log_int)
anova(allometry_perioticW_log_int)

#Plot diagnostics - look for no pattern in left plots, no outliers in bottom right and not too much deviation from the dotted line in top right
autoplot(allometry_bullaL_log_int, smooth.colour = NA)
autoplot(allometry_bullaW_log_int, smooth.colour = NA)
autoplot(allometry_perioticL_log_int, smooth.colour = NA)
autoplot(allometry_perioticW_log_int, smooth.colour = NA)

#Save regressions to file
sink("allometry/allometry_log_groups_int.txt", append = F)
print("Bulla length")
print(summary(allometry_bullaL_log_int))
print(anova(allometry_bullaL_log_int))

print("Bulla width")
print(summary(allometry_bullaW_log_int))
print(anova(allometry_bullaW_log_int))

print("Periotic length")
print(summary(allometry_perioticL_log_int))
print(anova(allometry_perioticL_log_int))

print("Periotic width")
print(summary(allometry_perioticW_log_int))
print(anova(allometry_perioticW_log_int))
sink()

##Check if groups can be different only in size (intercept) and not slope
allometry_bullaL_log_comb_taxa <- lm(bullaL_log ~ BZW_log + group, data  = bulla_meas)
allometry_bullaW_log_comb_taxa <- lm(bullaW_log ~ BZW_log + group, data  = bulla_meas)
allometry_perioticL_log_comb_taxa <- lm(perioticL_log ~ BZW_log + group, data  = periotic_meas)
allometry_perioticW_log_comb_taxa <- lm(perioticW_log ~ BZW_log + group, data  = periotic_meas)

#Check results
summary(allometry_bullaL_log_comb)
anova(allometry_bullaL_log_comb)
summary(allometry_bullaW_log_comb)
anova(allometry_bullaW_log_comb)
summary(allometry_perioticL_log_comb)
anova(allometry_perioticL_log_comb)
summary(allometry_perioticW_log_comb)
anova(allometry_perioticW_log_comb)

#Save regressions to file
sink("allometry/allometry_log_groups_comb.txt", append = F)
print("Bulla length")
print(summary(allometry_bullaL_log_comb))
print(anova(allometry_bullaL_log_comb))

print("Bulla width")
print(summary(allometry_bullaW_log_comb))
print(anova(allometry_bullaW_log_comb))

print("Periotic length")
print(summary(allometry_perioticL_log_comb))
print(anova(allometry_perioticL_log_comb))

print("Periotic width")
print(summary(allometry_perioticW_log_comb))
print(anova(allometry_perioticW_log_comb))
sink()

##Compare models with ANOVA and calculate AIC scores for natural values
models_bullaL_log <- list(allometry_bullaL_log, allometry_bullaL_log_comb, allometry_bullaL_log_int)
models_bullaW_log <- list(allometry_bullaW_log, allometry_bullaW_log_comb, allometry_bullaW_log_int)
models_perioticL_log <- list(allometry_perioticL_log, allometry_perioticL_log_comb, allometry_perioticL_log_int)
models_perioticW_log <- list(allometry_perioticW_log, allometry_perioticW_log_comb, allometry_perioticW_log_int)

model_names <- c("na", "comb", "int")

#ANOVAs - is a model significantly better than the others?
anova(allometry_bullaL_log, allometry_bullaL_log_comb, allometry_bullaL_log_int)
anova(allometry_bullaW_log, allometry_bullaW_log_comb, allometry_bullaW_log_int)
anova(allometry_perioticL_log, allometry_perioticL_log_comb, allometry_perioticL_log_int)
anova(allometry_perioticW_log, allometry_perioticW_log_comb, allometry_perioticW_log_int)

#AICc score - what is the best model? 
aictab(cand.set = models_bullaL_log, modnames = model_names)
aictab(cand.set = models_bullaW_log, modnames = model_names)
aictab(cand.set = models_perioticL_log, modnames = model_names)
aictab(cand.set = models_perioticW_log, modnames = model_names)

#Save results to file
sink("allometry/allometry_log_models.txt", append = F)
print("Bulla length")
anova(allometry_bullaL_log, allometry_bullaL_log_comb, allometry_bullaL_log_int)
aictab(cand.set = models_bullaL_log, modnames = model_names)

print("Bulla width")
anova(allometry_bullaW_log, allometry_bullaW_log_comb, allometry_bullaW_log_int)
aictab(cand.set = models_bullaW_log, modnames = model_names)

print("Periotic length")
anova(allometry_perioticL_log, allometry_perioticL_log_comb, allometry_perioticL_log_int)
aictab(cand.set = models_perioticL_log, modnames = model_names)

print("Periotic width")
anova(allometry_perioticW_log, allometry_perioticW_log_comb, allometry_perioticW_log_int)
aictab(cand.set = models_perioticW_log, modnames = model_names)
sink()

##Plots group ----
#Plot regression lines per group
#Add confidence intervals
#Create data for confidence intervals
bullaL_log_newX <- expand.grid(BZW_log = seq(from = min(bulla_meas$BZW_log), to = max(bulla_meas$BZW_log), length.out = 30), #use min and max of x values (BZW) as limits and use number of specimens as length of sequence
                               group = c("Mysticeti", "Odontoceti"))         #warp x_vals on values of x axis (BZW), add groupings
bullaL_log_newY <- predict(allometry_bullaL_log_int, newdata = bullaL_log_newX, interval="confidence",
                           level = 0.95)
#Make data frame of data for confidence intervals
bullaL_log_conf_intervals <- data.frame(bullaL_log_newX, bullaL_log_newY)
#Rename columns to match main plot tibble variables for x and y
bullaL_log_conf_intervals <- rename(bullaL_log_conf_intervals, bullaL_log = fit)
bullaL_log_conf_intervals

#Plot
allometry_bullaL_log_int_plot <- ggplot(bulla_meas, aes(y = bullaL_log, x = BZW_log, fill = group, color = group)) +
  geom_smooth(data = bullaL_log_conf_intervals, aes(ymin = lwr, ymax = upr, fill = group, colour = group, linetype = group), stat = 'identity',          #confidence intervals and reg line, before points
              size = 0.8, alpha = 0.3, show.legend = F)+      #put col and other graphics OUTSIDE of aes()!!!
  geom_point(size = 3, alpha = 0.2)+       #points after, so they are on top
  scale_color_manual(name = "Groups", labels  = c("Mysticeti", "Odontoceti"), values = c(mypalette_earbones[1],mypalette_earbones[3]), #select colors from palette from taxa
                     aesthetics = c("color","fill"))+         
  theme_classic(base_size = 12)+
  xlab("Log(BZW)")+
  ylab("Log(Bulla length) *")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11), legend.position = "none", legend.direction = "vertical", 
        axis.title.x = element_text(vjust = -1), 
        axis.title.y = element_text(face = "bold", vjust = 2, size = 14))
allometry_bullaL_log_int_plot <- move_layers(allometry_bullaL_log_int_plot, "GeomPoint", position = "top")
allometry_bullaL_log_int_plot <- allometry_bullaL_log_int_plot + 
  add_phylopic(B.bonaerensis, alpha = 1, x = 4.5, y = 4, ysize = 0.28, color = mypalette_earbones[1])+
  add_phylopic(St.attenuata, alpha = 1, x = 5.2, y = 3, ysize = 0.25, color = mypalette_earbones[3])
allometry_bullaL_log_int_plot

# ggtitle ("Log-transformed BZW vs Bulla length by group - p-value < 0.001***")+  #copy from model summary

#Add confidence intervals
#Create data for confidence intervals
bullaW_log_newX <- expand.grid(BZW_log = seq(from = min(bulla_meas$BZW_log), to = max(bulla_meas$BZW_log), length.out = 30), #use min and max of x values (BZW) as limits and use number of specimens as length of sequence
                               group = c("Mysticeti", "Odontoceti"))         #warp x_vals on values of x axis (BZW), add groupings
bullaW_log_newY <- predict(allometry_bullaW_log_int, newdata = bullaW_log_newX, interval="confidence",
                           level = 0.95)
#Make data frame of data for confidence intervals
bullaW_log_conf_intervals <- data.frame(bullaW_log_newX, bullaW_log_newY)
#Rename columns to match main plot tibble variables for x and y
bullaW_log_conf_intervals <- rename(bullaW_log_conf_intervals, bullaW_log = fit)
bullaW_log_conf_intervals

#Plot
allometry_bullaW_log_int_plot <- ggplot(bulla_meas, aes(y = bullaW_log, x = BZW_log, fill = group, color = group)) +
  geom_smooth(data = bullaW_log_conf_intervals, aes(ymin = lwr, ymax = upr, fill = group, colour = group, linetype = group), stat = 'identity',          #confidence intervals and reg line, before points
              size = 0.8, alpha = 0.3, show.legend = F)+      #put col and other graphics OUTSIDE of aes()!!!
  geom_point(size = 3, alpha = 0.2)+       #points after, so they are on top
  scale_color_manual(name = "Groups", labels  = c("Mysticeti", "Odontoceti"), values = c(mypalette_earbones[1],mypalette_earbones[3]), #select colors from palette from taxa
                     aesthetics = c("color","fill"))+         
  theme_classic(base_size = 12)+
  xlab("Log(BZW)")+
  ylab("Log(Bulla width) *")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11), legend.position = "none", legend.direction = "vertical", 
        axis.title.x = element_text(vjust = -1), 
        axis.title.y = element_text(face = "bold", vjust = 2, size = 14))
allometry_bullaW_log_int_plot <- move_layers(allometry_bullaW_log_int_plot, "GeomPoint", position = "top")
allometry_bullaW_log_int_plot <- allometry_bullaW_log_int_plot + 
  add_phylopic(B.bonaerensis, alpha = 1, x = 4.5, y = 3.8, ysize = 0.28, color = mypalette_earbones[1])+
  add_phylopic(St.attenuata, alpha = 1, x = 5.2, y = 2.7, ysize = 0.25, color = mypalette_earbones[3])
allometry_bullaW_log_int_plot

#ggtitle ("Log-transformed BZW vs Bulla width by group - p-value < 0.001***")+  #copy from model summary

#Add confidence intervals
#Create data for confidence intervals
perioticL_log_newX <- expand.grid(BZW_log = seq(from = min(periotic_meas$BZW_log), to = max(periotic_meas$BZW_log), length.out = 21), #use min and max of x values (BZW) as limits and use number of specimens as length of sequence
                                  group = c("Mysticeti", "Odontoceti"))         #warp x_vals on values of x axis (BZW), add groupings
perioticL_log_newY <- predict(allometry_perioticL_log_int, newdata = perioticL_log_newX, interval="confidence",
                              level = 0.95)
#Make data frame of data for confidence intervals
perioticL_log_conf_intervals <- data.frame(perioticL_log_newX, perioticL_log_newY)
#Rename columns to match main plot tibble variables for x and y
perioticL_log_conf_intervals <- rename(perioticL_log_conf_intervals, perioticL_log = fit)
perioticL_log_conf_intervals

#Plot
allometry_perioticL_log_int_plot <- ggplot(periotic_meas, aes(y = perioticL_log, x = BZW_log, fill = group, color = group)) +
  geom_smooth(data = perioticL_log_conf_intervals, aes(ymin = lwr, ymax = upr, fill = group, colour = group, linetype = group), stat = 'identity',          #confidence intervals and reg line, before points
              size = 0.8, alpha = 0.3, show.legend = F)+      #put col and other graphics OUTSIDE of aes()!!!
  geom_point(size = 3, alpha = 0.2)+       #points after, so they are on top
  scale_color_manual(name = "Groups", labels  = c("Mysticeti", "Odontoceti"), values = c(mypalette_earbones[1],mypalette_earbones[3]), #select colors from palette from taxa
                     aesthetics = c("color","fill"))+         
  theme_classic(base_size = 12)+
  xlab("Log(BZW)")+
  ylab("Log(Periotic length) *")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11), legend.position = "none", legend.direction = "vertical", 
        axis.title.x = element_text(vjust = -1), 
        axis.title.y = element_text(face = "bold", vjust = 2, size = 14))
allometry_perioticL_log_int_plot <- move_layers(allometry_perioticL_log_int_plot, "GeomPoint", position = "top")
allometry_perioticL_log_int_plot <- allometry_perioticL_log_int_plot + 
  add_phylopic(B.bonaerensis, alpha = 1, x = 4.7, y = 4.2, ysize = 0.28, color = mypalette_earbones[1])+
  add_phylopic(St.attenuata, alpha = 1, x = 5.5, y = 3.2, ysize = 0.25, color = mypalette_earbones[3])
allometry_perioticL_log_int_plot

#ggtitle ("Log-transformed BZW vs Periotic length by group - p-value < 0.001***")+  #copy from model summary

#Add confidence intervals
#Create data for confidence intervals
perioticW_log_newX <- expand.grid(BZW_log = seq(from = min(periotic_meas$BZW_log), to = max(periotic_meas$BZW_log), length.out = 21), #use min and max of x values (BZW) as limits and use number of specimens as length of sequence
                                  group = c("Mysticeti", "Odontoceti"))         #warp x_vals on values of x axis (BZW), add groupings
perioticW_log_newY <- predict(allometry_perioticW_log_int, newdata = perioticW_log_newX, interval="confidence",
                              level = 0.95)
#Make data frame of data for confidence intervals
perioticW_log_conf_intervals <- data.frame(perioticW_log_newX, perioticW_log_newY)
#Rename columns to match main plot tibble variables for x and y
perioticW_log_conf_intervals <- rename(perioticW_log_conf_intervals, perioticW_log = fit)
perioticW_log_conf_intervals

#Plot
allometry_perioticW_log_int_plot <- ggplot(periotic_meas, aes(y = perioticW_log, x = BZW_log, fill = group, color = group)) +
  geom_smooth(data = perioticW_log_conf_intervals, aes(ymin = lwr, ymax = upr, fill = group, colour = group, linetype = group), stat = 'identity',          #confidence intervals and reg line, before points
              size = 0.8, alpha = 0.3, show.legend = F)+      #put col and other graphics OUTSIDE of aes()!!!
  geom_point(size = 3, alpha = 0.2)+       #points after, so they are on top
  scale_color_manual(name = "Groups", labels  = c("Mysticeti", "Odontoceti"), values = c(mypalette_earbones[1],mypalette_earbones[3]), #select colors from palette from taxa
                     aesthetics = c("color","fill"))+         
  theme_classic(base_size = 12)+
  xlab("Log(BZW)")+
  ylab("Log(Periotic width) ***")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11), legend.position = "none", legend.direction = "vertical", 
        axis.title.x = element_text(vjust = -1), 
        axis.title.y = element_text(face = "bold", vjust = 2, size = 14))
allometry_perioticW_log_int_plot <- move_layers(allometry_perioticW_log_int_plot, "GeomPoint", position = "top")
allometry_perioticW_log_int_plot <- allometry_perioticW_log_int_plot + 
  add_phylopic(B.bonaerensis, alpha = 1, x = 4.7, y = 3.6, ysize = 0.18, color = mypalette_earbones[1])+
  add_phylopic(St.attenuata, alpha = 1, x = 5.5, y = 2.8, ysize = 0.16, color = mypalette_earbones[3])
allometry_perioticW_log_int_plot

#ggtitle ("Log-transformed BZW vs Periotic width by group - p-value < 0.001*** (group p < 0.001***)")+  #copy from model summary

#All best models plots together
grid.arrange(allometry_bullaL_log_int_plot,  allometry_perioticL_log_int_plot, allometry_bullaW_log_int_plot, allometry_perioticW_log_int_plot)

#ALLOMETRY OF BULLA AND PERIOTIC MEASURAMENTS FOR EACH TAXON (2-way ANOVA) ----
##Create basic models with no interactions with groups for comparison
allometry_bullaL_log_taxa <- lm(bullaL_log ~ BZW_log, data  = bulla_meas_taxa)
allometry_bullaW_log_taxa <- lm(bullaW_log ~ BZW_log, data  = bulla_meas_taxa)
allometry_perioticL_log_taxa <- lm(perioticL_log ~ BZW_log, data  = periotic_meas_taxa)
allometry_perioticW_log_taxa <- lm(perioticW_log ~ BZW_log, data  = periotic_meas_taxa)

#Save regressions to file
sink("allometry/allometry_log_taxa.txt", append = F)
print("Bulla length")
print(summary(allometry_bullaL_log_taxa))
print(anova(allometry_bullaL_log_taxa))

print("Bulla width")
print(summary(allometry_bullaW_log_taxa))
print(anova(allometry_bullaW_log_taxa))

print("Periotic length")
print(summary(allometry_perioticL_log_taxa))
print(anova(allometry_perioticL_log_taxa))

print("Periotic width")
print(summary(allometry_perioticW_log_taxa))
print(anova(allometry_perioticW_log_taxa))
sink()

##Regression of log values with taxon effect - check differences between taxa
allometry_bullaL_log_int_taxa <- lm(bullaL_log ~ BZW_log * taxon, data  = bulla_meas_taxa)
allometry_bullaW_log_int_taxa <- lm(bullaW_log ~ BZW_log * taxon, data  = bulla_meas_taxa)
allometry_perioticL_log_int_taxa <- lm(perioticL_log ~ BZW_log * taxon, data  = periotic_meas_taxa)
allometry_perioticW_log_int_taxa <- lm(perioticW_log ~ BZW_log * taxon, data  = periotic_meas_taxa)

#Check results
summary(allometry_bullaL_log_int_taxa)
anova(allometry_bullaL_log_int_taxa)
summary(allometry_bullaW_log_int_taxa)
anova(allometry_bullaW_log_int_taxa)
summary(allometry_perioticL_log_int_taxa)
anova(allometry_perioticL_log_int_taxa)
summary(allometry_perioticW_log_int_taxa)
anova(allometry_perioticW_log_int_taxa)

#Plot diagnostics - look for no pattern in left plots, no outliers in bottom right and not too much deviation from the dotted line in top right
autoplot(allometry_bullaL_log_int_taxa, smooth.colour = NA)
autoplot(allometry_bullaW_log_int_taxa, smooth.colour = NA)
autoplot(allometry_perioticL_log_int_taxa, smooth.colour = NA)
autoplot(allometry_perioticW_log_int_taxa, smooth.colour = NA)

#Save regressions to file
sink("allometry/allometry_log_taxa_int.txt", append = F)
print("Bulla length")
print(summary(allometry_bullaL_log_int_taxa))
print(anova(allometry_bullaL_log_int_taxa))

print("Bulla width")
print(summary(allometry_bullaW_log_int_taxa))
print(anova(allometry_bullaW_log_int_taxa))

print("Periotic length")
print(summary(allometry_perioticL_log_int_taxa))
print(anova(allometry_perioticL_log_int_taxa))

print("Periotic width")
print(summary(allometry_perioticW_log_int_taxa))
print(anova(allometry_perioticW_log_int_taxa))
sink()

##Check if taxa can be different only in size (intercept) and not slope
allometry_bullaL_log_comb_taxa <- lm(bullaL_log ~ BZW_log + taxon, data  = bulla_meas_taxa)
allometry_bullaW_log_comb_taxa <- lm(bullaW_log ~ BZW_log + taxon, data  = bulla_meas_taxa)
allometry_perioticL_log_comb_taxa <- lm(perioticL_log ~ BZW_log + taxon, data  = periotic_meas_taxa)
allometry_perioticW_log_comb_taxa <- lm(perioticW_log ~ BZW_log + taxon, data  = periotic_meas_taxa)

#Check results
summary(allometry_bullaL_log_comb_taxa)
anova(allometry_bullaL_log_comb_taxa)
summary(allometry_bullaW_log_comb_taxa)
anova(allometry_bullaW_log_comb_taxa)
summary(allometry_perioticL_log_comb_taxa)
anova(allometry_perioticL_log_comb_taxa)
summary(allometry_perioticW_log_comb_taxa)
anova(allometry_perioticW_log_comb_taxa)

#Save regressions to file
sink("allometry/allometry_log_taxa_comb.txt", append = F)
print("Bulla length")
print(summary(allometry_bullaL_log_comb_taxa))
print(anova(allometry_bullaL_log_comb_taxa))

print("Bulla width")
print(summary(allometry_bullaW_log_comb_taxa))
print(anova(allometry_bullaW_log_comb_taxa))

print("Periotic length")
print(summary(allometry_perioticL_log_comb_taxa))
print(anova(allometry_perioticL_log_comb_taxa))

print("Periotic width")
print(summary(allometry_perioticW_log_comb_taxa))
print(anova(allometry_perioticW_log_comb_taxa))
sink()

##Compare models with ANOVA and calculate AIC scores for natural values
models_bullaL_log_taxa <- list(allometry_bullaL_log_taxa, allometry_bullaL_log_comb_taxa, allometry_bullaL_log_int_taxa)
models_bullaW_log_taxa <- list(allometry_bullaW_log_taxa, allometry_bullaW_log_comb_taxa, allometry_bullaW_log_int_taxa)
models_perioticL_log_taxa <- list(allometry_perioticL_log_taxa, allometry_perioticL_log_comb_taxa, allometry_perioticL_log_int_taxa)
models_perioticW_log_taxa <- list(allometry_perioticW_log_taxa, allometry_perioticW_log_comb_taxa, allometry_perioticW_log_int_taxa)

#ANOVAs - is a model significantly better than the others?
anova(allometry_bullaL_log_taxa, allometry_bullaL_log_comb_taxa, allometry_bullaL_log_int_taxa)
anova(allometry_bullaW_log_taxa, allometry_bullaW_log_comb_taxa, allometry_bullaW_log_int_taxa)
anova(allometry_perioticL_log_taxa, allometry_perioticL_log_comb_taxa, allometry_perioticL_log_int_taxa)
anova(allometry_perioticW_log_taxa, allometry_perioticW_log_comb_taxa, allometry_perioticW_log_int_taxa)

#AICc score - what is the best model? 
aictab(cand.set = models_bullaL_log_taxa, modnames = model_names)
aictab(cand.set = models_bullaW_log_taxa, modnames = model_names)
aictab(cand.set = models_perioticL_log_taxa, modnames = model_names)
aictab(cand.set = models_perioticW_log_taxa, modnames = model_names)

#Save results to file
sink("allometry/allometry_log_models_taxa.txt", append = F)
print("Bulla length")
anova(allometry_bullaL_log_taxa, allometry_bullaL_log_comb_taxa, allometry_bullaL_log_int_taxa)
aictab(cand.set = models_bullaL_log_taxa, modnames = model_names)

print("Bulla width")
anova(allometry_bullaW_log_taxa, allometry_bullaW_log_comb_taxa, allometry_bullaW_log_int_taxa)
aictab(cand.set = models_bullaW_log_taxa, modnames = model_names)

print("Periotic length")
anova(allometry_perioticL_log_taxa, allometry_perioticL_log_comb_taxa, allometry_perioticL_log_int_taxa)
aictab(cand.set = models_perioticL_log_taxa, modnames = model_names)

print("Periotic width")
anova(allometry_perioticW_log_taxa, allometry_perioticW_log_comb_taxa, allometry_perioticW_log_int_taxa)
aictab(cand.set = models_perioticW_log_taxa, modnames = model_names)
sink()

##Plots taxa ----
#Plot regression lines per taxon
#Add confidence intervals
#Get taxa names to copy for new conf intervals dataset
as.factor(bulla_meas_taxa$taxon)

#Create group column of equal lenght of new con intervals data frame
myst <- data.frame(group = (rep("Mysticeti", len = 60)))
odont <- data.frame(group = (rep("Odontoceti", len = 60)))
groups_taxa <- rbind(myst, odont)

#Create data for confidence intervals
bullaL_log_newX_taxa <- expand.grid(BZW_log = seq(from = min(bulla_meas_taxa$BZW_log), to = max(bulla_meas_taxa$BZW_log), length.out = 30), #use min and max of x values (BZW) as limits and use number of specimens as length of sequence
                               taxon = c("B.bonaerensis", "B.physalus", "Ph.phocoena", "St.attenuata"))         #warp x_vals on values of x axis (BZW), add groupings
bullaL_log_newY_taxa <- predict(allometry_bullaL_log_int_taxa, newdata = bullaL_log_newX_taxa, interval="confidence",
                           level = 0.95)
#Make data frame of data for confidence intervals
bullaL_log_conf_intervals_taxa <- data.frame(bullaL_log_newX_taxa, bullaL_log_newY_taxa, groups_taxa)
#Rename columns to match main plot tibble variables for x and y
bullaL_log_conf_intervals_taxa <- rename(bullaL_log_conf_intervals_taxa, bullaL_log = fit)
bullaL_log_conf_intervals_taxa

#Plot
allometry_bullaL_log_int_taxa_plot <- ggplot(bulla_meas_taxa, aes(y = bullaL_log, x = BZW_log, fill = taxon, color = taxon, linetype = group)) +
  geom_smooth(data = bullaL_log_conf_intervals_taxa, aes(ymin = lwr, ymax = upr, fill = taxon, colour = taxon, linetype = group), stat = 'identity',          #confidence intervals and reg line, before points
              size = 0.8, alpha = 0.3, show.legend = F)+      #put col and other graphics OUTSIDE of aes()!!!
  geom_point(size = 3, alpha = 0.2)+       #points after, so they are on top
  scale_color_manual(name = "Taxa", labels  = c("B.bonaerensis", "B.physalus", "Ph.phocoena", "St.attenuata"), values = c(mypalette_earbones[1:4]), #select colors from palette from taxa
                     aesthetics = c("color","fill"))+         
  theme_classic(base_size = 12)+
  xlab("Log(BZW)")+
  ylab("Log(Bulla length) n.s.")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11), legend.position = "none", legend.direction = "vertical", 
        axis.title.x = element_text(vjust = -1), 
        axis.title.y = element_text(face = "bold", vjust = 2, size = 14))
allometry_bullaL_log_int_taxa_plot <- move_layers(allometry_bullaL_log_int_taxa_plot, "GeomPoint", position = "top")
allometry_bullaL_log_int_taxa_plot <- allometry_bullaL_log_int_taxa_plot + 
  add_phylopic(B.bonaerensis, alpha = 1, x = 5.5, y = 2.5, ysize = 0.25, color = mypalette_earbones[1])+
  add_phylopic(B.physalus, alpha = 1, x = 5.5, y = 2.1, ysize = 0.25, color = mypalette_earbones[2])+
  add_phylopic(Ph.phocoena, alpha = 1, x = 5.5, y = 1.7, ysize = 0.35, color = mypalette_earbones[3])+
  add_phylopic(St.attenuata, alpha = 1, x = 5.5, y = 1.3, ysize = 0.28, color = mypalette_earbones[4])
allometry_bullaL_log_int_taxa_plot

# ggtitle ("Log-transformed BZW vs Bulla length by group - p-value < 0.001***")+  #copy from model summary

#Add confidence intervals
#Create data for confidence intervals
bullaW_log_newX_taxa <- expand.grid(BZW_log = seq(from = min(bulla_meas_taxa$BZW_log), to = max(bulla_meas_taxa$BZW_log), length.out = 30), #use min and max of x values (BZW) as limits and use number of specimens as length of sequence
                                    taxon = c("B.bonaerensis", "B.physalus", "Ph.phocoena", "St.attenuata"))         #warp x_vals on values of x axis (BZW), add groupings
bullaW_log_newY_taxa <- predict(allometry_bullaW_log_int_taxa, newdata = bullaW_log_newX_taxa, interval="confidence",
                                level = 0.95)
#Make data frame of data for confidence intervals
bullaW_log_conf_intervals_taxa <- data.frame(bullaW_log_newX_taxa, bullaW_log_newY_taxa, groups_taxa)
#Rename columns to match main plot tibble variables for x and y
bullaW_log_conf_intervals_taxa <- rename(bullaW_log_conf_intervals_taxa, bullaW_log = fit)
bullaW_log_conf_intervals_taxa

#Plot
allometry_bullaW_log_int_taxa_plot <- ggplot(bulla_meas_taxa, aes(y = bullaW_log, x = BZW_log, fill = taxon, color = taxon, linetype = group)) +
  geom_smooth(data = bullaW_log_conf_intervals_taxa, aes(ymin = lwr, ymax = upr, fill = taxon, colour = taxon, linetype = group), stat = 'identity',          #confidence intervals and reg line, before points
              size = 0.8, alpha = 0.3, show.legend = F)+      #put col and other graphics OUTSIDE of aes()!!!
  geom_point(size = 3, alpha = 0.2)+       #points after, so they are on top
  scale_color_manual(name = "Taxa", labels  = c("B.bonaerensis", "B.physalus", "Ph.phocoena", "St.attenuata"), values = c(mypalette_earbones[1:4]), #select colors from palette from taxa
                     aesthetics = c("color","fill"))+         
  theme_classic(base_size = 12)+
  xlab("Log(BZW)")+
  ylab("Log(Bulla width) ***")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11), legend.position = "none", legend.direction = "vertical", 
        axis.title.x = element_text(vjust = -1), 
        axis.title.y = element_text(face = "bold", vjust = 2, size = 14))
allometry_bullaW_log_int_taxa_plot <- move_layers(allometry_bullaW_log_int_taxa_plot, "GeomPoint", position = "top")
allometry_bullaW_log_int_taxa_plot <- allometry_bullaW_log_int_taxa_plot + 
  add_phylopic(B.bonaerensis, alpha = 1, x = 5.5, y = 1.7, ysize = 0.3, color = mypalette_earbones[1])+
  add_phylopic(B.physalus, alpha = 1, x = 4.8, y = 4.1, ysize = 0.32, color = mypalette_earbones[2])+
  add_phylopic(Ph.phocoena, alpha = 1, x = 5.5, y = 1.2, ysize = 0.42, color = mypalette_earbones[3])+
  add_phylopic(St.attenuata, alpha = 1, x = 5.5, y = 0.7, ysize = 0.32, color = mypalette_earbones[4])
allometry_bullaW_log_int_taxa_plot

#ggtitle ("Log-transformed BZW vs Bulla width by group - p-value < 0.001***")+  #copy from model summary

#Add confidence intervals
#Create group column of equal lenght of new con intervals data frame
myst_p <- data.frame(group = (rep("Mysticeti", len = 30)))
groups_taxa_p <- rbind(myst_p, odont)

#Create data for confidence intervals
perioticL_log_newX_taxa <- expand.grid(BZW_log = seq(from = min(periotic_meas_taxa$BZW_log), to = max(periotic_meas_taxa$BZW_log), length.out = 30), #use min and max of x values (BZW) as limits and use number of specimens as length of sequence
                                    taxon = c("B.bonaerensis", "Ph.phocoena", "St.attenuata"))         #warp x_vals on values of x axis (BZW), add groupings
perioticL_log_newY_taxa <- predict(allometry_perioticL_log_int_taxa, newdata = perioticL_log_newX_taxa, interval="confidence",
                                level = 0.95)
#Make data frame of data for confidence intervals
perioticL_log_conf_intervals_taxa <- data.frame(perioticL_log_newX_taxa, perioticL_log_newY_taxa, groups_taxa_p)
#Rename columns to match main plot tibble variables for x and y
perioticL_log_conf_intervals_taxa <- rename(perioticL_log_conf_intervals_taxa, perioticL_log = fit)
perioticL_log_conf_intervals_taxa

#Plot
allometry_perioticL_log_int_taxa_plot <- ggplot(periotic_meas_taxa, aes(y = perioticL_log, x = BZW_log, fill = taxon, color = taxon, linetype = group)) +
  geom_smooth(data = perioticL_log_conf_intervals_taxa, aes(ymin = lwr, ymax = upr, fill = taxon, colour = taxon, linetype = group), stat = 'identity',          #confidence intervals and reg line, before points
              size = 0.8, alpha = 0.3, show.legend = F)+      #put col and other graphics OUTSIDE of aes()!!!
  geom_point(size = 3, alpha = 0.2)+       #points after, so they are on top
  scale_color_manual(name = "Taxa", labels  = c("B.bonaerensis", "Ph.phocoena", "St.attenuata"), values = c(mypalette_earbones[1], mypalette_earbones[3:4]), #select colors from palette from taxa
                     aesthetics = c("color","fill"))+         
  theme_classic(base_size = 12)+
  xlab("Log(BZW)")+
  ylab("Log(Periotic length) n.s.")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11), legend.position = "none", legend.direction = "vertical", 
        axis.title.x = element_text(vjust = -1), 
        axis.title.y = element_text(face = "bold", vjust = 2, size = 14))
allometry_perioticL_log_int_taxa_plot <- move_layers(allometry_perioticL_log_int_taxa_plot, "GeomPoint", position = "top")
allometry_perioticL_log_int_taxa_plot <- allometry_perioticL_log_int_taxa_plot + 
  add_phylopic(B.bonaerensis, alpha = 1, x = 5.5, y = 2.2, ysize = 0.26, color = mypalette_earbones[1])+
  add_phylopic(Ph.phocoena, alpha = 1, x = 5.5, y = 1.7, ysize = 0.37, color = mypalette_earbones[3])+
  add_phylopic(St.attenuata, alpha = 1, x = 5.5, y = 1.3, ysize = 0.28, color = mypalette_earbones[4])
allometry_perioticL_log_int_taxa_plot

# ggtitle ("Log-transformed BZW vs periotic length by group - p-value < 0.001***")+  #copy from model summary

#Add confidence intervals
#Create data for confidence intervals
perioticW_log_newX_taxa <- expand.grid(BZW_log = seq(from = min(periotic_meas_taxa$BZW_log), to = max(periotic_meas_taxa$BZW_log), length.out = 30), #use min and max of x values (BZW) as limits and use number of specimens as length of sequence
                                       taxon = c("B.bonaerensis", "Ph.phocoena", "St.attenuata"))         #warp x_vals on values of x axis (BZW), add groupings
perioticW_log_newY_taxa <- predict(allometry_perioticW_log_int_taxa, newdata = perioticW_log_newX_taxa, interval="confidence",
                                   level = 0.95)
#Make data frame of data for confidence intervals
perioticW_log_conf_intervals_taxa <- data.frame(perioticW_log_newX_taxa, perioticW_log_newY_taxa, groups_taxa_p)
#Rename columns to match main plot tibble variables for x and y
perioticW_log_conf_intervals_taxa <- rename(perioticW_log_conf_intervals_taxa, perioticW_log = fit)
perioticW_log_conf_intervals_taxa

#Plot
allometry_perioticW_log_int_taxa_plot <- ggplot(periotic_meas_taxa, aes(y = perioticW_log, x = BZW_log, fill = taxon, color = taxon, linetype = group)) +
  geom_smooth(data = perioticW_log_conf_intervals_taxa, aes(ymin = lwr, ymax = upr, fill = taxon, colour = taxon, linetype = group), stat = 'identity',          #confidence intervals and reg line, before points
              size = 0.8, alpha = 0.3, show.legend = F)+      #put col and other graphics OUTSIDE of aes()!!!
  geom_point(size = 3, alpha = 0.2)+       #points after, so they are on top
  scale_color_manual(name = "Taxa", labels  = c("B.bonaerensis", "Ph.phocoena", "St.attenuata"), values = c(mypalette_earbones[1], mypalette_earbones[3:4]), #select colors from palette from taxa
                     aesthetics = c("color","fill"))+         
  theme_classic(base_size = 12)+
  xlab("Log(BZW)")+
  ylab("Log(Periotic width) ***")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11), legend.position = "none", legend.direction = "vertical", 
        axis.title.x = element_text(vjust = -1), 
        axis.title.y = element_text(face = "bold", vjust = 2, size = 14))
allometry_perioticW_log_int_taxa_plot <- move_layers(allometry_perioticW_log_int_taxa_plot, "GeomPoint", position = "top")
allometry_perioticW_log_int_taxa_plot <- allometry_perioticW_log_int_taxa_plot + 
  add_phylopic(B.bonaerensis, alpha = 1, x = 4.8, y = 3.5, ysize = 0.17, color = mypalette_earbones[1])+
  add_phylopic(Ph.phocoena, alpha = 1, x = 5.65, y = 2.7, ysize = 0.21, color = mypalette_earbones[3])+
  add_phylopic(St.attenuata, alpha = 1, x = 5.65, y = 3.2, ysize = 0.17, color = mypalette_earbones[4])
allometry_perioticW_log_int_taxa_plot

#ggtitle ("Log-transformed BZW vs periotic width by group - p-value < 0.001***")+  #copy from model summary

#All best models plots together
grid.arrange(allometry_bullaL_log_int_taxa_plot,  allometry_perioticL_log_int_taxa_plot,
             allometry_bullaW_log_int_taxa_plot, allometry_perioticW_log_int_taxa_plot)


#USE LOG VALUES - ALLOMETRY OF BULLA AND PERIOTIC MEASURAMENTS FOR EACH GROUP - NATURAL VALUES (2-way ANOVA) ----
##Regression of natural values with group effect - check differences between groups
allometry_BZW_bullaL_gp <- lm(bullaL ~ BZW * group, data  = bulla_meas)
allometry_BZW_bullaW_gp <- lm(bullaW ~ BZW * group, data  = bulla_meas)
allometry_BZW_perioticL_gp <- lm(perioticL ~ BZW * group, data  = periotic_meas)
allometry_BZW_perioticW_gp <- lm(perioticW ~ BZW * group, data  = periotic_meas)

#Check results
summary(allometry_BZW_bullaL_gp)
anova(allometry_BZW_bullaL_gp)
summary(allometry_BZW_bullaW_gp)
anova(allometry_BZW_bullaW_gp)
summary(allometry_BZW_perioticL_gp)
anova(allometry_BZW_perioticL_gp)
summary(allometry_BZW_perioticW_gp)
anova(allometry_BZW_perioticW_gp)

#Plot diagnostics
autoplot(allometry_BZW_bullaL_gp, smooth.colour = NA)
autoplot(allometry_BZW_bullaW_gp, smooth.colour = NA)
autoplot(allometry_BZW_perioticL_gp, smooth.colour = NA)
autoplot(allometry_BZW_perioticW_gp, smooth.colour = NA)

#Save regressions to file
sink("Output/allometry_BZW_LW_gp-int.txt", append = F)
print("Bulla length")
print(summary(allometry_BZW_bullaL_gp))
print(anova(allometry_BZW_bullaL_gp))

print("Bulla width")
print(summary(allometry_BZW_bullaW_gp))
print(anova(allometry_BZW_bullaW_gp))

print("Periotic length")
print(summary(allometry_BZW_perioticL_gp))
print(anova(allometry_BZW_perioticL_gp))

print("Periotic width")
print(summary(allometry_BZW_perioticW_gp))
print(anova(allometry_BZW_perioticW_gp))
sink()

#Plot regression lines per group
#Add confidence intervals
#Create data for confidence intervals
bullaL_newX <- expand.grid(BZW = seq(from = min(bulla_meas$BZW), to = max(bulla_meas$BZW), length.out = 30), #use min and max of x values (BZW) as limits and use number of specimens as length of sequence
                           group = c("Mysticeti", "Odontoceti"))         #warp x_vals on values of x axis (BZW), add groupings
bullaL_newY <- predict(allometry_BZW_bullaL_gp, newdata = bullaL_newX, interval="confidence",
                       level = 0.95)
#Make data frame of data for confidence intervals
bullaL_conf_intervals <- data.frame(bullaL_newX, bullaL_newY)
#Rename columns to match main plot tibble variables for x and y
bullaL_conf_intervals <- rename(bullaL_conf_intervals, bullaL = fit)
bullaL_conf_intervals

#Plot
allometry_BZW_bullaL_gp_plot <- ggplot(bulla_meas, aes(y = bullaL, x = BZW, fill = group, color = group)) +
  geom_smooth(data = bullaL_conf_intervals, aes(ymin = lwr, ymax = upr, fill = group, colour = group, linetype = group), stat = 'identity',          #confidence intervals and reg line, before points
              size = 0.8, alpha = 0.2, show.legend = F)+      #put col and other graphics OUTSIDE of aes()!!!
  geom_point(size = 3)+       #points after, so they are on top
  scale_color_manual(name = "Groups", labels  = c("Mysticeti", "Odontoceti"), values = mypalette_earbones[1:2], aesthetics = c("color","fill"))+         
  theme_classic(base_size = 12)+
  xlab("BZW (mm)")+
  ylab("Bulla length (mm)")+
  ggtitle ("BZW vs Bulla length by group - p-value < 0.001***")+  #copy from model summary
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12), legend.position = "none", legend.direction = "vertical", 
        axis.title.x = element_text(vjust = -1), 
        axis.title.y = element_text(vjust = 2))
allometry_BZW_bullaL_gp_plot <- move_layers(allometry_BZW_bullaL_gp_plot, "GeomPoint", position = "top")
allometry_BZW_bullaL_gp_plot <- allometry_BZW_bullaL_gp_plot + 
  add_phylopic(myst, alpha = 1, x = 100, y = 60, ysize = 30, color = mypalette_earbones[1])+
  add_phylopic(odont, alpha = 1, x = 300, y = 20, ysize = 25, color = mypalette_earbones[2])

#Add confidence intervals
#Create data for confidence intervals
bullaW_newX <- expand.grid(BZW = seq(from = min(bulla_meas$BZW), to = max(bulla_meas$BZW), length.out = 30), #use min and max of x values (BZW) as limits and use number of specimens as length of sequence
                           group = c("Mysticeti", "Odontoceti"))         #warp x_vals on values of x axis (BZW), add groupings
bullaW_newY <- predict(allometry_BZW_bullaW_gp, newdata = bullaW_newX, interval="confidence",
                       level = 0.95)
#Make data frame of data for confidence intervals
bullaW_conf_intervals <- data.frame(bullaW_newX, bullaW_newY)
#Rename columns to match main plot tibble variables for x and y
bullaW_conf_intervals <- rename(bullaW_conf_intervals, bullaW = fit)
bullaW_conf_intervals

#Plot
allometry_BZW_bullaW_gp_plot <- ggplot(bulla_meas, aes(y = bullaW, x = BZW, fill = group, color = group)) +
  geom_smooth(data = bullaW_conf_intervals, aes(ymin = lwr, ymax = upr, fill = group, colour = group, linetype = group), stat = 'identity',          #confidence intervals and reg line, before points
              size = 0.8, alpha = 0.2, show.legend = F)+      #put col and other graphics OUTSIDE of aes()!!!
  geom_point(size = 3)+       #points after, so they are on top
  scale_color_manual(name = "Groups", labels  = c("Mysticeti", "Odontoceti"), values = mypalette_earbones[1:2], aesthetics = c("color","fill"))+         
  theme_classic(base_size = 12)+
  xlab("BZW (mm)")+
  ylab("Bulla width (mm)")+
  ggtitle ("BZW vs Bulla width by group - p-value < 0.001***")+  #copy from model summary
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12), legend.position = "none", legend.direction = "vertical", 
        axis.title.x = element_text(vjust = -1), 
        axis.title.y = element_text(vjust = 2))
allometry_BZW_bullaW_gp_plot <- move_layers(allometry_BZW_bullaW_gp_plot, "GeomPoint", position = "top")
allometry_BZW_bullaW_gp_plot + 
  add_phylopic(myst, alpha = 1, x = 100, y = 60, ysize = 30, color = mypalette_earbones[1])+
  add_phylopic(odont, alpha = 1, x = 300, y = 20, ysize = 25, color = mypalette_earbones[2])

#Add confidence intervals
#Create data for confidence intervals
perioticL_newX <- expand.grid(BZW = seq(from = min(periotic_meas$BZW), to = max(periotic_meas$BZW), length.out = 21), #use min and max of x values (BZW) as limits and use number of specimens as length of sequence
                              group = c("Mysticeti", "Odontoceti"))         #warp x_vals on values of x axis (BZW), add groupings
perioticL_newY <- predict(allometry_BZW_perioticL_gp, newdata = perioticL_newX, interval="confidence",
                          level = 0.95)
#Make data frame of data for confidence intervals
perioticL_conf_intervals <- data.frame(perioticL_newX, perioticL_newY)
#Rename columns to match main plot tibble variables for x and y
perioticL_conf_intervals <- rename(perioticL_conf_intervals, perioticL = fit)
perioticL_conf_intervals

#Plot
allometry_BZW_perioticL_gp_plot <- ggplot(periotic_meas, aes(y = perioticL, x = BZW, fill = group, color = group)) +
  geom_smooth(data = perioticL_conf_intervals, aes(ymin = lwr, ymax = upr, fill = group, colour = group, linetype = group), stat = 'identity',          #confidence intervals and reg line, before points
              size = 0.8, alpha = 0.2, show.legend = F)+      #put col and other graphics OUTSIDE of aes()!!!
  geom_point(size = 3)+       #points after, so they are on top
  scale_color_manual(name = "Groups", labels  = c("Mysticeti", "Odontoceti"), values = mypalette_earbones[1:2], aesthetics = c("color","fill"))+         
  theme_classic(base_size = 12)+
  xlab("BZW (mm)")+
  ylab("Periotic length (mm)")+
  ggtitle ("BZW vs Periotic length by group - p-value < 0.001***")+  #copy from model summary
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12), legend.position = "none", legend.direction = "vertical", 
        axis.title.x = element_text(vjust = -1), 
        axis.title.y = element_text(vjust = 2))
allometry_BZW_perioticL_gp_plot <- move_layers(allometry_BZW_perioticL_gp_plot, "GeomPoint", position = "top")
allometry_BZW_perioticL_gp_plot + 
  add_phylopic(myst, alpha = 1, x = 100, y = 60, ysize = 30, color = mypalette_earbones[1])+
  add_phylopic(odont, alpha = 1, x = 300, y = 20, ysize = 25, color = mypalette_earbones[2])

#Add confidence intervals
#Create data for confidence intervals
perioticW_newX <- expand.grid(BZW = seq(from = min(periotic_meas$BZW), to = max(periotic_meas$BZW), length.out = 21), #use min and max of x values (BZW) as limits and use number of specimens as length of sequence
                              group = c("Mysticeti", "Odontoceti"))         #warp x_vals on values of x axis (BZW), add groupings
perioticW_newY <- predict(allometry_BZW_perioticW_gp, newdata = perioticW_newX, interval="confidence",
                          level = 0.95)
#Make data frame of data for confidence intervals
perioticW_conf_intervals <- data.frame(perioticW_newX, perioticW_newY)
#Rename columns to match main plot tibble variables for x and y
perioticW_conf_intervals <- rename(perioticW_conf_intervals, perioticW = fit)
perioticW_conf_intervals

#Plot
allometry_BZW_perioticW_gp_plot <- ggplot(periotic_meas, aes(y = perioticW, x = BZW, fill = group, color = group)) +
  geom_smooth(data = perioticW_conf_intervals, aes(ymin = lwr, ymax = upr, fill = group, colour = group, linetype = group), stat = 'identity',          #confidence intervals and reg line, before points
              size = 0.8, alpha = 0.2, show.legend = F)+      #put col and other graphics OUTSIDE of aes()!!!
  geom_point(size = 3)+       #points after, so they are on top
  scale_color_manual(name = "Groups", labels  = c("Mysticeti", "Odontoceti"), values = mypalette_earbones[1:2], aesthetics = c("color","fill"))+         
  theme_classic(base_size = 12)+
  xlab("BZW (mm)")+
  ylab("Periotic width (mm)")+
  ggtitle ("BZW vs Periotic width by group - p-value < 0.001*** (group p < 0.01 **)")+  #copy from model summary
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11), legend.position = "none", legend.direction = "vertical", 
        axis.title.x = element_text(vjust = -1), 
        axis.title.y = element_text(vjust = 2))
allometry_BZW_perioticW_gp_plot <- move_layers(allometry_BZW_perioticW_gp_plot, "GeomPoint", position = "top")
allometry_BZW_perioticW_gp_plot <- allometry_BZW_perioticW_gp_plot + 
  add_phylopic(myst, alpha = 1, x = 100, y = 40, ysize = 30, color = mypalette_earbones[1])+
  add_phylopic(odont, alpha = 1, x = 300, y = 18, ysize = 25, color = mypalette_earbones[2])

#Check if groups can be different only in size (intercept) and not slope
allometry_BZW_bullaL_gp2 <- lm(bullaL ~ BZW + group, data  = bulla_meas)
allometry_BZW_bullaW_gp2 <- lm(bullaW ~ BZW + group, data  = bulla_meas)
allometry_BZW_perioticL_gp2 <- lm(perioticL ~ BZW + group, data  = periotic_meas)
allometry_BZW_perioticW_gp2 <- lm(perioticW ~ BZW + group, data  = periotic_meas)

#Check results
summary(allometry_BZW_bullaL_gp2)
anova(allometry_BZW_bullaL_gp2)
summary(allometry_BZW_bullaW_gp2)
anova(allometry_BZW_bullaW_gp2)
summary(allometry_BZW_perioticL_gp2)
anova(allometry_BZW_perioticL_gp2)
summary(allometry_BZW_perioticW_gp2)
anova(allometry_BZW_perioticW_gp2)

#Save regressions to file
sink("Output/allometry_BZW_LW_gp+.txt", append = F)
print("Bulla length")
print(summary(allometry_BZW_bullaL_gp2))
print(anova(allometry_BZW_bullaL_gp2))

print("Bulla width")
print(summary(allometry_BZW_bullaW_gp2))
print(anova(allometry_BZW_bullaW_gp2))

print("Periotic length")
print(summary(allometry_BZW_perioticL_gp2))
print(anova(allometry_BZW_perioticL_gp2))

print("Periotic width")
print(summary(allometry_BZW_perioticW_gp2))
print(anova(allometry_BZW_perioticW_gp2))
sink()



#USE 2-way ANOVA - ALLOMETRY OF BULLA AND PERIOTIC MEASURAMENTS FOR EACH TAXON (LM) ----
##Regression of natural values by group
#Perform LM analysis to get p-values and coefficients of slope/intercept
allometry_BZW_bullaL_mysticeti <- lm(bullaL ~ BZW, data  = bulla_meas_mysticeti)
allometry_BZW_bullaL_odontoceti <- lm(bullaL ~ BZW, data = bulla_meas_odontoceti)

allometry_BZW_bullaW_mysticeti <- lm(bullaW ~ BZW, data = bulla_meas_mysticeti)
allometry_BZW_bullaW_odontoceti <- lm(bullaW ~ BZW, data = bulla_meas_odontoceti)

allometry_BZW_perioticL_mysticeti <- lm(perioticL ~ BZW, data = periotic_meas_mysticeti)
allometry_BZW_perioticL_odontoceti <- lm(perioticL ~ BZW, data  = periotic_meas_odontoceti)

allometry_BZW_perioticW_mysticeti <- lm(perioticW ~ BZW, data  = periotic_meas_mysticeti)
allometry_BZW_perioticW_odontoceti <- lm(perioticW ~ BZW, data  = periotic_meas_odontoceti)

#Check results
summary(allometry_BZW_bullaL_mysticeti)
summary(allometry_BZW_bullaL_odontoceti)

summary(allometry_BZW_bullaW_mysticeti)
summary(allometry_BZW_bullaW_odontoceti)

summary(allometry_BZW_perioticL_mysticeti)
summary(allometry_BZW_perioticL_odontoceti)

summary(allometry_BZW_perioticW_mysticeti)
summary(allometry_BZW_perioticW_odontoceti)

#Save regressions to file
sink("Output/allometry_BZW_LW_myst+odont.txt", append = F)
print("Bulla length")
print(summary(allometry_BZW_bullaL_mysticeti))
print(summary(allometry_BZW_bullaL_odontoceti))

print("Bulla width")
print(summary(allometry_BZW_bullaW_mysticeti))
print(summary(allometry_BZW_bullaW_odontoceti))

print("Periotic length")
print(summary(allometry_BZW_perioticL_mysticeti))
print(summary(allometry_BZW_perioticL_odontoceti))

print("Periotic width")
print(summary(allometry_BZW_perioticW_mysticeti))
print(summary(allometry_BZW_perioticW_odontoceti))
sink()

#Plot diagnostic plots
autoplot(allometry_BZW_bullaL_mysticeti, smooth.colour = NA)
autoplot(allometry_BZW_bullaL_odontoceti, smooth.colour = NA)

#Plot regressions for each bone-meas
allometry_BZW_bullaL_plot <- ggplot(bulla_meas, aes(y = bullaL, x = BZW, fill = group, color = group)) + #xend useful to make sure graphs goes close to 0
  geom_point(size = 3)+ 
  stat_smooth(data = bulla_meas_mysticeti, aes(y = bullaL, x = BZW), method = lm, formula = y ~ x, se = T, 
              inherit.aes = F, fullrange = T, fill = mypalette_earbones[1], alpha = 0.2,
              linetype = 1, colour = mypalette_earbones[1])+
  stat_smooth(data = bulla_meas_odontoceti, aes(y = bullaL, x = BZW), method = lm, formula = y ~ x, se = T, 
              inherit.aes = F, fullrange = T, fill = mypalette_earbones[2], alpha = 0.2,
              linetype = 2, colour = mypalette_earbones[2])+ 
  scale_color_manual(name = "Groups", labels  = c("Mysticeti", "Odontoceti"), values = mypalette_earbones[1:2], aesthetics = c("color","fill"))+
  theme_bw()+
  xlab("BZW (mm)")+
  ylab("Bulla length (mm)")+
  theme(legend.direction = "vertical", axis.title.x = element_text(vjust = -1), axis.title.y = element_text(vjust = 2))
allometry_BZW_bullaL_plot <- move_layers(allometry_BZW_bullaL, "GeomPoint", position = "top")
allometry_BZW_bullaL_plot

allometry_BZW_bullaW_plot <- ggplot(bulla_meas, aes(y = bullaW, x = BZW, fill = group, color = group)) + #xend useful to make sure graphs goes close to 0
  geom_point(size = 3)+ 
  stat_smooth(data = bulla_meas_mysticeti, aes(y = bullaW, x = BZW), method = lm, formula = y ~ x, se = T, 
              inherit.aes = F, fullrange = T, fill = mypalette_earbones[1], alpha = 0.2,
              linetype = 1, colour = mypalette_earbones[1])+
  stat_smooth(data = bulla_meas_odontoceti, aes(y = bullaW, x = BZW), method = lm, formula = y ~ x, se = T, 
              inherit.aes = F, fullrange = T, fill = mypalette_earbones[2], alpha = 0.2,
              linetype = 2, colour = mypalette_earbones[2])+ 
  scale_color_manual(name = "Groups", labels  = c("Mysticeti", "Odontoceti"), values = mypalette_earbones[1:2], aesthetics = c("color","fill"))+
  theme_bw()+
  xlab("BZW (mm)")+
  ylab("Bulla width (mm)")+
  theme(legend.direction = "vertical", axis.title.x = element_text(vjust = -1), axis.title.y = element_text(vjust = 2))
allometry_BZW_bullaW_plot <- move_layers(allometry_BZW_bullaW, "GeomPoint", position = "top")
allometry_BZW_bullaW_plot

allometry_BZW_perioticL_plot <- ggplot(periotic_meas, aes(y = perioticL, x = BZW, fill = group, color = group)) + #xend useful to make sure graphs goes close to 0
  geom_point(size = 3)+ 
  stat_smooth(data = periotic_meas_mysticeti, aes(y = perioticL, x = BZW), method = lm, formula = y ~ x, se = T, 
              inherit.aes = F, fullrange = T, fill = mypalette_earbones[1], alpha = 0.2,
              linetype = 1, colour = mypalette_earbones[1])+
  stat_smooth(data = periotic_meas_odontoceti, aes(y = perioticL, x = BZW), method = lm, formula = y ~ x, se = T, 
              inherit.aes = F, fullrange = T, fill = mypalette_earbones[2], alpha = 0.2,
              linetype = 2, colour = mypalette_earbones[2])+ 
  scale_color_manual(name = "Groups", labels  = c("Mysticeti", "Odontoceti"), values = mypalette_earbones[1:2], aesthetics = c("color","fill"))+
  theme_bw()+
  xlab("BZW (mm)")+
  ylab("Periotic length (mm)")+
  theme(legend.direction = "vertical", axis.title.x = element_text(vjust = -1), axis.title.y = element_text(vjust = 2))
allometry_BZW_perioticL_plot <- move_layers(allometry_BZW_perioticL, "GeomPoint", position = "top")
allometry_BZW_perioticL_plot

allometry_BZW_perioticW_plot <- ggplot(periotic_meas, aes(y = perioticW, x = BZW, fill = group, color = group)) + #xend useful to make sure graphs goes close to 0
  geom_point(size = 3)+ 
  stat_smooth(data = periotic_meas_mysticeti, aes(y = perioticW, x = BZW), method = lm, formula = y ~ x, se = T, 
              inherit.aes = F, fullrange = T, fill = mypalette_earbones[1], alpha = 0.2,
              linetype = 1, colour = mypalette_earbones[1])+
  stat_smooth(data = periotic_meas_odontoceti, aes(y = perioticW, x = BZW), method = lm, formula = y ~ x, se = T, 
              inherit.aes = F, fullrange = T, fill = mypalette_earbones[2], alpha = 0.2,
              linetype = 2, colour = mypalette_earbones[2])+ 
  scale_color_manual(name = "Groups", labels  = c("Mysticeti", "Odontoceti"), values = mypalette_earbones[1:2], aesthetics = c("color","fill"))+
  theme_bw()+
  xlab("BZW (mm)")+
  ylab("Periotic width (mm)")+
  theme(legend.direction = "vertical", axis.title.x = element_text(vjust = -1), axis.title.y = element_text(vjust = 2))
allometry_BZW_perioticW_plot <- move_layers(allometry_BZW_perioticW, "GeomPoint", position = "top")
allometry_BZW_perioticW_plot

##Regression of log values by group
#Perform LM analysis to get p-values and coefficients of slope/intercept
allometry_BZW_bullaL_log_mysticeti <- lm(bullaL_log ~ BZW_log, data  = bulla_meas_mysticeti)
allometry_BZW_bullaL_log_odontoceti <- lm(bullaL_log ~ BZW_log, data = bulla_meas_odontoceti)

allometry_BZW_bullaW_log_mysticeti <- lm(bullaW_log ~ BZW_log, data = bulla_meas_mysticeti)
allometry_BZW_bullaW_log_odontoceti <- lm(bullaW_log ~ BZW_log, data = bulla_meas_odontoceti)

allometry_BZW_perioticL_log_mysticeti <- lm(perioticL_log ~ BZW_log, data = periotic_meas_mysticeti)
allometry_BZW_perioticL_log_odontoceti <- lm(perioticL_log ~ BZW_log, data  = periotic_meas_odontoceti)

allometry_BZW_perioticW_log_mysticeti <- lm(perioticW_log ~ BZW_log, data  = periotic_meas_mysticeti)
allometry_BZW_perioticW_log_odontoceti <- lm(perioticW_log ~ BZW_log, data  = periotic_meas_odontoceti)

#Check results
summary(allometry_BZW_bullaL_log_mysticeti)
summary(allometry_BZW_bullaL_log_odontoceti)

summary(allometry_BZW_bullaW_log_mysticeti)
summary(allometry_BZW_bullaW_log_odontoceti)

summary(allometry_BZW_perioticL_log_mysticeti)
summary(allometry_BZW_perioticL_log_odontoceti)

summary(allometry_BZW_perioticW_log_mysticeti)
summary(allometry_BZW_perioticW_log_odontoceti)

#Save regressions to file
sink("Output/allometry_BZW_LW_log_myst+odont.txt", append = F)
print("Bulla length")
print(summary(allometry_BZW_bullaL_log_mysticeti))
print(summary(allometry_BZW_bullaL_log_odontoceti))

print("Bulla width")
print(summary(allometry_BZW_bullaW_log_mysticeti))
print(summary(allometry_BZW_bullaW_log_odontoceti))

print("Periotic length")
print(summary(allometry_BZW_perioticL_log_mysticeti))
print(summary(allometry_BZW_perioticL_log_odontoceti))

print("Periotic width")
print(summary(allometry_BZW_perioticW_log_mysticeti))
print(summary(allometry_BZW_perioticW_log_odontoceti))
sink()

#Plot diagnostic plots
autoplot(allometry_BZW_bullaL_log_mysticeti, smooth.colour = NA)
autoplot(allometry_BZW_bullaL_log_odontoceti, smooth.colour = NA)

#Plot regressions for each bone-meas - no SE
allometry_BZW_bullaL_log_plot <- ggplot(bulla_meas, aes(y = bullaL_log, x = BZW_log, fill = group, color = group)) + #xend useful to make sure graphs goes close to 0
  geom_point(size = 3)+ 
  stat_smooth(data = bulla_meas_mysticeti, aes(y = bullaL_log, x = BZW_log), method = lm, formula = y ~ x, se = F, 
              inherit.aes = F, fullrange = T,
              linetype = 1, colour = mypalette_earbones[1])+
  stat_smooth(data = bulla_meas_odontoceti, aes(y = bullaL_log, x = BZW_log), method = lm, formula = y ~ x, se = F, 
              inherit.aes = F, fullrange = T,
              linetype = 2, colour = mypalette_earbones[2])+ 
  scale_color_manual(name = "Groups", labels  = c("Mysticeti", "Odontoceti"), values = mypalette_earbones[1:2], aesthetics = c("color","fill"))+
  theme_bw()+
  xlab("Log(BZW)")+
  ylab("Log(Bulla length)")+
  theme(legend.direction = "vertical", axis.title.x = element_text(vjust = -1), axis.title.y = element_text(vjust = 2))
allometry_BZW_bullaL_log_plot <- move_layers(allometry_BZW_bullaL_log, "GeomPoint", position = "top")
allometry_BZW_bullaL_log_plot

allometry_BZW_bullaW_log_plot <- ggplot(bulla_meas, aes(y = bullaW_log, x = BZW_log, fill = group, color = group)) + #xend useful to make sure graphs goes close to 0
  geom_point(size = 3)+ 
  stat_smooth(data = bulla_meas_mysticeti, aes(y = bullaL_log, x = BZW_log), method = lm, formula = y ~ x, se = F, 
              inherit.aes = F, fullrange = T,
              linetype = 1, colour = mypalette_earbones[1])+
  stat_smooth(data = bulla_meas_odontoceti, aes(y = bullaL_log, x = BZW_log), method = lm, formula = y ~ x, se = F, 
              inherit.aes = F, fullrange = T,
              linetype = 2, colour = mypalette_earbones[2])+ 
  scale_color_manual(name = "Groups", labels  = c("Mysticeti", "Odontoceti"), values = mypalette_earbones[1:2], aesthetics = c("color","fill"))+
  theme_bw()+
  xlab("Log(BZW)")+
  ylab("Log(Bulla width)")+
  theme(legend.direction = "vertical", axis.title.x = element_text(vjust = -1), axis.title.y = element_text(vjust = 2))
allometry_BZW_bullaW_log_plot <- move_layers(allometry_BZW_bullaW_log, "GeomPoint", position = "top")
allometry_BZW_bullaW_log_plot

allometry_BZW_perioticL_log <- ggplot(periotic_meas, aes(y = perioticL_log, x = BZW_log, fill = group, color = group)) + #xend useful to make sure graphs goes close to 0
  geom_point(size = 3)+ 
  stat_smooth(data = bulla_meas_mysticeti, aes(y = bullaL_log, x = BZW_log), method = lm, formula = y ~ x, se = F, 
              inherit.aes = F, fullrange = T,
              linetype = 1, colour = mypalette_earbones[1])+
  stat_smooth(data = bulla_meas_odontoceti, aes(y = bullaL_log, x = BZW_log), method = lm, formula = y ~ x, se = F, 
              inherit.aes = F, fullrange = T,
              linetype = 2, colour = mypalette_earbones[2])+ 
  scale_color_manual(name = "Groups", labels  = c("Mysticeti", "Odontoceti"), values = mypalette_earbones[1:2], aesthetics = c("color","fill"))+
  theme_bw()+
  xlab("Log(BZW)")+
  ylab("Log(Periotic length)")+
  theme(legend.direction = "vertical", axis.title.x = element_text(vjust = -1), axis.title.y = element_text(vjust = 2))
allometry_BZW_perioticL_log_plot <- move_layers(allometry_BZW_perioticL_log, "GeomPoint", position = "top")
allometry_BZW_perioticL_log_plot

allometry_BZW_perioticW_log_plot <- ggplot(periotic_meas, aes(y = perioticW_log, x = BZW_log, fill = group, color = group)) + #xend useful to make sure graphs goes close to 0
  geom_point(size = 3)+ 
  stat_smooth(data = bulla_meas_mysticeti, aes(y = bullaL_log, x = BZW_log), method = lm, formula = y ~ x, se = F, 
              inherit.aes = F, fullrange = T,
              linetype = 1, colour = mypalette_earbones[1])+
  stat_smooth(data = bulla_meas_odontoceti, aes(y = bullaL_log, x = BZW_log), method = lm, formula = y ~ x, se = F, 
              inherit.aes = F, fullrange = T,
              linetype = 2, colour = mypalette_earbones[2])+ 
  scale_color_manual(name = "Groups", labels  = c("Mysticeti", "Odontoceti"), values = mypalette_earbones[1:2], aesthetics = c("color","fill"))+
  theme_bw()+
  xlab("Log(BZW)")+
  ylab("Log(Periotic width)")+
  theme(legend.direction = "vertical", axis.title.x = element_text(vjust = -1), axis.title.y = element_text(vjust = 2))
allometry_BZW_perioticW_log_plot <- move_layers(allometry_BZW_perioticW_log, "GeomPoint", position = "top")
allometry_BZW_perioticW_log_plot
