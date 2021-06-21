## GPSA DATA ANALYSIS ##

#CH.1 - LOAD LIBRARIES ----
#always do this first!!
library(geomorph) 
library(geiger)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(ggrepel)
library(gginnards)
library(ggphylomorpho)
library(ggfortify)
library(borealis)
library(ggthemes)
library(ggpubr)
library(ggplotify)
library(Morpho)

#CH. 7 - PCA ALLOMETRY RESIDUALS ----

#Run PCA on complete dataset
PCA_residuals <- gm.prcomp(allometry_residuals)

#List of PC components and proportion of variations
PCA_residuals

##View plot
PCA_res_plot <- plot(PCA_residuals, main = "PCA residuals",  pch = 21, #title and type of point to be used
                     col = "deeppink", bg = "deeppink", cex = 1, font.main = 2)  #improve graphics
#Add quick labels to plot
text(x = PCA_residuals$x[,1], y = PCA_residuals$x[,2], labels = rownames(PCA_residuals$x), 
     pos = 1, offset = 0.5, cex = 0.75)    #improve graphics

#Save PC scores as object to use later
pcscores_res <- PCA_residuals$x

##Make better plot using ggplot
#Read PC scores as tibble
pcscores_res <- as_tibble(pcscores_res)
#Add labels and other attributes to tibble as columns
pcscores_res <- pcscores_res %>% mutate(specimens = gdf$code, stage = gdf$stage, group = gdf$group, category = gdf$category)
glimpse(pcscores_res)

#Nice PCA plot with stages and groups
PCA_res_ggplot <- ggplot(pcscores_res, aes(x = Comp1, y = Comp2, label = specimens, colour = stage, shape = group))+
  geom_point(size = 3)+
  geom_text_repel(colour = "black", size = 3.5, max.overlaps = 40)+
  scale_colour_manual(name = "Growth stage", labels =  c("Early Fetus", "Late Fetus", "Neonate", "Juvenile", "Adult"), #to be ordered as they appear in tibble
                      values = mypalette_category)+            #legend and color adjustments
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("PC 1 (16.7%)")+ #copy this from standard PCA plot (PCA_res_plot)
  ylab("PC 2 (11.8%)")+
  ggtitle("PCA res data")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_res_ggplot

#Make hulls for PCA plot with hulls around growth stages
hulls_res <- pcscores_res %>%
  group_by(stage) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)
glimpse(hulls_res)

#Nice PCA plot with hulls around stages
PCA_res_hulls_ggplot <- ggplot(pcscores_res, aes(x = Comp1, y = Comp2, label = specimens, colour = stage, fill = stage))+
  geom_point(size = 3, aes(shape = group))+
  scale_colour_manual(name = "Growth stage", labels =  c("Early Fetus", "Late Fetus", "Neonate", "Juvenile", "Adult"), #to be ordered as they appear in tibble
                      values = mypalette_category)+            #legend and color adjustments
  geom_polygon(data = hulls_res, aes(x = x, y = y, fill = stage), alpha = .5, show.legend = FALSE)+ #colored hulls with transparency
  scale_fill_manual(name = "Growth stage", labels = c("Early Fetus", "Late Fetus", "Neonate", "Juvenile", "Adult"),
                    values =  mypalette_category)+ #must match scale_colour_manual
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("PC 1 (16.7%)")+ #copy this from standard PCA plot (PCA_res_plot)
  ylab("PC 2 (11.8%)")+
  ggtitle("PCA res data")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
  geom_text_repel(colour = "black", size = 3.5, max.overlaps = 40)

#Visualize plot and save as PDF using menu in bar on the right
PCA_res_hulls_ggplot  

#Make hulls for PCA plot with hulls around categories (3 stages: early fetus, late fetus, postnatal)
hulls_res_category <- pcscores_res %>%
  group_by(category) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)
glimpse(hulls_res_category)

#Nice PCA plot with hulls around categories (3 stages: early fetus, late fetus, postnatal)
PCA_res_category_ggplot <- ggplot(pcscores_res, aes(x = Comp1, y = Comp2, label = specimens, colour = category, fill = category))+
  geom_point(size = 3, aes(shape = group))+
  scale_colour_manual(name = "Growth stage", labels =  c("Early Fetus", "Late Fetus", "Postnatal"), #to be ordered as they appear in tibble
                      values = c(mypalette_category[1:2], mypalette_category[4]))+            #legend and color adjustments
  geom_polygon(data = hulls_res_category, aes(x = x, y = y, fill = category), alpha = .5, show.legend = FALSE)+ #colored hulls with transparency
  scale_fill_manual(name = "Growth stage", labels = c("Early Fetus", "Late Fetus", "Postnatal"),
                    values =  c(mypalette_category[1:2], mypalette_category[4]))+ #must match scale_colour_manual
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("PC 1 (16.7%)")+ #copy this from standard PCA plot (PCA_res_plot)
  ylab("PC 2 (11.8%)")+
  ggtitle("PCA res data")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  geom_text_repel(colour = "black", size = 3.5, max.overlaps = 40)

#Visualize plot and save as PDF using menu in bar on the right
PCA_res_category_ggplot

##PCA plot for each stage
#Create one tibble for each stage
#First split between stages
pcscores_res_stages <- pcscores_res %>% group_by(category) %>% group_split()

#Save as separate tibbles - 1 for early fetus, 1 for late fetus, 1 for all 3 postnatal stages (neonate, juvenile, adult)
pcscores_res_earlyfetus <- pcscores_res_stages[[1]]
pcscores_res_latefetus <- pcscores_res_stages[[2]]
pcscores_res_postnatal <- pcscores_res_stages[[3]]

#Nice PCA plot with groups early fetus
PCA_res_earlyfetus_ggplot <- ggplot(pcscores_res_earlyfetus, aes(x = Comp1, y = Comp2, label = specimens, shape = group))+
  geom_point(size = 3, color = mypalette_category[1])+
  geom_text_repel(colour = "black", size = 3.5, max.overlaps = 40)+
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("PC 1 (16.7%)")+ #copy this from standard PCA plot (PCA_res_plot)
  ylab("PC 2 (11.8%)")+
  ggtitle("PCA early fetus")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_res_earlyfetus_ggplot

#Nice PCA plot with groups late fetus
PCA_res_latefetus_ggplot <- ggplot(pcscores_res_latefetus, aes(x = Comp1, y = Comp2, label = specimens, shape = group))+
  geom_point(size = 3, color = mypalette_category[2])+
  geom_text_repel(colour = "black", size = 3.5, max.overlaps = 40)+
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("PC 1 (16.7%)")+ #copy this from standard PCA plot (PCA_res_plot)
  ylab("PC 2 (11.8%)")+
  ggtitle("PCA late fetus")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_res_latefetus_ggplot

#Nice PCA plot with groups postnatal
PCA_res_postnatal_ggplot <- ggplot(pcscores_res_postnatal, aes(x = Comp1, y = Comp2, label = specimens, colour = stage, shape = group))+
  geom_point(size = 3)+
  geom_text_repel(colour = "black", size = 3.5, max.overlaps = 40)+
  scale_colour_manual(name = "Growth stage", labels =  c("Neonate", "Juvenile", "Adult"), #to be ordered as they appear in tibble
                      values = mypalette_category[3:5])+            #legend and color adjustments
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("PC 1 (16.7%)")+ #copy this from standard PCA plot (PCA_res_plot)
  ylab("PC 2 (11.8%)")+
  ggtitle("PCA postnatal")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_res_postnatal_ggplot

#CH. 8 - MORPHOLOGICAL DISPARITY -----

#Calculate Procrustes variances and distances between groups, with p-value for each pair of groups
#How different are each group shapes compared to other group shapes?
#Disparity between stages, considering entire dataset and not different groups
stage_disparity <- morphol.disparity(coords ~ 1, groups = ~ stage, iter = 999, data = gdf)

#Results and significance
summary(stage_disparity)

#Save results to file
sink("Output/stage_disparity.txt")
print(summary(stage_disparity))
sink() 

#Disparity between categories, considering entire dataset and not different groups
category_disparity <- morphol.disparity(coords ~ 1, groups = ~ category, iter = 999, data = gdf)

#Results and significance
summary(category_disparity)

#Save results to file
sink("Output/category_disparity.txt")
print(summary(category_disparity))
sink() 

#Disparity between categories, considering entire dataset AND groups
category_group_disparity <- morphol.disparity(coords ~ 1, groups = ~ category*group, iter = 999, data = gdf)

#Results and significance
summary(category_group_disparity)

#Save results to file
sink("Output/category_group_disparity.txt")
print(summary(category_group_disparity))
sink() 

#Disparity between groups for each stage - use previously defined dataframe for each stage (early fetus, late fetus, postnatal)
#Early fetus
group_disparity_earlyfetus <- morphol.disparity(coords ~ 1, groups = ~ group, iter = 999, data = gdf_earlyfetus)

#Results and significance
summary(group_disparity_earlyfetus)

#Save results to file
sink("Output/group_disparity_earlyfetus.txt")
summary(group_disparity_earlyfetus)
sink() 

#Late fetus
group_disparity_latefetus <- morphol.disparity(coords ~ 1, groups = ~ group, iter = 999, data = gdf_latefetus)

#Results and significance
summary(group_disparity_latefetus)

#Save results to file
sink("Output/group_disparity_latefetus.txt")
summary(group_disparity_latefetus)
sink() 

#Postnatal
group_disparity_postnatal <- morphol.disparity(coords ~ 1, groups = ~ group, iter = 999, data = gdf_postnatal)

#Results and significance
summary(group_disparity_postnatal)

#Save results to file
sink("Output/group_disparity_postnatal.txt")
summary(group_disparity_postnatal)
sink() 

#CH. 9 - TRAJECTORY ANALYSIS ----
#Shows trajectories of shape variation, possible to compare trajectories of growth bewteen odontocetes and mysticetes

##BEFORE doing this make sure to re-import the classifiers file and re-order it (problems with previous category names) and add it again to the gdf
#Import classifiers for analyses, make sure they are in the same order as coordinates
classifiers <- read_csv("Data/bulla_classifiers.csv")

##Order classifiers by stage, useful for plot legend
#Make factor for variable
classifiers$stage <- factor(classifiers$stage, 
                            levels = c("earlyFetus", "lateFetus", "neonate", "juvenile", "adult")) #copy from string printed with the code above
#Order
classifiers <- classifiers[order(classifiers$stage),]

##Make data frame for analyses in geomorph
gdf <- geomorph.data.frame(coords = shape_array, 
                           code = classifiers$code, group = classifiers$group, 
                           stage = classifiers$stage, category = classifiers$category,
                           bulla_log = classifiers$bullaL_log, periotic_log = classifiers$perioticL_log)


#First perform procD.lm to create linear model that describes what we are trying to test - shape changes at each stage (category) considering the 2 groups
fit_shape_stage_group <- procD.lm(coords ~ group * category, iter = 499, data = gdf, RRPP = F)

#Check that there is a significant correlation
summary(fit_shape_stage_group)

#Use fit to calculate trajectories
group_trajectory <- trajectory.analysis(fit_shape_stage_group, groups = gdf$group, traj.pts = gdf$category, 
                                        pca = TRUE, print.progress = TRUE) 

#View results
#Magnitude differences between trajectories, standard summary - are trajectories different?
summary(group_trajectory, show.trajectories = TRUE, attribute = "MD") 
#Trajectory correlations - are the trajectories of the 2 groups correlated?
summary(group_trajectory, show.trajectories = TRUE, attribute = "TC", angle.type = "deg")
#Trajectory shape differences - what is the distance between trajectories?
summary(group_trajectory, show.trajectories = TRUE, attribute = "SD") 

#Save results to file
sink("Output/group_trajectory.txt")
print("Magnitude difference (absolute difference between path distances)")
summary(group_trajectory, show.trajectories = TRUE, attribute = "MD") 
print("Correlations (angles) between trajectories")
summary(group_trajectory, show.trajectories = TRUE, attribute = "TC", angle.type = "deg")
print("Shape differences between trajectory vectors")
summary(group_trajectory, show.trajectories = TRUE, attribute = "SD") 
sink() 

#Plot results - PCA of fitted values
group_trajectory_plot <- plot(group_trajectory, main = "Trajectories of growth by group",  pch = c(21,22), #title and type of point to be used
                              col = c("darkgray","lightgray"), bg = c("darkgray","lightgray"), cex = 1, font.main = 2) #improve graphics
#Add line between groups
add.trajectories(group_trajectory_plot, 
                 traj.pch = c(21, 22), traj.col = 1, traj.lty = 1, traj.lwd = 1, traj.cex = 1.5, traj.bg = 1, 
                 start.bg = "green", end.bg = "red") #trajectory line graphics
#Add legend to see which trajectory belongs to each group
legend(x= -20, y = 10, legend = c("Mysticeti","Odontoceti"), pch =  c(21, 22), pt.bg = 1)


##Make better PCA plot using ggplot
#Read PC scores as tibble
group_trajectory_pcscores <- as_tibble(group_trajectory_plot[["pc.points"]])

#Add group names and other attributes to tibble as columns
group_trajectory_pcscores <- group_trajectory_pcscores %>% mutate(category = gdf$category, group = gdf$group)
glimpse(group_trajectory_pcscores)

#Calculate means of PC1 and PC2 at each category per group to draw trajectories
group_trajectory_pcscores_means <- group_trajectory_pcscores %>% group_by(group, category) %>%
  summarise_at(vars(PC1, PC2), list(mean = mean))              #function for means, both columns
glimpse(group_trajectory_pcscores_means)

#Rename columns so they are easier to use for plot
group_trajectory_pcscores_means <- group_trajectory_pcscores_means %>% rename(x = PC1_mean, y = PC2_mean)
group_trajectory_pcscores_means

#Nice plot
group_trajectory_ggplot <- ggplot(group_trajectory_pcscores, aes(x = PC1, y = PC2, colour = category, shape = group))+
  geom_point(size = 3, alpha = 0.4)+
  geom_point(data = group_trajectory_pcscores_means, aes(x = x, y = y, colour = category, shape = group), size = 5)+
  #add trajectory lines, one line for each, write row number from tibble, should be in order as legend of plot
  geom_segment(data = group_trajectory_pcscores_means, aes(x = x[1], y = y[1],  #earlyfetusMyst
                                                           xend =  x[2], yend = y[2], colour = age), #latefetusMyst
               colour = "gray17", size = 1, linejoin = 'mitre', linetype = 1, arrow = arrow(angle = 30, length = unit(0.03, "npc"), ends = "last", type = "closed"))+
  #add arrow at end  
  geom_segment(data = group_trajectory_pcscores_means, aes(x = x[2], y = y[2], #latefetusMyst
                                                           xend = x[3], yend = y[3], colour = age), #postnatalMyst
               colour = "gray17", size = 1, linejoin = 'mitre', linetype = 1, arrow = arrow(angle = 30, length = unit(0.03, "npc"), ends = "last", type = "closed"))+
  #add trajectory lines, one line for each, write row number from tibble, should be in order as legend of plot
  geom_segment(data = group_trajectory_pcscores_means, aes(x = x[4], y = y[4],  #earlyfetusOdont
                                                           xend =  x[5], yend = y[5], colour = age), #latefetusOdont
               colour = "gray17", size = 1, linejoin = 'mitre', linetype = 2, arrow = arrow(angle = 30, length = unit(0.03, "npc"), ends = "last", type = "closed"))+
  #add arrow at end  
  geom_segment(data = group_trajectory_pcscores_means, aes(x = x[5], y = y[5], #latefetusOdont
                                                           xend = x[6], yend = y[6], colour = age), #postnatalOdont
               colour = "gray17", size = 1, linejoin = 'mitre', linetype = 2, arrow = arrow(angle = 30, length = unit(0.03, "npc"), ends = "last", type = "closed"))+
  scale_colour_manual(name = "Growth stage", labels = c("Early Fetus", "Late Fetus", "Postnatal"), #to be ordered as they appear in tibble
                      values = c(mypalette_category[1:2], mypalette_category[4]))+            #legend and color adjustments
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("PC 1 (58.42%)")+ #copy this from standard trajectory plot
  ylab("PC 2 (16.05%)")+
  ggtitle("Trajectories of growth by group")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))  #title font and position

#Visualize plot and save as PDF using menu in bar on the right
group_trajectory_ggplot
