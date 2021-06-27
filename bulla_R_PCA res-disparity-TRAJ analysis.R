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
library(RColorBrewer)
library(ggthemes)
library(ggpubr)
library(ggplotify)
library(Morpho)
library(rphylopic)
library(png)
library(gridExtra)

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
pcscores_res <- pcscores_res %>% mutate(specimens = gdf$code, taxon = gdf$taxon, group = gdf$group, category = gdf$category)
glimpse(pcscores_res)

#Nice plot with labels and specimens colored by category
PCA_res_ggplot <- ggplot(pcscores_res, aes(x = Comp1, y = Comp2, label = specimens, colour = category, shape = group))+
  geom_point(size = 3)+
  geom_text_repel(colour = "black", size = 3.5, max.overlaps = 60)+
  scale_colour_manual(name = "Growth category", labels =  c("Early Fetus", "Late Fetus", "Postnatal"), #to be ordered as they appear in tibble
                      values = mypalette_category)+            #legend and color adjustments
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("PC 1 (11.46%)")+ #copy this from standard PCA plot (PCA_res_plot)
  ylab("PC 2 (7.45%)")+
  ggtitle("PCA residuals")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_res_ggplot

#Make hulls for PCA plot with hulls around growth category
hulls_res <- pcscores_res %>%
  group_by(category) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)
glimpse(hulls_res)

#Nice PCA plot with hulls around categories
PCA_res_hulls_ggplot <- ggplot(pcscores_res, aes(x = Comp1, y = Comp2, label = specimens, colour = category, fill = category))+
  geom_point(size = 3, aes(shape = group))+
  scale_colour_manual(name = "Growth category", labels =  c("Early Fetus", "Late Fetus", "Postnatal"), #to be ordered as they appear in tibble
                      values = mypalette_category)+            #legend and color adjustments
  geom_polygon(data = hulls_res, aes(x = x, y = y, fill = category), alpha = .5, show.legend = FALSE)+ #colored hulls with transparency
  scale_fill_manual(name = "Growth category", labels = c("Early Fetus", "Late Fetus", "Postnatal"),
                    values =  mypalette_category)+ #must match scale_colour_manual
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("PC 1 (11.46%)")+ #copy this from standard PCA plot (PCA_res_plot)
  ylab("PC 2 (7.45%)")+
  ggtitle("PCA residuals")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

#Visualize plot and save as PDF using menu in bar on the right
PCA_res_hulls_ggplot  

##Plots for each category ----
#PCA plots for each category
#Create one tibble for each category
#First split between categories
pcscores_res_category <- pcscores_res %>% group_by(category) %>% group_split()
#Check order
View(pcscores_res_category)

#Save as separate tibbles - 1 for early fetus, 1 for late fetus, 1 for postnatal
pcscores_res_earlyfetus <- pcscores_res_category[[1]]
pcscores_res_latefetus <- pcscores_res_category[[2]]
pcscores_res_postnatal <- pcscores_res_category[[3]]

#Nice PCA plot with groups early fetus
PCA_res_earlyfetus_ggplot <- ggplot(pcscores_res_earlyfetus, aes(x = Comp1, y = Comp2, label = specimens, shape = group))+
  geom_point(size = 3, color = mypalette_category[1])+
  geom_text_repel(colour = "black", size = 3.5, max.overlaps = 40)+
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("PC 1 (11.46%)")+ #copy this from standard PCA plot (PCA_res_plot)
  ylab("PC 2 (7.45%)")+
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
  xlab("PC 1 (11.46%)")+ #copy this from standard PCA plot (PCA_res_plot)
  ylab("PC 2 (7.45%)")+
  ggtitle("PCA late fetus")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_res_latefetus_ggplot

#Nice PCA plot with groups postnatal
PCA_res_postnatal_ggplot <- ggplot(pcscores_res_postnatal, aes(x = Comp1, y = Comp2, label = specimens, shape = group))+
  geom_point(size = 3, color = mypalette_category[3])+
  geom_text_repel(colour = "black", size = 3.5, max.overlaps = 40)+
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("PC 1 (11.46%)")+ #copy this from standard PCA plot (PCA_res_plot)
  ylab("PC 2 (7.45%)")+
  ggtitle("PCA postnatal")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_res_postnatal_ggplot

##Plots for well sampled taxa only ----
##PCA plot including only the 4 well sampled taxa

#Make 1 tibble to include only the selected taxa
pcscores_taxa_res <- pcscores_res %>% filter(taxon %in% best_taxa)
#Replace to have only 1 taxon name for minkes
pcscores_taxa_res[pcscores_taxa_res == "B.acutorostrata"] <- "B.bonaerensis"

#Nice PCA plot with only  well sampled taxa
PCA_taxa_res_ggplot <- ggplot(pcscores_taxa_res, aes(x = Comp1, y = Comp2, shape = group, colour = taxon, alpha = category, label = specimens))+
  geom_point(size = 3)+
  geom_text_repel(colour = "black", size = 3.5, max.overlaps = 60, show.legend=FALSE)+
  scale_colour_manual(name = "Taxa", labels =  c("B.bonaerensis", "B.physalus", "Ph.phocoena", "St.attenuata"), #to be ordered as they appear in tibble
                      values = mypalette_taxa)+            #legend and color adjustments
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  scale_alpha_manual(name = "Growth category", labels =  c("Early Fetus", "Late Fetus", "Postnatal"), values = c(0.4, 0.7, 1))+
  theme_bw()+
  xlab("PC 1 (11.46%)")+ #copy this from standard PCA plot (PCA_res_plot)
  ylab("PC 2 (7.45%)")+
  ggtitle("PCA residuals selected taxa")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_taxa_res_ggplot

#Make hulls for PCA plot with hulls around taxa
hulls_taxa_res <- pcscores_taxa_res %>%
  group_by(taxon) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)
glimpse(hulls_taxa_res)

#Nice PCA plot with hulls around categories
PCA_taxa_res_hulls_ggplot <- ggplot(pcscores_taxa_res, aes(x = Comp1, y = Comp2, colour = taxon, alpha = category))+
  geom_point(size = 3, aes(shape = group))+
  scale_colour_manual(name = "Taxa", labels =  c("B.bonaerensis", "B.physalus", "Ph.phocoena", "St.attenuata"), #to be ordered as they appear in tibble
                      values = mypalette_taxa)+            #legend and color adjustments
  scale_alpha_manual(name = "Growth category", labels =  c("Early Fetus", "Late Fetus", "Postnatal"), values = c(0.4, 0.7, 1))+
  geom_polygon(data = hulls_taxa_res, aes(x = x, y = y, fill = taxon), alpha = .2, show.legend = FALSE)+ #colored hulls with transparency
  scale_fill_manual(name = "Taxa", labels =  c("B.bonaerensis", "B.physalus", "Ph.phocoena", "St.attenuata"), #to be ordered as they appear in tibble
                    values = mypalette_taxa)+ #must match scale_colour_manual
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("PC 1 (11.46%)")+ #copy this from standard PCA plot (PCA_res_plot)
  ylab("PC 2 (7.45%)")+
  ggtitle("PCA residuals selected taxa")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

#Visualize plot and save as PDF using menu in bar on the right
PCA_taxa_res_hulls_ggplot  

##Regression PC1 and PC2 res vs group and category ----
#Check if main components may be associated with group (Mysticeti/Odontoceti)
reg_PC1res_group <- lm(Comp1 ~ group, data = pcscores_res)
reg_PC2res_group <- lm(Comp2 ~ group, data = pcscores_res)

#View results and p-value
summary(reg_PC1res_group)
summary(reg_PC2res_group)

#Save results of significant regression to file
sink("bulla_R/PC1-PC2_res_group_lm.txt")
print("PC1")
summary(reg_PC1res_group)

print("PC2")
summary(reg_PC2res_group)
sink() 

#Check if main components may be associated with category (early, late, postnatal)
reg_PC1res_category <- lm(Comp1 ~ category, data = pcscores_res)
reg_PC2res_category <- lm(Comp2 ~ category, data = pcscores_res)

#View results and p-value
summary(reg_PC1res_category)
summary(reg_PC2res_category)

#Save results of significant regression to file
sink("bulla_R/PC1-PC2_res_category_lm.txt")
print("PC1")
summary(reg_PC1res_category)

print("PC2")
summary(reg_PC2res_category)
sink() 

#CH. 8 - MORPHOLOGICAL DISPARITY (pc scores) FOR GROUPS, TAXA, GROWTH CATEGORY - only well sampled taxa used for analysis -----
#Calculate Procrustes variances and distances between groups, with p-value for each pair of groups
#How different are each group shapes compared to other group shapes?

#Disparity between categories, not considering groups or taxa
category_disparity <- morphol.disparity(pcscores_taxa[1:75] ~ 1, groups = ~ categories_taxa, iter = 999)

#Results and significance
summary(category_disparity)

#Save results to file
sink("bulla_R/category_disparity.txt")
print(summary(category_disparity))
sink() 

#Disparity between categories, considering groups
category_group_disparity <- morphol.disparity(pcscores_taxa[1:75] ~ 1, groups = ~ categories_taxa*groups_taxa, iter = 999)

#Results and significance
summary(category_group_disparity)

#Save results to file
sink("bulla_R/category_group_disparity.txt")
print(summary(category_group_disparity))
sink() 

#Disparity between categories, considering taxa
category_taxa_disparity <- morphol.disparity(pcscores_taxa[1:75] ~ 1, groups = ~ categories_taxa*taxa_taxa, iter = 999)

#Results and significance
summary(category_taxa_disparity)

#Save results to file
sink("bulla_R/category_taxa_disparity.txt")
print(summary(category_taxa_disparity))
sink() 

#Disparity between groups for each stage - use objects with pc scores for each category (early fetus, late fetus) from ANOVAs
#Create objects for groups for each category (early and late) from pc scores dataset
groups_taxa_early <- pcscores_taxa_size_early$group

groups_taxa_late <- pcscores_taxa_size_late$group

#Early fetus
group_disparity_early <- morphol.disparity(pcscores_taxa_size_early[1:75] ~ 1, groups = ~ groups_taxa_early, iter = 999)

#Results and significance
summary(group_disparity_early)

#Save results to file
sink("bulla_R/group_disparity_early.txt")
summary(group_disparity_early)
sink() 

#Late fetus
group_disparity_late <- morphol.disparity(pcscores_taxa_size_late[1:75] ~ 1, groups = ~ groups_taxa_late, iter = 999)

#Results and significance
summary(group_disparity_late)

#Save results to file
sink("bulla_R/group_disparity_latefetus.txt")
summary(group_disparity_late)
sink() 

#CH. 9 - TRAJECTORY ANALYSIS - entire dataset used for analysis ----
#Shows trajectories of shape variation, possible to compare trajectories of growth between odontocetes and mysticetes
#Use original shape data, entire dataset to account for uneven sampling
#Uneven sampling doesn't allow to look at the best sampled taxa separately

#First perform procD.lm to create linear model that describes what we are trying to test - shape changes at each stage (category) considering the 2 groups
fit_category_group <- procD.lm(coords ~ group * category, iter = 499, data = gdf, RRPP = F)

#Check that there is a significant correlation
summary(fit_category_group)

#Use fit to calculate trajectories
group_trajectory <- trajectory.analysis(fit_category_group, groups = gdf$group, traj.pts = gdf$category, 
                                        pca = TRUE, print.progress = TRUE) 

#View results
#Magnitude differences between trajectories, standard summary -  difference in path lengths of trajectories between the 2 groups
summary(group_trajectory, show.trajectories = TRUE, attribute = "MD") 
#Trajectory correlations - angular differences between trajectory principal axes for the 2 groups
summary(group_trajectory, show.trajectories = TRUE, attribute = "TC", angle.type = "deg")
#Trajectory shape differences - "Procrustes" distance between trajectories of the 2 groups
summary(group_trajectory, show.trajectories = TRUE, attribute = "SD") 

#Save results to file
sink("bulla_R/group_trajectory.txt")
print("Magnitude difference (difference in path lengths of trajectories)")
summary(group_trajectory, show.trajectories = TRUE, attribute = "MD") 
print("Correlations (angles) between trajectories (angular differences between trajectory principal axes)")
summary(group_trajectory, show.trajectories = TRUE, attribute = "TC", angle.type = "deg")
print("Shape differences between trajectory vectors (Procrustes distance between trajectories)")
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
legend(x= -15, y = -10, legend = c("Mysticeti","Odontoceti"), pch =  c(21, 22), pt.bg = 1)

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
group_trajectory_ggplot <- ggplot(group_trajectory_pcscores, aes(x = PC1, y = PC2, shape = group, alpha = category))+
  geom_point(size = 3, colour = "darkgray")+
  geom_point(data = group_trajectory_pcscores_means, aes(x = x, y = y, colour = group, shape = group, alpha = category), size = 5, inherit.aes = F)+
  #add trajectory lines, one line for each, write row number from tibble, should be in order as legend of plot
  geom_segment(data = group_trajectory_pcscores_means, aes(x = x[1], y = y[1],  #earlyfetusMyst
                                                           xend =  x[2], yend = y[2], colour = age), #latefetusMyst
               colour = mypalette_taxa[2], size = 1.2, linejoin = 'mitre', linetype = 1, 
               arrow = arrow(angle = 30, length = unit(0.03, "npc"), ends = "last", type = "closed"), show.legend = F)+
  #add arrow at end  
  geom_segment(data = group_trajectory_pcscores_means, aes(x = x[2], y = y[2], #latefetusMyst
                                                           xend = x[3], yend = y[3], colour = age), #postnatalMyst
               colour = mypalette_taxa[2], size = 1.2, linejoin = 'mitre', linetype = 1, 
               arrow = arrow(angle = 30, length = unit(0.03, "npc"), ends = "last", type = "closed"), show.legend = F)+
  #add trajectory lines, one line for each, write row number from tibble, should be in order as legend of plot
  geom_segment(data = group_trajectory_pcscores_means, aes(x = x[4], y = y[4],  #earlyfetusOdont
                                                           xend =  x[5], yend = y[5], colour = age), #latefetusOdont
               colour = mypalette_taxa[4], size = 1.2, linejoin = 'mitre', linetype = 2, 
               arrow = arrow(angle = 30, length = unit(0.03, "npc"), ends = "last", type = "closed"), show.legend = F)+
  #add arrow at end  
  geom_segment(data = group_trajectory_pcscores_means, aes(x = x[5], y = y[5], #latefetusOdont
                                                           xend = x[6], yend = y[6], colour = age), #postnatalOdont
               colour = mypalette_taxa[4], size = 1.2, linejoin = 'mitre', linetype = 2, 
               arrow = arrow(angle = 30, length = unit(0.03, "npc"), ends = "last", type = "closed"), show.legend = F)+
  scale_alpha_manual(name = "Growth stage", labels = c("Early Fetus", "Late Fetus", "Postnatal"), #to be ordered as they appear in tibble
                     values = c(0.4, 0.7, 1))+            #legend and color adjustments
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  scale_colour_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = c(mypalette_taxa[2], mypalette_taxa[4]))+
  theme_bw()+
  xlab("PC 1 (60.24%)")+ #copy this from standard trajectory plot
  ylab("PC 2 (14.8%)")+
  ggtitle("Tympanic bulla")+
  guides(shape = guide_legend(label = F, title = NULL, override.aes = list(shape = NA)), 
         colour = guide_legend(label = F, title = NULL))+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16), legend.position = "bottom", legend.direction = "horizontal")
#Add silhouettes groups
group_trajectory_ggplot <- group_trajectory_ggplot   + 
  add_phylopic(B.physalus, alpha = 1, x = -5, y = -7, ysize = 1.8, color = mypalette_taxa[2])+
  add_phylopic(St.attenuata, alpha = 1, x = 0, y = 2.5, ysize = 1.7, color = mypalette_taxa[4])
#Visualize plot and save as PDF using menu in bar on the right
group_trajectory_ggplot

