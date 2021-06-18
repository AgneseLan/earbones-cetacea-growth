
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

#CH.2 - DATA IMPORT AND PREP ----

#Import shape data into R
#Takes a while, big file
raw_dat <- read.table("Data/_gpsa_homologized_points_bulla_R.dat")
#Do not check the object, too big. It will display an error if something went wrong

#Import classifiers for analyses, make sure they are in the same order as coordinates
classifiers <- read_csv("Data/bulla_classifiers.csv")
#Check
glimpse(classifiers)

#Add column with specimen code to shape data, needed to assign row names later
raw_dat$code <- classifiers$code

#Import ordination projections (PCOORD from GPSA software) into R - similar to PCA scores
ordination_scores <- read.table("Data/_gpsa_ordination_projections_bulla_R.RData")

#Add column with specimen code to shape data, needed to assign row names later
ordination_scores$code <- classifiers$code

#Import ordination values (PCOORD axis values from GPSA software) into R - amount of variation represented by each axis, for plot
ordination_values <- read.table("Data/_gpsa_ordination_values_bulla_R.RData")

#Make new columns with % or variation for each axis and cumulative variation - useful for plots and if axes need to be excluded later
ordination_values <- ordination_values %>% mutate(total_values = colSums(ordination_values[1]))  %>% 
  mutate(values_100 = (singular_values*100/total_values)) %>% #percentage of variation described by each axis
  mutate(values_cum = cumsum(values_100)) #cumulative variation of that axis plus all the variation described by the previous ones

##Order classifiers by stage, useful for plot legend
#Make factor for variable
classifiers$stage <- factor(classifiers$stage, 
                            levels = c("earlyFetus", "lateFetus", "neonate", "juvenile", "adult")) #copy from string printed with the code above
#Order
classifiers <- classifiers[order(classifiers$stage),]
#Check
glimpse(classifiers)

#Order shape data by stage to match classifiers
classifiers$code

#Make factor for variable
raw_dat$code <- factor(raw_dat$code, #check levels below from string printed with the code above
                       levels = c("Bf1",   "Ff1",    "Gf1" , "Hf3"  , "Hf4"  , "Mf1" ,  "Sf1" ,  "Sf2"  , "Delf1" ,"Delf2", "Delf3", "Lagf1","Lagf2", "Monf1" ,"Phof1" ,"Af3" ,  "Af4" , "Af5" ,  "Af6",  
                                  "Af7",   "Af8",   "Mf2",   "Ddef1",  "Ddef2",  "Lagf3", "Neof1", "Staf1", 
                                  "Gn1",   "Mn1",   "Pn1",   "Glon1", "Stan1", "Phoj1", "Staj1", "Staj2", "Staj3", "Laga1", "Phoa1")) 
#Order
raw_dat <- raw_dat[order(raw_dat$code),]

#Order ordination scores to match classifiers - same list as above
ordination_scores$code <- factor(ordination_scores$code, #check levels below from string printed with the code above
                                 levels = c("Bf1",   "Ff1",    "Gf1" , "Hf3"  , "Hf4"  , "Mf1" ,  "Sf1" ,  "Sf2"  , "Delf1" ,"Delf2", "Delf3", "Lagf1","Lagf2", "Monf1" ,"Phof1" ,"Af3" ,  "Af4" , "Af5" ,  "Af6",  
                                            "Af7",   "Af8",   "Mf2",   "Ddef1",  "Ddef2",  "Lagf3", "Neof1", "Staf1", 
                                            "Gn1",   "Mn1",   "Pn1",   "Glon1", "Stan1", "Phoj1", "Staj1", "Staj2", "Staj3", "Laga1", "Phoa1")) 
#Order
ordination_scores <- ordination_scores[order(ordination_scores$code),]

#Create data frame with only numerical values for shape and set row names as specimen code
raw_dat <- remove_rownames(raw_dat) #remove automatic row names
raw_dat <- data.frame(column_to_rownames(raw_dat, var = "code")) #adding row names to data frame
has_rownames(raw_dat) #check that now it has rownames

#Create data frame with only numerical values for ordination scores and set row names as specimen code
ordination_scores <- remove_rownames(ordination_scores) #remove automatic row names
ordination_scores<- data.frame(column_to_rownames(ordination_scores, var = "code")) #adding row names to data frame
has_rownames(ordination_scores) #check that now it has rownames

#Save specimens code as object
specimens <- row.names(raw_dat)

#Save growth stages as factor, useful for later analysis
stages <- classifiers$stage

#Save groups as factor, useful for later analysis
groups <- classifiers$group

#Transform in 3D array, first number is number of landmarks, second is dimensions (3)
#To find the number of landmark points, divide the number of variables in raw_dat - visible in Environment - by 3
shape_array <- arrayspecs(raw_dat, 16747 , #number in raw_dat divided by 3
                          3) #number of dimensions

#Calculate mean shape coordinates
mean_shape <- mshape(shape_array) 

##Make data frame for analyses in geomorph
gdf <- geomorph.data.frame(coords = shape_array, 
                           code = classifiers$code, group = classifiers$group, 
                           stage = classifiers$stage, category = classifiers$category,
                           bulla_log = classifiers$bullaL_log, periotic_log = classifiers$perioticL_log)
glimpse(gdf)

#Check for specimens that are outliers - if they are very young specimens it is expected
plotOutliers(gdf$coords)
#Save PDF to output folder using menu in bar on the right

##Make palette with ggthemes - color and/or shapes
#Palette from ggthemes_data
mypalette <- ggthemes_data[["tableau"]][["color-palettes"]][["ordered-sequential"]][["Blue"]][["value"]]

#Choose colors for growth stages and image to check
mypalette <- as.matrix(mypalette)
project_palette <- c(mypalette[3,], mypalette[7,], mypalette[11,], mypalette[15,], mypalette[19,])
image(1:5, 1, as.matrix(1:5), col = project_palette, xlab = "Project palette", ylab = "", yaxt = "n")

#Create shape palette for groups
shapes <- c(15,19) #these are a square and a circle, use ?pch to see more shapes

#Save the work
#Press the "Save all documents open" double floppy disk above and run this code
save.image("D:/Agnese/ear bones project/R/bulla_R/.RData")
savehistory("D:/Agnese/ear bones project/R/bulla_R/.Rhistory")

#CH. 3 - PCA COMPLETE DATASET ----

#Run PCA on complete dataset
PCA_all <- gm.prcomp(gdf$coords)

#List of PC components and proportion of variation
PCA_all 

##View basic R plot to check - use to copy values for PC1 of PC2 axes
PCA_all_plot <- plot(PCA_all, main = "PCA all data",  pch = 21, #title and type of point to be used
                     col = "deeppink",   #border of points
                     bg = "deeppink",    #fill of points
                     cex = 1,            #size of points (1=regular)
                     font.main = 2)       #bold font for title
#Add quick labels to plot
text(x = PCA_all$x[,1], y = PCA_all$x[,2], labels = rownames(PCA_all$x), 
     pos = 1,       #position relative to data point
     offset = 0.5,  #distance from data point
     cex = 0.75)    #font size (1=regular)

#Save PC scores as object to use later
pcscores_all <- PCA_all$x 

##Make better plot using ggplot
#Read PC scores as tibble
pcscores_all <- as_tibble(pcscores_all)
#Add labels and other attributes to tibble as columns
pcscores_all <- pcscores_all %>% mutate(specimens = gdf$code, stage = gdf$stage, group = gdf$group, category = gdf$category)
glimpse(pcscores_all)

#Nice PCA plot with stages and groups
PCA_all_ggplot <- ggplot(pcscores_all, aes(x = Comp1, y = Comp2, label = specimens, colour = stage, shape = group))+
  geom_point(size = 3)+
  geom_text_repel(colour = "black", size = 3.5, max.overlaps = 40)+
  scale_colour_manual(name = "Growth stage", labels =  c("Early Fetus", "Late Fetus", "Neonate", "Juvenile", "Adult"), #to be ordered as they appear in tibble
                      values = project_palette)+            #legend and color adjustments
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("PC 1 (20.16%)")+ #copy this from standard PCA plot (PCA_all_plot)
  ylab("PC 2 (12.8%)")+
  ggtitle("PCA all data")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_ggplot

#Make hulls for PCA plot with hulls around growth stages
hulls_all <- pcscores_all %>%
  group_by(stage) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)
glimpse(hulls_all)

#Nice PCA plot with hulls around stages
PCA_all_hulls_ggplot <- ggplot(pcscores_all, aes(x = Comp1, y = Comp2, label = specimens, colour = stage, fill = stage))+
  geom_point(size = 3, aes(shape = group))+
  scale_colour_manual(name = "Growth stage", labels =  c("Early Fetus", "Late Fetus", "Neonate", "Juvenile", "Adult"), #to be ordered as they appear in tibble
                      values = project_palette)+            #legend and color adjustments
  geom_polygon(data = hulls_all, aes(x = x, y = y, fill = stage), alpha = .5, show.legend = FALSE)+ #colored hulls with transparency
  scale_fill_manual(name = "Growth stage", labels = c("Early Fetus", "Late Fetus", "Neonate", "Juvenile", "Adult"),
                    values =  project_palette)+ #must match scale_colour_manual
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("PC 1 (20.16%)")+ #copy this from standard PCA plot (PCA_all_plot)
  ylab("PC 2 (12.8%)")+
  ggtitle("PCA all data")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
  geom_text_repel(colour = "black", size = 3.5, max.overlaps = 40)

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_hulls_ggplot  

#Make hulls for PCA plot with hulls around categories (3 stages: early fetus, late fetus, postnatal)
hulls_all_category <- pcscores_all %>%
  group_by(category) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)
glimpse(hulls_all_category)

#Nice PCA plot with hulls around categories (3 stages: early fetus, late fetus, postnatal)
PCA_all_category_ggplot <- ggplot(pcscores_all, aes(x = Comp1, y = Comp2, label = specimens, colour = category, fill = category))+
  geom_point(size = 3, aes(shape = group))+
  scale_colour_manual(name = "Growth stage", labels =  c("Early Fetus", "Late Fetus", "Postnatal"), #to be ordered as they appear in tibble
                      values = c(project_palette[1:2], project_palette[4]))+            #legend and color adjustments
  geom_polygon(data = hulls_all_category, aes(x = x, y = y, fill = category), alpha = .5, show.legend = FALSE)+ #colored hulls with transparency
  scale_fill_manual(name = "Growth stage", labels = c("Early Fetus", "Late Fetus", "Postnatal"),
                    values =  c(project_palette[1:2], project_palette[4]))+ #must match scale_colour_manual
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("PC 1 (20.16%)")+ #copy this from standard PCA plot (PCA_all_plot)
  ylab("PC 2 (12.8%)")+
  ggtitle("PCA all data")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  geom_text_repel(colour = "black", size = 3.5, max.overlaps = 40)

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_category_ggplot

##PCA plot for each stage
#Create one tibble for each stage
#First split between stages
pcscores_all_stages <- pcscores_all %>% group_by(category) %>% group_split()

#Save as separate tibbles - 1 for early fetus, 1 for late fetus, 1 for all 3 postnatal stages (neonate, juvenile, adult)
pcscores_all_earlyfetus <- pcscores_all_stages[[2]]
pcscores_all_latefetus <- pcscores_all_stages[[3]]
pcscores_all_postnatal <- pcscores_all_stages[[1]]

#Nice PCA plot with groups early fetus
PCA_all_earlyfetus_ggplot <- ggplot(pcscores_all_earlyfetus, aes(x = Comp1, y = Comp2, label = specimens, shape = group))+
  geom_point(size = 3, color = project_palette[1])+
  geom_text_repel(colour = "black", size = 3.5, max.overlaps = 40)+
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("PC 1 (20.16%)")+ #copy this from standard PCA plot (PCA_all_plot)
  ylab("PC 2 (12.8%)")+
  ggtitle("PCA early fetus")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_earlyfetus_ggplot

#Nice PCA plot with groups late fetus
PCA_all_latefetus_ggplot <- ggplot(pcscores_all_latefetus, aes(x = Comp1, y = Comp2, label = specimens, shape = group))+
  geom_point(size = 3, color = project_palette[2])+
  geom_text_repel(colour = "black", size = 3.5, max.overlaps = 40)+
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("PC 1 (20.16%)")+ #copy this from standard PCA plot (PCA_all_plot)
  ylab("PC 2 (12.8%)")+
  ggtitle("PCA late fetus")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_latefetus_ggplot

#Nice PCA plot with groups postnatal
PCA_all_postnatal_ggplot <- ggplot(pcscores_all_postnatal, aes(x = Comp1, y = Comp2, label = specimens, colour = stage, shape = group))+
  geom_point(size = 3)+
  geom_text_repel(colour = "black", size = 3.5, max.overlaps = 40)+
  scale_colour_manual(name = "Growth stage", labels =  c("Neonate", "Juvenile", "Adult"), #to be ordered as they appear in tibble
                      values = project_palette[3:5])+            #legend and color adjustments
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("PC 1 (20.16%)")+ #copy this from standard PCA plot (PCA_all_plot)
  ylab("PC 2 (12.8%)")+
  ggtitle("PCA postnatal")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_postnatal_ggplot


##Regression PC1 and PC2 vs bulla length - check if size an important factor in results overall
#Create data frame with data
pcscores_all_size <- pcscores_all %>% mutate(size = gdf$bulla_log)

#Calculate regression for each component
reg_PC1all_size <- lm(Comp1 ~ size, data = pcscores_all_size)
reg_PC2all_size <- lm(Comp2 ~ size, data = pcscores_all_size)

#View results and p-value
summary(reg_PC1all_size)
summary(reg_PC2all_size)

#Save results of significant regression to file
sink("Output/PC1all_size_lm.txt")
print(summary(reg_PC1all_size))
sink() 
#Repeat for each significance regression
sink("Output/PC2all_size_lm.txt")
print(summary(reg_PC2all_size))
sink() 


#Save the work
#Press the "Save all documents open" double floppy disk above and run this code
save.image("D:/Agnese/ear bones project/R/bulla_R/.RData")
savehistory("D:/Agnese/ear bones project/R/bulla_R/.Rhistory")

#CH. 4 - ALLOMETRY CORRECTION ----
##Evaluate allometry and get the allometry-free shapes using ln bulla length
#Important if PCs highly correlated with size

#Regression shape on size
allometry <- procD.lm(coords ~ bulla_log, data = gdf, iter=999, print.progress = TRUE) 

#Check results and p-value for significance
summary(allometry)

#Save results of significant regression to file
sink("Output/allometry_shape_size.txt")
print(summary(allometry))
sink() 

#Create residuals array to then save as coordinates for analyses
allometry_array <- arrayspecs(allometry$residuals,p = dim(gdf$coords)[1], k = dim(gdf$coords)[2]) 

#New shapes adjusted for allometry with CS to use in analyses
allometry_residuals <- allometry_array + array(mean_shape, dim(allometry_array)) 

#Save mean shape of allometry-adjusted shapes to use in a analyses
mean_shape_residuals <- mshape(allometry_residuals)

#Regression score of shape vs bulla length or periotic length  - regression method with "RegScore" plotting
allometry_plot_regscore <- plot(allometry, type = "regression",predictor = gdf$bulla_log, reg.type = "RegScore",
                                main = "Shape vs Ln(bulla length)",xlab = "Ln(bulla length)", pch = 21, col = "chartreuse4", bg = "chartreuse4", cex = 1.2, font.main = 2)   #improve graphics

##Add regression line with confidence intervals to plot
#Create object to use for linear model
allometry_regscores <- allometry_plot_regscore[["RegScore"]] 

##Make better allometry plot with ggplot
#Create data frame object that ggplot can read - use data from plot object you want to improve
allometry_plot <- data.frame(size = allometry_plot_regscore[["plot.args"]][["x"]], RegScores = allometry_plot_regscore[["plot.args"]][["y"]])
allometry_plot

#Convert data frame to tibble
allometry_plot <- as_tibble(allometry_plot)
#Add labels and other attributes to tibble as columns
allometry_plot <- allometry_plot %>% mutate(specimens = gdf$code, stage = gdf$stage, group = gdf$group, category = gdf$category)
glimpse(allometry_plot)

#Nice plot with specimens colored by age AND regression line with confidence intervals
allometry_ggplot <- ggplot(allometry_plot, aes(x = size, y = RegScores, label = specimens, colour = stage, shape = group))+
  geom_smooth(aes(x = size, y = RegScores), method = 'lm', inherit.aes = F,         #confidence intervals and reg line, before points
              colour = "darkblue", fill = 'gainsboro', linetype = "dashed", size = 0.5)+ #should be straight regression line with confidence interval in grey
  geom_point(size = 3)+       #points after, so they are on top
  scale_colour_manual(name = "Growth stage", labels =  c("Early Fetus", "Late Fetus", "Neonate", "Juvenile", "Adult"), #to be ordered as they appear in tibble
                      values = project_palette)+            #legend and color adjustments          
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_classic(base_size = 12)+
  xlab("Ln(Bulla length)")+
  ylab("Regression Score")+
  ggtitle ("Allometry plot - p-value = 0.001**")+  #copy from model summary
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  geom_text_repel(colour = "black", size = 3.5,          #label last so that they are on top of fill
                  force_pull = 3, point.padding = 1)     #position of tables relative to point (proximity and distance) 

#Visualize plot and save as PDF using menu in bar on the right
allometry_ggplot

#Save the work
#Press the "Save all documents open" double floppy disk above and run this code
save.image("D:/Agnese/ear bones project/R/bulla_R/.RData")
savehistory("D:/Agnese/ear bones project/R/bulla_R/.Rhistory")

#CH. 5 - ORDINATION PLOT ----

#Plot the ordination scores imported from the GPSA
#Transform ordination scores as tibble
ordination_scores <- as_tibble(ordination_scores)
#Add labels and other attributes to tibble as columns
ordination_scores <- ordination_scores %>% mutate(specimens = gdf$code, stage = gdf$stage, group = gdf$group)
glimpse(ordination_scores)

#Recall values to copy for plot axis labels - use values_100 column for the first 2 axes
ordination_values

#Nice plot - axis 1 and 2
ordination_scores_ggplot <- ggplot(ordination_scores, aes(x = axis1, y = axis2, label = specimens, colour = stage, shape = group))+
  geom_point(size = 3)+
  geom_text_repel(colour = "black", size = 3.5, max.overlaps = 40)+
  scale_colour_manual(name = "Growth stage", labels =  c("Early Fetus", "Late Fetus", "Neonate", "Juvenile", "Adult"), #to be ordered as they appear in tibble
                      values = project_palette)+            #legend and color adjustments
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("I (29.3%)")+ #copy this from ordination values_100 column printed before
  ylab("II (17.8%)")+
  ggtitle("PCOORD GPSA")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12))  #title font and position

#Visualize plot and save as PDF using menu in bar on the right
ordination_scores_ggplot

#CH. 6 - ANOVA SHAPE GROUPS  ----

##Conduct ANOVA to test if there is significant difference in allometry between groups
allometry_group_anova <- procD.lm(coords ~ bulla_log * group, iter = 999, data = gdf, RRPP = F)

#Check results
summary(allometry_group_anova)

#Save results to file
sink("Output/allometry_group_anova.txt")
summary(allometry_group_anova)
sink() 

#Plot of allometry with group differences - use previous reg scores
allometry_group_ggplot <- ggplot(allometry_plot, aes(x = size, y = RegScores, label = specimens, colour = stage, shape = group))+
  geom_smooth(aes(x = size, y = RegScores, linetype = group), method = 'lm', inherit.aes = F,         #confidence intervals and reg line, before points
              colour = "darkblue", fill = 'gainsboro', size = 0.5, show.legend = F)+ #should be straight regression line with confidence interval in grey
  geom_point(size = 3)+       #points after, so they are on top
  scale_colour_manual(name = "Growth stage", labels =  c("Early Fetus", "Late Fetus", "Neonate", "Juvenile", "Adult"), #to be ordered as they appear in tibble
                      values = project_palette)+            #legend and color adjustments          
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_classic(base_size = 12)+
  xlab("Ln(Bulla length)")+
  ylab("Regression Score")+
  ggtitle ("Allometry plot - p-value = 0.001**")+  #copy from model summary
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  geom_text_repel(colour = "black", size = 3.5,          #label last so that they are on top of fill
                  force_pull = 3, point.padding = 1)     #position of tables relative to point (proximity and distance) 

#Visualize plot and save as PDF using menu in bar on the right
allometry_group_ggplot

##Test if there is a significant difference between allometry general model and allometry considering groups
#Create model with group as additive with no interaction
allometry_add_group_anova <- procD.lm(coords ~ bulla_log + group, iter = 999, data = gdf, RRPP = F)

#Check results
summary(allometry_add_group_anova)

#ANOVA - is a model significantly better than the others?

allometry_models_anova <- anova(allometry, allometry_add_group_anova, allometry_group_anova)
allometry_models_anova

#Save results to file
sink("Output/allometry_models_anova.txt")
allometry_models_anova
sink() 

##Conduct ANOVA to test if there is significant shape variation among groups for each category
#Are shapes in each group different from the other groups at different growth stages?
#2 analyses, use shape data and pc scores
#Shapes
#Make 3 data frames, 1 for each category - the numbers are the rows corresponding to specimens at that growth stage
gdf_earlyfetus <- geomorph.data.frame(coords = shape_array[,,1:15], 
                                      code = classifiers$code[1:15], group = classifiers$group[1:15], 
                                      stage = classifiers$stage[1:15], category = classifiers$category[1:15],
                                      bulla_log = classifiers$bullaL_log[1:15], periotic_log = classifiers$perioticL_log[1:15])
glimpse(gdf_earlyfetus)

gdf_latefetus <- geomorph.data.frame(coords = shape_array[,,16:27], 
                                     code = classifiers$code[16:27], group = classifiers$group[16:27], 
                                     stage = classifiers$stage[16:27], category = classifiers$category[16:27],
                                     bulla_log = classifiers$bullaL_log[16:27], periotic_log = classifiers$perioticL_log[16:27])
glimpse(gdf_latefetus)

gdf_postnatal <- geomorph.data.frame(coords = shape_array[,,28:38], 
                                     code = classifiers$code[28:38], group = classifiers$group[28:38], 
                                     stage = classifiers$stage[28:38], category = classifiers$category[28:38],
                                     bulla_log = classifiers$bullaL_log[28:38], periotic_log = classifiers$perioticL_log[28:38])
glimpse(gdf_postnatal)

#Make separate dataframe to check shape is not different postnatal in odontocetes
gdf_postnatal_odont <- geomorph.data.frame(coords = shape_array[,,31:38], 
                                           code = classifiers$code[31:38], group = classifiers$group[31:38], 
                                           stage = classifiers$stage[31:38], category = classifiers$category[31:38],
                                           bulla_log = classifiers$bullaL_log[31:38], periotic_log = classifiers$perioticL_log[31:38])
glimpse(gdf_postnatal_odont)

#Start by making sure it is ok to lump postnatal stages with ANOVA of shape per stage in postnatal odontocetes
shape_stage_odont_anova <- procD.lm(coords ~ stage, iter = 999, data = gdf_postnatal_odont, RRPP = F)

#Check results - should not be significant
summary(shape_stage_odont_anova)

#Save results to file
sink("Output/shape_stage_odont_anova.txt")
summary(shape_stage_odont_anova)
anova(shape_stage_odont_anova)
sink() 

#Proceed with shape differences between groups for each stage
#Early fetus
shape_group_earlyfetus_anova <- procD.lm(coords ~ group, iter = 999, data = gdf_earlyfetus, RRPP = F)

#Check results
summary(shape_group_earlyfetus_anova)

#Save results to file
sink("Output/shape_group_earlyfetus_anova.txt")
summary(shape_group_earlyfetus_anova)
anova(shape_group_earlyfetus_anova)
sink() 

#Late fetus
shape_group_latefetus_anova <- procD.lm(coords ~ group, iter = 999, data = gdf_latefetus, RRPP = F)

#Check results
summary(shape_group_latefetus_anova)

#Save results to file
sink("Output/shape_group_latefetus_anova.txt")
summary(shape_group_latefetus_anova)
anova(shape_group_latefetus_anova)
sink() 

#Postnatal
shape_group_postnatal_anova <- procD.lm(coords ~ group, iter = 999, data = gdf_postnatal, RRPP = F)

#Check results
summary(shape_group_postnatal_anova)

#Save results to file
sink("Output/shape_group_postnatal_anova.txt")
summary(shape_group_postnatal_anova)
anova(shape_group_postnatal_anova)
sink() 

#PC scores
#Repeat same analyses but using pc scores for each stage
#Early fetus
pcscores_group_earlyfetus_anova <- procD.lm(pcscores_all_earlyfetus[1:38] #this number is the number of Components from the PCA, make sure it is correct by calling PCA_all$x
                                            ~ group, iter = 999, data = pcscores_all_earlyfetus, RRPP = F)

#Check results
summary(pcscores_group_earlyfetus_anova)

#Save results to file
sink("Output/pcscores_group_earlyfetus_anova.txt")
summary(pcscores_group_earlyfetus_anova)
anova(pcscores_group_earlyfetus_anova)
sink() 

#Late fetus
pcscores_group_latefetus_anova <- procD.lm(pcscores_all_latefetus[1:38] #this number is the number of Components from the PCA, make sure it is correct by calling PCA_all$x
                                           ~ group, iter = 999, data = pcscores_all_latefetus, RRPP = F)

#Check results
summary(pcscores_group_latefetus_anova)

#Save results to file
sink("Output/pcscores_group_latefetus_anova.txt")
summary(pcscores_group_latefetus_anova)
anova(pcscores_group_latefetus_anova)
sink() 

#Postnatal
pcscores_group_postnatal_anova <- procD.lm(pcscores_all_postnatal[1:29] #this number is the number of Components from the PCA, make sure it is correct by calling PCA_all$x
                                           ~ group, iter = 999, data = pcscores_all_postnatal, RRPP = F)

#Check results
summary(pcscores_group_postnatal_anova)

#Save results to file
sink("Output/pcscores_group_postnatal_anova.txt")
summary(pcscores_group_postnatal_anova)
anova(pcscores_group_postnatal_anova)
sink() 

