## 3 - SHAPE DATA ANALYSIS (GPSA) ##
# Periotic dataset

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

#CH.2 - DATA IMPORT AND PREP ----

#Import shape data into R
#Takes a while, big file
raw_dat <- read.table("Data/_gpsa_homologized_points_periotic_R.dat")
#Do not check the object, too big. It will display an error if something went wrong

#Import classifiers for analyses, make sure they are in the same order as coordinates
classifiers <- read_csv("Data/periotic_classifiers.csv")
#Check
glimpse(classifiers)

#Add column with specimen code to shape data, needed to assign row names later
raw_dat$code <- classifiers$code

#Import ordination projections (PCOORD from GPSA software) into R - similar to PCA scores
ordination_scores <- read.table("Data/_gpsa_ordination_projections_periotic_R.RData")

#Add column with specimen code to shape data, needed to assign row names later
ordination_scores$code <- classifiers$code

#Import ordination values (PCOORD axis values from GPSA software) into R - amount of variation represented by each axis, for plot
ordination_values <- read.table("Data/_gpsa_ordination_values_periotic_R.RData", header = T)
#Check
glimpse(ordination_values)

#Make new columns with % of variation for each axis and cumulative variation - useful for plots and if axes need to be excluded later
ordination_values <- ordination_values %>% 
  mutate(values_100 = (pro*100)) %>% #percentage of variation described by each axis
  mutate(values_cum = (cum*100)) #cumulative variation of that axis plus all the variation described by the previous ones
#Check
glimpse(ordination_values)

#Sink values to file in more readable format
sink("periotic_R/ordination_values.txt", append = F)
print(ordination_values)
sink()

##Order classifiers by category, useful for plot legend
#Make factor for variable
classifiers$category <- factor(classifiers$category, 
                               levels = c("early", "late", "postnatal")) 
#Order
classifiers <- classifiers[order(classifiers$category),]
#Check
glimpse(classifiers)

#Order shape data by category to match classifiers
classifiers$code

#Make factor for variable
raw_dat$code <- factor(raw_dat$code, #check levels below from string printed with the code above
                       levels = c("Bf1" ,  "Ff6"  , "Ff7" ,  "Hf8"  , "Hf9"  , "Sf2" ,  "Wf1",   "Kogf2" ,
                                  "Lagf1" ,"Lagf2", "Monf1", "Phof1", "Phof5", "Phof6", "Staf4" ,"Staf5", "Af3",
                                  "Af4" ,  "Af5"  , "Af6" ,  "Af7"  , "Af8"  , "Mf2" ,  "Ddef1", "Ddef2",
                                  "Glof3" ,"Glof4", "Lagf3" ,"Neof1", "Phof7" ,"Phof8", "Phof9", "Staf1", "Staf6",
                                  "Ttrf2" ,"Gn1" ,  "Mn1"  , "Pn1" ,  "Laga1", "Phoa1", "Phoa2" ,"Phoa3", 
                                  "Staj1" ,"Stan1" ,"Stln1" ,"Stln2")) 
#Order
raw_dat <- raw_dat[order(raw_dat$code),]

#Order ordination scores to match classifiers - same list as above
ordination_scores$code <- factor(ordination_scores$code, #check levels below from string printed with the code above
                                 levels = c("Bf1" ,  "Ff6"  , "Ff7" ,  "Hf8"  , "Hf9"  , "Sf2" ,  "Wf1",   "Kogf2" ,
                                            "Lagf1" ,"Lagf2", "Monf1", "Phof1", "Phof5", "Phof6", "Staf4" ,"Staf5", "Af3",
                                            "Af4" ,  "Af5"  , "Af6" ,  "Af7"  , "Af8"  , "Mf2" ,  "Ddef1", "Ddef2",
                                            "Glof3" ,"Glof4", "Lagf3" ,"Neof1", "Phof7" ,"Phof8", "Phof9", "Staf1", "Staf6",
                                            "Ttrf2" ,"Gn1" ,  "Mn1"  , "Pn1" ,  "Laga1", "Phoa1", "Phoa2" ,"Phoa3", 
                                            "Staj1" ,"Stan1" ,"Stln1" ,"Stln2"))  
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

##Save objects for later analysis
#Save specimens code as object
specimens <- row.names(raw_dat)

#Save growth categories as factor, useful for later analysis
categories <- classifiers$category

#Save groups as factor, useful for later analysis
groups <- classifiers$group

#Save taxa as factor, useful for later analysis
taxa <- classifiers$taxon

##Make natural log relevant columns
classifiers <- classifiers %>% mutate(BZW_log = log(BZW), bullaL_log  = log(bullaL), bullaW_log  = log(bullaW), 
                                      perioticL_log = log(perioticL), perioticW_log = log(perioticW))

##Create object with classifiers only for the best sampled taxa - useful later in analyses
#Create vector with the names of the best sampled taxa
#B.bonaerensis/acutorostata (Mysticeti), B.physalus (Mysticeti), Ph. phocoena (Odontoceti), St.attenuata (Odontoceti)
best_taxa <- c("B.acutorostrata", "B.bonaerensis", "B.physalus", "Ph.phocoena", "St.attenuata")

#Create new classifiers object
classifiers_taxa <- classifiers %>% filter(taxon %in% best_taxa)
#Replace to have only 1 taxon name for minkes
classifiers_taxa[classifiers_taxa == "B.acutorostrata"] <- "B.bonaerensis"
#Replace B.physalus with B.boanerensis - same genus, no B.bonaerensis for early fetus and needed for analysis of trajectory to have it
classifiers_taxa[classifiers_taxa == "B.physalus"] <- "B.bonaerensis"

##Save objects for later analysis
#Save growth categories as factor, useful for later analysis
categories_taxa <- classifiers_taxa$category

#Save groups as factor, useful for later analysis
groups_taxa <- classifiers_taxa$group

#Save taxa as factor, useful for later analysis
taxa_taxa <- classifiers_taxa$taxon

#Transform in 3D array, first number is number of landmarks, second is dimensions (3)
#To find the number of landmark points, divide the number of variables in raw_dat - visible in Environment - by 3
shape_array <- arrayspecs(raw_dat, 24894 , #number in raw_dat divided by 3
                          3) #number of dimensions

#Calculate mean shape coordinates
mean_shape <- mshape(shape_array) 

##Make data frame for analyses in geomorph
gdf <- geomorph.data.frame(coords = shape_array, 
                           code = classifiers$code, group = classifiers$group, 
                           category = classifiers$category, taxon = classifiers$taxon, BZW_log = classifiers$BZW_log,
                           bulla_log = classifiers$bullaL_log, periotic_log = classifiers$perioticL_log)
glimpse(gdf)

#Check for specimens that are outliers - if they are very young specimens it is expected
plotOutliers(gdf$coords)
#Save PDF/image to output folder using menu in bar on the right

##Make data frame for only best sampled taxa
#Search for row numbers for each taxon
best_taxa_rows <- c(
  which(classifiers[,4] == "B.acutorostrata"),
  which(classifiers[,4] == "B.bonaerensis"),
  which(classifiers[,4] == "B.physalus"),
  which(classifiers[,4] == "Ph.phocoena"),
  which(classifiers[,4] == "St.attenuata"))
#Display list of row numbers that do not contain these taxa in order
omit <- as.vector(setdiff(1:46, best_taxa_rows)) #1 to max number of observations from raw_dat

#Make new shape array eliminating rows of other taxa from shape array
shape_array_taxa <- shape_array[,,-omit]

#Calculate mean shape coordinates
mean_shape_taxa <- mshape(shape_array_taxa) 

#Make new gdf with classifiers for these taxa only
gdf_taxa <- geomorph.data.frame(coords = shape_array_taxa, code = classifiers_taxa$code, group = classifiers_taxa$group, 
                                category = classifiers_taxa$category, taxon = classifiers_taxa$taxon, BZW_log = classifiers_taxa$BZW_log,
                                bulla_log = classifiers_taxa$bullaL_log, periotic_log = classifiers_taxa$perioticL_log)
glimpse(gdf_taxa)

##Make palette with ggthemes-RColorBrewer - color and/or shapes
mypalette_blue <- as.matrix(ggthemes_data[["tableau"]][["color-palettes"]][["ordered-sequential"]][["Blue"]][["value"]])
image(1:20, 1, as.matrix(1:20), col = mypalette_blue, xlab = "Blue",
      ylab = "", xaxt = "n", yaxt = "n", bty = "n")

mypalette_paired <- brewer.pal(12,"Paired")
image(1:12, 1, as.matrix(1:12), col = mypalette_paired, xlab = "Paired",
      ylab = "", xaxt = "n", yaxt = "n", bty = "n")

#Palette for 4 best sampled taxa - B.bonaerensis/acutorostrata (Mysticeti), B.physalus (Mysticeti), Ph. phocoena (Odontoceti), St.attenuata (Odontoceti)
#same colors/taxa as growth and allometry project
mypalette_taxa <- c(mypalette_paired[4], mypalette_paired[3], mypalette_paired[9], mypalette_paired[10])
image(1:4, 1, as.matrix(1:4), col = mypalette_taxa, xlab = "taxa colors - B.bona, B.phys, Phoc., Sten.",
      ylab = "", xaxt = "n", yaxt = "n", bty = "n")

#Palette for categories - early, late, postnatal
mypalette_category <- c(mypalette_blue[3,], mypalette_blue[9,], mypalette_blue[16,])
mypalette_category_image <- image(1:3, 1, as.matrix(1:3), col = mypalette_category, xlab = "categories colors - early, late, postnatal", ylab = "", yaxt = "n")

#Create shape palette for groups
shapes <- c(15,19) #these are a square and a circle, use ?pch to see more shapes
#Create shape palette for best taxa
shapes_taxa <- c(15,19,17) #these are a square and a circle, use ?pch to see more shapes

##Images for plots
B.bonaerensis <- readPNG("Data/b.bona.png")
B.physalus <- readPNG("Data/b.physalus.png") #use for Mysticeti data only
Ph.phocoena <- readPNG("Data/phocoena.png")
St.attenuata <- readPNG("Data/stenella.png")

#CH. 3 - PCA COMPLETE DATASET ----

#Run PCA on complete dataset
PCA_all <- gm.prcomp(gdf$coords)

#List of PC components and proportion of variation
PCA_all 

#Sink list of PC component values to file and make it more readable in Word/Excel
sink("periotic_R/PCA_all_values.txt", append = F)
print(PCA_all)
sink()
#Word: replace all spaces as tabs, replace "Proportion of Variance" and Cumulative Proportion" as "Prop" and "Cum", save as new txt file
#Excel: check file read well as table, copy and paste transposed data, save
#Copy new edited file in Data folder after editing

#Re-import PC component values in more readable format to make columns like ordination values
PCA_all_values <- read.table("Data/PCA_all_values2.txt", header = T)
PCA_all_values

#Make new columns with % of variation for each axis and cumulative variation - useful for plots and if axes need to be excluded later
PCA_all_values <- PCA_all_values %>% 
  mutate(values_100 = (Prop*100)) %>% #percentage of variation described by each axis
  mutate(values_cum = (Cum*100)) #cumulative variation of that axis plus all the variation described by the previous ones
#Check
glimpse(PCA_all_values)

##Compare values between PCA and PCOORD (ordination) - make sure proportion of the first few axes is very similar
glimpse(ordination_values)
glimpse(PCA_all_values)

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
pcscores_all <- pcscores_all %>% mutate(specimens = gdf$code, taxon = gdf$taxon, group = gdf$group, category = gdf$category)
glimpse(pcscores_all)

#Nice plot with labels and specimens colored by category
PCA_all_ggplot <- ggplot(pcscores_all, aes(x = Comp1, y = Comp2, label = specimens, colour = category, shape = group))+
  geom_point(size = 3)+
  geom_text_repel(colour = "black", size = 3.5, max.overlaps = 60)+
  scale_colour_manual(name = "Growth category", labels =  c("Early Fetus", "Late Fetus", "Postnatal"), #to be ordered as they appear in tibble
                      values = mypalette_category)+            #legend and color adjustments
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("PC 1 (14.24%)")+ #copy this from standard PCA plot (PCA_all_plot)
  ylab("PC 2 (9.49%)")+
  ggtitle("PCA all data")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_ggplot

#Make hulls for PCA plot with hulls around growth category
hulls_all <- pcscores_all %>%
  group_by(category) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)
glimpse(hulls_all)

#Nice PCA plot with hulls around categories
PCA_all_hulls_ggplot <- ggplot(pcscores_all, aes(x = Comp1, y = Comp2, label = specimens, colour = category, fill = category))+
  geom_point(size = 3, aes(shape = group))+
  scale_colour_manual(name = "Growth category", labels =  c("Early Fetus", "Late Fetus", "Postnatal"), #to be ordered as they appear in tibble
                      values = mypalette_category)+            #legend and color adjustments
  geom_polygon(data = hulls_all, aes(x = x, y = y, fill = category), alpha = .5, show.legend = FALSE)+ #colored hulls with transparency
  scale_fill_manual(name = "Growth category", labels = c("Early Fetus", "Late Fetus", "Postnatal"),
                    values =  mypalette_category)+ #must match scale_colour_manual
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("PC 1 (14.24%)")+ #copy this from standard PCA plot (PCA_all_plot)
  ylab("PC 2 (9.49%)")+
  ggtitle("Periotic")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_hulls_ggplot  

##Regression PC1 and PC2 vs size (log periotic  length), group, category ----
#Check if size an important factor in results overall
#Create data frame with data
pcscores_all_size <- pcscores_all %>% mutate(size = gdf$periotic_log)

#Calculate regression for each component
reg_PC1all_size <- lm(Comp1 ~ size, data = pcscores_all_size)
reg_PC2all_size <- lm(Comp2 ~ size, data = pcscores_all_size)

#View results
summary(reg_PC1all_size)
summary(reg_PC2all_size)

#ANOVA tables
anova(reg_PC1all_size)
anova(reg_PC2all_size)

#Check if main components may be associated with group (Mysticeti/Odontoceti)
reg_PC1all_group <- lm(Comp1 ~ group, data = pcscores_all_size)
reg_PC2all_group <- lm(Comp2 ~ group, data = pcscores_all_size)

#View results
summary(reg_PC1all_group)
summary(reg_PC2all_group)

#ANOVA tables
anova(reg_PC1all_group)
anova(reg_PC2all_group)

#Check if main components may be associated with category (early, late, postnatal)
reg_PC1all_category <- lm(Comp1 ~ category, data = pcscores_all_size)
reg_PC2all_category <- lm(Comp2 ~ category, data = pcscores_all_size)

#View results
summary(reg_PC1all_category)
summary(reg_PC2all_category)

#ANOVA tables
anova(reg_PC1all_category)
anova(reg_PC2all_category)

#Save results of regressions to file
sink("periotic_R/PC1-PC2_lms.txt")
print("PC1 size")
summary(reg_PC1all_size)
anova(reg_PC1all_size)

print("PC2 size")
summary(reg_PC2all_size)
anova(reg_PC2all_size)

print("PC1 group")
summary(reg_PC1all_group)
anova(reg_PC1all_group)

print("PC2 group")
summary(reg_PC2all_group)
anova(reg_PC2all_group)

print("PC1 category")
summary(reg_PC1all_category)
anova(reg_PC1all_category)

print("PC2 category")
summary(reg_PC2all_category)
anova(reg_PC2all_category)
sink() 

#CH. 4 - PCOORD PLOT ----
#Plot the ordination scores from PCOORD imported from the GPSA

#Transform ordination scores as tibble
ordination_scores <- as_tibble(ordination_scores)
#Add labels and other attributes to tibble as columns
ordination_scores <- ordination_scores %>% mutate(specimens = gdf$code, taxon = gdf$taxon, group = gdf$group, category = gdf$category)
glimpse(ordination_scores)

#Recall values to copy for plot axis labels - use values_100 column for the first 2 axes
glimpse(ordination_values)

#Nice plot - axis 1 and 2
ordination_scores_ggplot <- ggplot(ordination_scores, aes(x = axis1, y = axis2, label = specimens, colour = category, shape = group))+
  geom_point(size = 3)+
  geom_text_repel(colour = "black", size = 3.5, max.overlaps = 50)+
  scale_colour_manual(name = "Growth category", labels =  c("Early Fetus", "Late Fetus", "Postnatal"), #to be ordered as they appear in tibble
                      values = mypalette_category)+            #legend and color adjustments
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("I (14.24%)")+ #copy this from ordination values_100 column printed before
  ylab("II (9.49%)")+
  ggtitle("PCOORD GPSA - Periotic")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12))  #title font and position

#Visualize plot and save as PDF using menu in bar on the right
ordination_scores_ggplot

##Regression Axis1 and Axis2 vs size (log periotic length), group, category ----
#Check if size an important factor in results overall
#Create data frame with data
ordination_scores_size <- ordination_scores %>% mutate(size = gdf$periotic_log)

#Calculate regression for each component
reg_axis1_size <- lm(axis1 ~ size, data = ordination_scores_size)
reg_axis2_size <- lm(axis2 ~ size, data = ordination_scores_size)

#View results
summary(reg_axis1_size)
summary(reg_axis2_size)

#ANOVA tables
anova(reg_axis1_size)
anova(reg_axis2_size)

#Check if main components may be associated with group (Mysticeti/Odontoceti)
reg_axis1_group <- lm(axis1 ~ group, data = ordination_scores_size)
reg_axis2_group <- lm(axis2 ~ group, data = ordination_scores_size)

#View results
summary(reg_axis1_group)
summary(reg_axis2_group)

#ANOVA tables
anova(reg_axis1_group)
anova(reg_axis2_group)

#Check if main components may be associated with category (early, late, postnatal)
reg_axis1_category <- lm(axis1 ~ category, data = ordination_scores_size)
reg_axis2_category <- lm(axis2 ~ category, data = ordination_scores_size)

#View results and p-value
summary(reg_axis1_category)
summary(reg_axis2_category)

#ANOVA tables
anova(reg_axis1_category)
anova(reg_axis2_category)

#Save results of regressions to file
sink("periotic_R/axis1-axis2_lms.txt")
print("I size")
summary(reg_axis1_size)
anova(reg_axis1_size)

print("II size")
summary(reg_axis2_size)
anova(reg_axis2_size)

print("I group")
summary(reg_axis1_group)
anova(reg_axis1_group)

print("II group")
summary(reg_axis2_group)
anova(reg_axis2_group)

print("I category")
summary(reg_axis1_category)
anova(reg_axis1_category)

print("II category")
summary(reg_axis2_category)
anova(reg_axis2_category)
sink() 

#CH. 5 - ANOVA OF SHAPE (pc scores) AND SIZE FOR TAXA - only well sampled taxa used for analysis  ----
#Conduct ANOVA to test if there is significant difference in allometry between each taxon
#Only well sampled taxa used for more even sampling

#Make 1 tibble to include only the selected taxa
pcscores_taxa <- pcscores_all %>% filter(taxon %in% best_taxa)
#Replace to have only 1 taxon name for minkes
pcscores_taxa[pcscores_taxa == "B.acutorostrata"] <- "B.bonaerensis"

#Create data frame of pc scores for well sampled taxa with size included
pcscores_taxa_size <- pcscores_taxa %>% mutate(size = gdf_taxa$periotic_log)

#Create basic model first with no interactions for comparison
allometry_pcscores_taxa <- procD.lm(pcscores_taxa_size[1:45] #number of PC components (1 per column)
                                    ~ pcscores_taxa_size$size, print.progress = FALSE, iter = 999)

#Check results
summary(allometry_pcscores_taxa)

##Pairwise comparison of regression models between taxa
#Test if there is a significant difference between allometry general model and allometry considering taxa separately

#Create 2 models, 1 with only intercept varying (combination) and one with slope and intercept varying (interaction)
allometry_pcscores_taxa_comb <- procD.lm(pcscores_taxa_size[1:45] ~ pcscores_taxa_size$size + pcscores_taxa_size$taxon,
                                         print.progress = FALSE, iter = 999)

allometry_pcscores_taxa_int <- procD.lm(pcscores_taxa_size[1:45] ~ pcscores_taxa_size$size * pcscores_taxa_size$taxon,
                                        print.progress = FALSE, iter = 999)

#Check results
summary(allometry_pcscores_taxa_comb)
summary(allometry_pcscores_taxa_int)

#Save results to file
sink("periotic_R/allometry_taxa.txt")
print("No group")
summary(allometry_pcscores_taxa)

print("Comb")
summary(allometry_pcscores_taxa_comb)

print("Int")
summary(allometry_pcscores_taxa_int)
sink() 

#Anova for difference between models - should find the comb and int to be significantly better than null and they ca be used in pairwise comparisons
anova(allometry_pcscores_taxa,allometry_pcscores_taxa_comb,allometry_pcscores_taxa_int)

#Pairwise comparison for the combination and interaction model
#Helps determine if there is a significant difference in slope (int model) in the allometry trajectory on top of difference in intercept (comb model)
pairwise_allometry_pcscores_taxa <- pairwise(allometry_pcscores_taxa_int, fit.null = allometry_pcscores_taxa_comb,
                                             groups = pcscores_taxa_size$taxon, 
                                             covariate =  pcscores_taxa_size$size, print.progress = FALSE) 
pairwise_allometry_pcscores_taxa

#Distances between slope vectors (end-points) - absolute difference between slopes of groups, if significant means int model better than comb
summary(pairwise_allometry_pcscores_taxa, confidence = 0.95, test.type = "dist") 

#Correlation between slope vectors (and angles) - similarity of vector orientation or angle,
#if significant means the vectors of the groups are oriented in different ways 
summary(pairwise_allometry_pcscores_taxa, confidence = 0.95, test.type = "VC",
        angle.type = "deg") 

#Save results to file
sink("periotic_R/pairwise_allometry_pcscores_taxa.txt")
print("ANOVA models")
anova(allometry_pcscores_taxa,allometry_pcscores_taxa_comb,allometry_pcscores_taxa_int)

print("1-Pairwise distances slopes")
summary(pairwise_allometry_pcscores_taxa, confidence = 0.95, test.type = "dist") 

print("2-Distance between angles")
summary(pairwise_allometry_pcscores_taxa, confidence = 0.95, test.type = "VC",
        angle.type = "deg") 


#Plot to obtain regression score of pc scores vs periotic length or periotic length, use best model - regression method with "RegScore" plotting
allometry_pcscores_taxa_plot <- plot(allometry_pcscores_taxa_int, type = "regression",predictor = pcscores_taxa_size$size, 
                                     reg.type = "RegScore", xlab = "Ln(periotic length)", 
                                     pch = 21, col = "chartreuse4", bg = "chartreuse4", cex = 1.2, font.main = 2)   #improve graphics


##Make better allometry plot with ggplot
#Create data frame object that ggplot can read - use data from plot object you want to improve
allometry_pcscores_taxa_tibble <- data.frame(size = allometry_pcscores_taxa_plot[["plot.args"]][["x"]], 
                                             RegScores = allometry_pcscores_taxa_plot[["plot.args"]][["y"]])

#Convert data frame to tibble
allometry_pcscores_taxa_tibble<- as_tibble(allometry_pcscores_taxa_tibble)
#Add labels and other attributes to tibble as columns
allometry_pcscores_taxa_tibble <- allometry_pcscores_taxa_tibble %>% mutate(specimens = gdf_taxa$code, taxon = gdf_taxa$taxon, 
                                                                            group = gdf_taxa$group, category = gdf_taxa$category)
glimpse(allometry_pcscores_taxa_tibble)


#Nice plot with specimens colored by group AND regression lines for each group
allometry_pcscores_taxa_ggplot <- ggplot(allometry_pcscores_taxa_tibble, aes(x = size, y = RegScores, colour = taxon, shape = group))+
  geom_smooth(aes(x = size, y = RegScores, colour =  taxon, fill = taxon, linetype = group), method = 'lm', inherit.aes = F,         #confidence intervals and reg line, before points
              alpha = 0.2, size = 1, show.legend = F)+ #should be straight regression line with confidence interval in grey
  geom_point(size = 3, alpha = 0.3)+       #points after, so they are on top
  scale_color_manual(values = c(mypalette_taxa[1], mypalette_taxa[3:4]), #select colors from palette from taxa
                     aesthetics = c("color","fill"))+          
  scale_shape_manual(values = shapes)+
  theme_classic(base_size = 12)+
  xlab("Log(Periotic length)")+
  ylab("Regression Score**")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5), 
        legend.position = "none", legend.direction = "horizontal",
        axis.title.x = element_text(face = "bold", size = 14))
#Add silhouettes groups
allometry_pcscores_taxa_ggplot <- allometry_pcscores_taxa_ggplot   + 
  add_phylopic(B.bonaerensis, alpha = 1, x = 3.8, y = -10, ysize = 3, color = mypalette_taxa[1])+
  add_phylopic(Ph.phocoena, alpha = 1, x = 2.8, y = -6, ysize = 3.4, color = mypalette_taxa[3])+
  add_phylopic(St.attenuata, alpha = 1, x = 2.3, y = 4, ysize = 2.7, color = mypalette_taxa[4])
#Visualize plot and save as PDF using menu in bar on the right
allometry_pcscores_taxa_ggplot

#CH. 6 - TRAJECTORY ANALYSIS - entire dataset used for analysis ----
#Shows trajectories of shape variation, possible to compare trajectories of growth between odontocetes and mysticetes
#Use original shape data, entire dataset to account for uneven sampling

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
sink("periotic_R/group_trajectory.txt")
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
  xlab("PC 1 (44.32%)")+ #copy this from standard trajectory plot
  ylab("PC 2 (24.62%)")+
  ggtitle("Periotic")+
  guides(shape = guide_legend(label = F, title = NULL, override.aes = list(shape = NA)), 
         colour = guide_legend(label = F, title = NULL))+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16), legend.position = "bottom", legend.direction = "horizontal")
#Add silhouettes groups
group_trajectory_ggplot <- group_trajectory_ggplot   + 
  add_phylopic(B.physalus, alpha = 1, x = -8, y = 5, ysize = 1.6, color = mypalette_taxa[2])+
  add_phylopic(St.attenuata, alpha = 1, x = 2.8, y = -2.5, ysize = 1.5, color = mypalette_taxa[4])
#Visualize plot and save as PDF using menu in bar on the right
group_trajectory_ggplot
