
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
library(rphylopic)
library(png)
library(gridExtra)


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
ordination_values <- read.table("Data/_gpsa_ordination_values_bulla_R.RData", header = T)
#Check
glimpse(ordination_values)

#Make new columns with % of variation for each axis and cumulative variation - useful for plots and if axes need to be excluded later
ordination_values <- ordination_values %>% 
  mutate(values_100 = (pro*100)) %>% #percentage of variation described by each axis
  mutate(values_cum = (cum*100)) #cumulative variation of that axis plus all the variation described by the previous ones
#Check
glimpse(ordination_values)

#Sink values to file in more readable format
sink("bulla_R/ordination_values.txt", append = F)
print(ordination_values )
sink()

##Order classifiers by category, useful for plot legend
#Make factor for variable
classifiers$category <- factor(classifiers$category, 
                            levels = c("early", "late", "born")) 
#Order
classifiers <- classifiers[order(classifiers$category),]
#Check
glimpse(classifiers)

#Order shape data by category to match classifiers
classifiers$code

#Make factor for variable
raw_dat$code <- factor(raw_dat$code, #check levels below from string printed with the code above
                       levels = c("Be3" ,  "Bf1"  , "Ff1"  , "Ff3" ,  "Ff4" ,  "Ff5"  , "Ff6" ,  "Ff7" ,  "Gf1"  , "Hf3" ,  "Hf4"  , "Hf6"  ,
                                  "Hf8",   "Hf9" ,  "Mf1",   "Sf1"  , "Sf2"  , "Wf1" ,  "Delf1",
                                  "Delf2" ,"Delf3", "Gloe1" ,"Glof1", "Glof2" ,"Kogf1" ,"Kogf2" ,"Lagf1", "Lagf2", "Monf1", "Monf2", "Phof1" ,
                                  "Phof3" ,"Phof4", "Phof5", "Phof6", "Phye1" ,"Phye2", "Phye3",
                                  "Phyf2", "Psef1", "Staf2", "Staf3", "Staf4", "Staf5", "Ttrf1" ,"Af3" ,  "Af4" ,  "Af5" ,  "Af6" , "Af7" ,  
                                  "Af8"  , "Mf2",   "Ddef1" ,"Ddef2", "Glof3" ,"Glof4" ,"Lagf3",
                                  "Neof1", "Phof7", "Phof8" ,"Phof9", "Staf1", "Staf6" ,"Ttrf2" ,"Gn1" ,  "Mn1" ,  "Pn1",   "Laga1" ,"Phoa1", 
                                  "Phoa2" ,"Phoa3", "Staj1", "Stan1" ,"Stln1", "Stln2")) 
#Order
raw_dat <- raw_dat[order(raw_dat$code),]

#Order ordination scores to match classifiers - same list as above
ordination_scores$code <- factor(ordination_scores$code, #check levels below from string printed with the code above
                                 levels = c("Be3" ,  "Bf1"  , "Ff1"  , "Ff3" ,  "Ff4" ,  "Ff5"  , "Ff6" ,  "Ff7" ,  "Gf1"  , "Hf3" ,  "Hf4"  , "Hf6"  ,
                                            "Hf8",   "Hf9" ,  "Mf1",   "Sf1"  , "Sf2"  , "Wf1" ,  "Delf1",
                                            "Delf2" ,"Delf3", "Gloe1" ,"Glof1", "Glof2" ,"Kogf1" ,"Kogf2" ,"Lagf1", "Lagf2", "Monf1", "Monf2", "Phof1" ,
                                            "Phof3" ,"Phof4", "Phof5", "Phof6", "Phye1" ,"Phye2", "Phye3",
                                            "Phyf2", "Psef1", "Staf2", "Staf3", "Staf4", "Staf5", "Ttrf1" ,"Af3" ,  "Af4" ,  "Af5" ,  "Af6" , "Af7" ,  
                                            "Af8"  , "Mf2",   "Ddef1" ,"Ddef2", "Glof3" ,"Glof4" ,"Lagf3",
                                            "Neof1", "Phof7", "Phof8" ,"Phof9", "Staf1", "Staf6" ,"Ttrf2" ,"Gn1" ,  "Mn1" ,  "Pn1",   "Laga1" ,"Phoa1", 
                                            "Phoa2" ,"Phoa3", "Staj1", "Stan1" ,"Stln1", "Stln2")) 
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

## save objects for later analysis
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
#B.bonaerensis (Mysticeti), B.physalus (Mysticeti), Ph. phocoena (Odontoceti), St.attenuata (Odontoceti)
best_taxa <- c("B.acutorostrata", "B.bonaerensis", "B.physalus", "Ph.phocoena", "St.attenuata")

#Create new classifiers object
classifiers_taxa <- classifiers %>% filter(taxon %in% best_taxa)
#Replace to have only 1 taxon name for minkes
classifiers_taxa[classifiers_taxa == "B.acutorostrata"] <- "B.bonaerensis"

#Transform in 3D array, first number is number of landmarks, second is dimensions (3)
#To find the number of landmark points, divide the number of variables in raw_dat - visible in Environment - by 3
shape_array <- arrayspecs(raw_dat, 24953 , #number in raw_dat divided by 3
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
omit <- as.vector(setdiff(1:75, best_taxa_rows)) 

#Make new shape array eliminating rows of other taxa from shape array
shape_array_taxa <- shape_array[,,-omit]

#Calculate mean shape coordinates
mean_shape_taxa <- mshape(shape_array_taxa) 

#Make new gdf with classifiers for these taxa only
gdf_taxa <- geomorph.data.frame(coords = shape_array_taxa, code = classifiers_taxa$code, group = classifiers_taxa$group, 
                                category = classifiers_taxa$category, taxon = classifiers_taxa$taxon,
                                bulla_log = classifiers_taxa$bullaL_log, periotic_log = classifiers_taxa$perioticL_log)
glimpse(gdf_taxa)


##Make palette with ggthemes - color and/or shapes
#Palettes from ggthemes_data
mypalette_blue <- as.matrix(ggthemes_data[["tableau"]][["color-palettes"]][["ordered-sequential"]][["Blue"]][["value"]])
image(1:20, 1, as.matrix(1:20), col = mypalette_blue, xlab = "Blue",
      ylab = "", xaxt = "n", yaxt = "n", bty = "n")
mypalette_tableau20 <- as.matrix(ggthemes_data[["tableau"]][["color-palettes"]][["regular"]][["Tableau 20"]][["value"]])
image(1:20, 1, as.matrix(1:20), col = mypalette_tableau20, xlab = "Tableau20",
      ylab = "", xaxt = "n", yaxt = "n", bty = "n")

#Palette for 4 best sampled taxa - B.bonaerensis (Mysticeti), B.physalus (Mysticeti), Ph. phocoena (Odontoceti), St.attenuata (Odontoceti)
#same colors/taxa as growth and allometry project
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

mypalette_taxa <- c(mypalette_Mysticeti[2], mypalette_Mysticeti[3], mypalette_Odontoceti[1], mypalette_Odontoceti[6])
mypalette_taxa_image <- image(1:4, 1, as.matrix(1:4), col = mypalette_taxa, xlab = "taxa colors - B.bona, B.phys, Phoc., Sten.",
                                  ylab = "", xaxt = "n", yaxt = "n", bty = "n")

#Palette for categories - early, late, born
mypalette_category <- c(mypalette_blue[3,], mypalette_blue[9,], mypalette_blue[16,])
mypalette_category_image <- image(1:3, 1, as.matrix(1:3), col = mypalette_category, xlab = "categories colors - early, late, born", ylab = "", yaxt = "n")

#Create shape palette for groups
shapes <- c(15,19) #these are a square and a circle, use ?pch to see more shapes

##Images for plots
B.bonaerensis <- readPNG("Data/b.bona.png")
B.physalus <- readPNG("Data/b.physalus.png")
Ph.phocoena <- readPNG("Data/phocoena.png")
St.attenuata <- readPNG("Data/stenella.png")

#CH. 3 - PCA COMPLETE DATASET ----

#Run PCA on complete dataset
PCA_all <- gm.prcomp(gdf$coords)

#List of PC components and proportion of variation
PCA_all 

#Sink list of PC component values to file and make it more readable in Word/Excel
sink("bulla_R/PCA_all_values.txt", append = F)
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

#Nice PCA plot with category and groups
PCA_all_ggplot <- ggplot(pcscores_all, aes(x = Comp1, y = Comp2, label = specimens, colour = category, shape = group))+
  geom_point(size = 3)+
  geom_text_repel(colour = "black", size = 3.5, max.overlaps = 60)+
  scale_colour_manual(name = "Growth category", labels =  c("Early Fetus", "Late Fetus", "Postnatal"), #to be ordered as they appear in tibble
                      values = mypalette_category)+            #legend and color adjustments
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("PC 1 (15.79%)")+ #copy this from standard PCA plot (PCA_all_plot)
  ylab("PC 2 (7.14%)")+
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
  xlab("PC 1 (15.79%)")+ #copy this from standard PCA plot (PCA_all_plot)
  ylab("PC 2 (7.14%)")+
  ggtitle("PCA all data")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_hulls_ggplot  

##Plots for each category ----
#PCA plots for each category
#Create one tibble for each category
#First split between categories
pcscores_all_category <- pcscores_all %>% group_by(category) %>% group_split()
#Check order
View(pcscores_all_category)

#Save as separate tibbles - 1 for early fetus, 1 for late fetus, 1 for all 3 postnatal category (neonate, juvenile, adult)
pcscores_all_earlyfetus <- pcscores_all_category[[1]]
pcscores_all_latefetus <- pcscores_all_category[[2]]
pcscores_all_postnatal <- pcscores_all_category[[3]]

#Nice PCA plot with groups early fetus
PCA_all_earlyfetus_ggplot <- ggplot(pcscores_all_earlyfetus, aes(x = Comp1, y = Comp2, label = specimens, shape = group))+
  geom_point(size = 3, color = mypalette_category[1])+
  geom_text_repel(colour = "black", size = 3.5, max.overlaps = 40)+
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("PC 1 (15.79%)")+ #copy this from standard PCA plot (PCA_all_plot)
  ylab("PC 2 (7.14%)")+
  ggtitle("PCA early fetus")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_earlyfetus_ggplot

#Nice PCA plot with groups late fetus
PCA_all_latefetus_ggplot <- ggplot(pcscores_all_latefetus, aes(x = Comp1, y = Comp2, label = specimens, shape = group))+
  geom_point(size = 3, color = mypalette_category[2])+
  geom_text_repel(colour = "black", size = 3.5, max.overlaps = 40)+
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("PC 1 (15.79%)")+ #copy this from standard PCA plot (PCA_all_plot)
  ylab("PC 2 (7.14%)")+
  ggtitle("PCA late fetus")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_latefetus_ggplot

#Nice PCA plot with groups postnatal
PCA_all_postnatal_ggplot <- ggplot(pcscores_all_postnatal, aes(x = Comp1, y = Comp2, label = specimens, shape = group))+
  geom_point(size = 3, color = mypalette_category[3])+
  geom_text_repel(colour = "black", size = 3.5, max.overlaps = 40)+
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("PC 1 (15.79%)")+ #copy this from standard PCA plot (PCA_all_plot)
  ylab("PC 2 (7.14%)")+
  ggtitle("PCA postnatal")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_postnatal_ggplot

##Plots for well sampled taxa only ----
##PCA plot including only the 4 well sampled taxa
#Create one tibble with the 4 taxa
#First split between all taxa
pcscores_all_taxon <- pcscores_all %>% group_by(taxon) %>% group_split()
#Check order
levels(as.factor(pcscores_all$taxon))

#Save as separate tibbles - 1 for each taxon
pcscores_all_B.acutorostrata <- pcscores_all_taxon[[1]]
pcscores_all_B.bonaerensis <- pcscores_all_taxon[[2]]
pcscores_all_B.physalus <- pcscores_all_taxon[[5]]
pcscores_all_Ph.phocoena <- pcscores_all_taxon[[17]]
pcscores_all_St.attenuata <- pcscores_all_taxon[[20]]
#Combine 2 minkes in 1 tibble
pcscores_all_minke <- bind_rows(pcscores_all_B.acutorostrata, pcscores_all_B.bonaerensis)
#Replace to have only 1 taxon name
pcscores_all_minke[pcscores_all_minke == "B.acutorostrata"] <- "B.bonaerensis"
glimpse(pcscores_all_minke)

#Make 1 tibble to include only the selected taxa
pcscores_taxa <- bind_rows(pcscores_all_minke, pcscores_all_B.physalus, pcscores_all_Ph.phocoena, pcscores_all_St.attenuata)

#Nice PCA plot with groups early fetus
PCA_taxa_ggplot <- ggplot(pcscores_taxa, aes(x = Comp1, y = Comp2, shape = group, colour = taxon, alpha = category, label = specimens))+
  geom_point(size = 3)+
  geom_text_repel(colour = "black", size = 3.5, max.overlaps = 60, show.legend=FALSE)+
  scale_colour_manual(name = "Taxa", labels =  c("B.bonaerensis", "B.physalus", "Ph.phocoena", "St.attenuata"), #to be ordered as they appear in tibble
                      values = mypalette_taxa)+            #legend and color adjustments
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  scale_alpha_manual(name = "Growth category", labels =  c("Early Fetus", "Late Fetus", "Postnatal"), values = c(0.4, 0.7, 1))+
  theme_bw()+
  xlab("PC 1 (15.79%)")+ #copy this from standard PCA plot (PCA_all_plot)
  ylab("PC 2 (7.14%)")+
  ggtitle("PCA selected taxa")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_taxa_ggplot

#Make hulls for PCA plot with hulls around taxa
hulls_taxa <- pcscores_taxa %>%
  group_by(taxon) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)
glimpse(hulls_taxa)

#Nice PCA plot with hulls around categories
PCA_taxa_hulls_ggplot <- ggplot(pcscores_taxa, aes(x = Comp1, y = Comp2, colour = taxon, alpha = category))+
  geom_point(size = 3, aes(shape = group))+
  scale_colour_manual(name = "Taxa", labels =  c("B.bonaerensis", "B.physalus", "Ph.phocoena", "St.attenuata"), #to be ordered as they appear in tibble
                      values = mypalette_taxa)+            #legend and color adjustments
  scale_alpha_manual(name = "Growth category", labels =  c("Early Fetus", "Late Fetus", "Postnatal"), values = c(0.4, 0.7, 1))+
  geom_polygon(data = hulls_taxa, aes(x = x, y = y, fill = taxon), alpha = .2, show.legend = FALSE)+ #colored hulls with transparency
  scale_fill_manual(name = "Taxa", labels =  c("B.bonaerensis", "B.physalus", "Ph.phocoena", "St.attenuata"), #to be ordered as they appear in tibble
                    values = mypalette_taxa)+ #must match scale_colour_manual
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_bw()+
  xlab("PC 1 (15.79%)")+ #copy this from standard PCA plot (PCA_all_plot)
  ylab("PC 2 (7.14%)")+
  ggtitle("PCA selected taxa")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

#Visualize plot and save as PDF using menu in bar on the right
PCA_taxa_hulls_ggplot  

##Regression PC1 and PC2 vs bulla length ----
#Check if size an important factor in results overall
#Create data frame with data
pcscores_all_size <- pcscores_all %>% mutate(size = gdf$bulla_log)

#Calculate regression for each component
reg_PC1all_size <- lm(Comp1 ~ size, data = pcscores_all_size)
reg_PC2all_size <- lm(Comp2 ~ size, data = pcscores_all_size)

#View results and p-value
summary(reg_PC1all_size)
summary(reg_PC2all_size)

#Save results of significant regression to file
sink("bulla_R/PC1-PC2_size_lm.txt")
print("PC1")
summary(reg_PC1all_size)

print("PC2")
summary(reg_PC2all_size)
sink() 

#CH. 4 - ALLOMETRY ----
##Evaluate allometry and get the allometry-free shapes using log bulla length

##Allometry entire daatset ----
#Regression shape on size
allometry <- procD.lm(coords ~ bulla_log, data = gdf, iter=999, print.progress = TRUE) 

#Check results and p-value for significance
summary(allometry)

#Save results of significant regression to file
sink("bulla_R/allometry_all.txt")
summary(allometry)
sink() 

#Create residuals array to then save as coordinates for analyses
allometry_array <- arrayspecs(allometry$residuals,p = dim(gdf$coords)[1], k = dim(gdf$coords)[2]) 

#New shapes adjusted for allometry with CS to use in analyses
allometry_residuals <- allometry_array + array(mean_shape, dim(allometry_array)) 

#Save mean shape of allometry-adjusted shapes to use in a analyses
mean_shape_residuals <- mshape(allometry_residuals)

#Regression score of shape vs bulla length or periotic length  - regression method with "RegScore" plotting
allometry_plot_regscore <- plot(allometry, type = "regression",predictor = gdf$bulla_log, reg.type = "RegScore",
                                main = "Shape vs Ln(bulla length)",xlab = "Ln(bulla length)", 
                                pch = 21, col = "chartreuse4", bg = "chartreuse4", cex = 1.2, font.main = 2)   #improve graphics

##Add regression line with confidence intervals to plot
#Create object to use for linear model
allometry_regscores <- allometry_plot_regscore[["RegScore"]] 

##Make better allometry plot with ggplot
#Create data frame object that ggplot can read - use data from plot object you want to improve
allometry_plot <- data.frame(size = allometry_plot_regscore[["plot.args"]][["x"]], RegScores = allometry_plot_regscore[["plot.args"]][["y"]])

#Convert data frame to tibble
allometry_plot <- as_tibble(allometry_plot)
#Add labels and other attributes to tibble as columns
allometry_plot <- allometry_plot %>% mutate(specimens = gdf$code, taxon = gdf$taxon, group = gdf$group, category = gdf$category)
glimpse(allometry_plot)

#Nice plot with specimens colored by age AND regression line with confidence intervals
allometry_ggplot <- ggplot(allometry_plot, aes(x = size, y = RegScores, colour = category, shape = group))+
  geom_smooth(aes(x = size, y = RegScores), method = 'lm', inherit.aes = F,         #confidence intervals and reg line, before points
              colour = "darkblue", fill = 'gainsboro', linetype = "dashed", size = 0.5)+ #should be straight regression line with confidence interval in grey
  geom_point(size = 3)+       #points after, so they are on top
  scale_colour_manual(name = "Growth category", labels =  c("Early Fetus", "Late Fetus", "Postnatal"), #to be ordered as they appear in tibble
                      values = mypalette_category)+            #legend and color adjustments          
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_classic(base_size = 12)+
  xlab("Log(Bulla length)")+
  ylab("Regression Score")+
  ggtitle ("Allometry plot - p-value = 0.001**")+  #copy from model summary
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

#Visualize plot and save as PDF using menu in bar on the right
allometry_ggplot

##Allometry well sampled taxa only ----
#Regression shape on size
allometry_taxa <- procD.lm(coords ~ bulla_log, data = gdf_taxa, iter=999, print.progress = TRUE) 

#Check results and p-value for significance
summary(allometry_taxa)

#Save results of significant regression to file
sink("bulla_R/allometry_taxa.txt")
summary(allometry_taxa)
sink() 

#Create residuals array to then save as coordinates for analyses
allometry_array_taxa <- arrayspecs(allometry_taxa$residuals,p = dim(gdf_taxa$coords)[1], k = dim(gdf_taxa$coords)[2]) 

#New shapes adjusted for allometry with CS to use in analyses
allometry_residuals_taxa <- allometry_array_taxa + array(mean_shape_taxa, dim(allometry_array_taxa)) 

#Save mean shape of allometry-adjusted shapes to use in a analyses
mean_shape_residuals <- mshape(allometry_residuals)

#Regression score of shape vs bulla length or periotic length  - regression method with "RegScore" plotting
allometry_plot_regscore_taxa <- plot(allometry_taxa, type = "regression",predictor = gdf_taxa$bulla_log, reg.type = "RegScore",
                                main = "Shape vs Ln(bulla length)",xlab = "Ln(bulla length)", 
                                pch = 21, col = "chartreuse4", bg = "chartreuse4", cex = 1.2, font.main = 2)   #improve graphics

##Add regression line with confidence intervals to plot
#Create object to use for linear model
allometry_regscores_taxa <- allometry_plot_regscore_taxa[["RegScore"]] 

##Make better allometry plot with ggplot
#Create data frame object that ggplot can read - use data from plot object you want to improve
allometry_plot_taxa <- data.frame(size = allometry_plot_regscore_taxa[["plot.args"]][["x"]], 
                                  RegScores = allometry_plot_regscore_taxa[["plot.args"]][["y"]])

#Convert data frame to tibble
allometry_plot_taxa <- as_tibble(allometry_plot_taxa)
#Add labels and other attributes to tibble as columns
allometry_plot_taxa <- allometry_plot_taxa %>% mutate(specimens = gdf_taxa$code, taxon = gdf_taxa$taxon, 
                                                      group = gdf_taxa$group, category = gdf_taxa$category)
glimpse(allometry_plot_taxa)

#Nice plot with specimens colored by age AND regression line with confidence intervals
allometry_taxa_ggplot <- ggplot(allometry_plot_taxa, aes(x = size, y = RegScores, colour = category, shape = group))+
  geom_smooth(aes(x = size, y = RegScores), method = 'lm', inherit.aes = F,         #confidence intervals and reg line, before points
              colour = "darkblue", fill = 'gainsboro', linetype = "dashed", size = 0.5)+ #should be straight regression line with confidence interval in grey
  geom_point(size = 3)+       #points after, so they are on top
  scale_colour_manual(name = "Growth category", labels =  c("Early Fetus", "Late Fetus", "Postnatal"), #to be ordered as they appear in tibble
                      values = mypalette_category)+            #legend and color adjustments          
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_classic(base_size = 12)+
  xlab("Log(Bulla length)")+
  ylab("Regression Score")+
  ggtitle ("Allometry plot - p-value = 0.001**")+  #copy from model summary
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

#Visualize plot and save as PDF using menu in bar on the right
allometry_taxa_ggplot

#Nice plot with specimens colored by taxon AND regression lines for each taxon with confidence intervals
allometry_taxa_ggplot_taxon <- ggplot(allometry_plot_taxa, aes(x = size, y = RegScores, colour = taxon, shape = group))+
  geom_smooth(aes(x = size, y = RegScores, colour =  taxon, fill = taxon, linetype = group), method = 'lm', inherit.aes = F,         #confidence intervals and reg line, before points
              alpha = 0.3, size = 1, show.legend = F)+ #should be straight regression line with confidence interval in grey
  geom_point(size = 3, alpha = 0.3)+       #points after, so they are on top
  scale_color_manual(name = "Taxa", labels  = c("B.bonaerensis", "B.physalus", "Ph.phocoena", "St.attenuata"), values = c(mypalette_taxa), #select colors from palette from taxa
                     aesthetics = c("color","fill"))+          
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_classic(base_size = 12)+
  xlab("Log(Bulla length)")+
  ylab("Regression Score")+
  ggtitle ("Allometry plot - p-value = 0.001**")+  #copy from model summary
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

#Visualize plot and save as PDF using menu in bar on the right
allometry_taxa_ggplot_taxon

#Nice plot with specimens colored by taxon AND regression lines for each taxon with confidence intervals no legend but phylopic
allometry_taxa_ggplot_taxon_nolegend <- ggplot(allometry_plot_taxa, aes(x = size, y = RegScores, colour = taxon, shape = group))+
  geom_smooth(aes(x = size, y = RegScores, colour =  taxon, fill = taxon, linetype = group), method = 'lm', inherit.aes = F,         #confidence intervals and reg line, before points
              alpha = 0.3, size = 1, show.legend = F)+ #should be straight regression line with confidence interval in grey
  geom_point(size = 3, alpha = 0.3)+       #points after, so they are on top
  scale_color_manual(name = "Taxa", labels  = c("B.bonaerensis", "B.physalus", "Ph.phocoena", "St.attenuata"), values = c(mypalette_taxa), #select colors from palette from taxa
                     aesthetics = c("color","fill"))+          
  scale_shape_manual(name = "Group", labels = c("Mysticeti", "Odontoceti"), values = shapes)+
  theme_classic(base_size = 12)+
  xlab("Log(Bulla length)")+
  ylab("Regression Score")+
  ggtitle ("Allometry plot - p-value = 0.001**")+  #copy from model summary
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "none")
#Add silhouettes taxa
allometry_taxa_ggplot_taxon_nolegend <- allometry_taxa_ggplot_taxon_nolegend  + 
  add_phylopic(B.bonaerensis, alpha = 1, x = 3.9, y = -4, ysize = 2.6, color = mypalette_taxa[1])+
  add_phylopic(B.physalus, alpha = 1, x = 3.3, y = 9.5, ysize = 2.6, color = mypalette_taxa[2])+
  add_phylopic(Ph.phocoena, alpha = 1, x = 2.2, y = 0, ysize = 3.2, color = mypalette_taxa[3])+
  add_phylopic(St.attenuata, alpha = 1, x = 2.3, y = -15, ysize = 2.7, color = mypalette_taxa[4])

#Visualize plot and save as PDF using menu in bar on the right
allometry_taxa_ggplot_taxon_nolegend


#CH. 5 - ORDINATION PLOT ----

#Plot the ordination scores imported from the GPSA
#Transform ordination scores as tibble
ordination_scores <- as_tibble(ordination_scores)
#Add labels and other attributes to tibble as columns
ordination_scores <- ordination_scores %>% mutate(specimens = gdf$code, taxon = gdf$taxon, group = gdf$group)
glimpse(ordination_scores)

#Recall values to copy for plot axis labels - use values_100 column for the first 2 axes
ordination_values

#Nice plot - axis 1 and 2
ordination_scores_ggplot <- ggplot(ordination_scores, aes(x = axis1, y = axis2, label = specimens, colour = category, shape = group))+
  geom_point(size = 3)+
  geom_text_repel(colour = "black", size = 3.5, max.overlaps = 40)+
  scale_colour_manual(name = "Growth category", labels =  c("Early Fetus", "Late Fetus", "Neonate", "Juvenile", "Adult"), #to be ordered as they appear in tibble
                      values = mypalette_category)+            #legend and color adjustments
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
allometry_group_ggplot <- ggplot(allometry_plot, aes(x = size, y = RegScores, label = specimens, colour = category, shape = group))+
  geom_smooth(aes(x = size, y = RegScores, linetype = group), method = 'lm', inherit.aes = F,         #confidence intervals and reg line, before points
              colour = "darkblue", fill = 'gainsboro', size = 0.5, show.legend = F)+ #should be straight regression line with confidence interval in grey
  geom_point(size = 3)+       #points after, so they are on top
  scale_colour_manual(name = "Growth category", labels =  c("Early Fetus", "Late Fetus", "Neonate", "Juvenile", "Adult"), #to be ordered as they appear in tibble
                      values = mypalette_category)+            #legend and color adjustments          
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
#Are shapes in each group different from the other groups at different growth category?
#2 analyses, use shape data and pc scores
#Shapes
#Make 3 data frames, 1 for each category - the numbers are the rows corresponding to specimens at that growth category
gdf_earlyfetus <- geomorph.data.frame(coords = shape_array[,,1:15], 
                                      code = classifiers$code[1:15], group = classifiers$group[1:15], 
                                      category = classifiers$category[1:15], category = classifiers$category[1:15],
                                      bulla_log = classifiers$bullaL_log[1:15], periotic_log = classifiers$perioticL_log[1:15])
glimpse(gdf_earlyfetus)

gdf_latefetus <- geomorph.data.frame(coords = shape_array[,,16:27], 
                                     code = classifiers$code[16:27], group = classifiers$group[16:27], 
                                     category = classifiers$category[16:27], category = classifiers$category[16:27],
                                     bulla_log = classifiers$bullaL_log[16:27], periotic_log = classifiers$perioticL_log[16:27])
glimpse(gdf_latefetus)

gdf_postnatal <- geomorph.data.frame(coords = shape_array[,,28:38], 
                                     code = classifiers$code[28:38], group = classifiers$group[28:38], 
                                     category = classifiers$category[28:38], category = classifiers$category[28:38],
                                     bulla_log = classifiers$bullaL_log[28:38], periotic_log = classifiers$perioticL_log[28:38])
glimpse(gdf_postnatal)

#Make separate dataframe to check shape is not different postnatal in odontocetes
gdf_postnatal_odont <- geomorph.data.frame(coords = shape_array[,,31:38], 
                                           code = classifiers$code[31:38], group = classifiers$group[31:38], 
                                           category = classifiers$category[31:38], category = classifiers$category[31:38],
                                           bulla_log = classifiers$bullaL_log[31:38], periotic_log = classifiers$perioticL_log[31:38])
glimpse(gdf_postnatal_odont)

#Start by making sure it is ok to lump postnatal category with ANOVA of shape per category in postnatal odontocetes
shape_category_odont_anova <- procD.lm(coords ~ category, iter = 999, data = gdf_postnatal_odont, RRPP = F)

#Check results - should not be significant
summary(shape_category_odont_anova)

#Save results to file
sink("Output/shape_category_odont_anova.txt")
summary(shape_category_odont_anova)
anova(shape_category_odont_anova)
sink() 

#Proceed with shape differences between groups for each category
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
#Repeat same analyses but using pc scores for each category
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

