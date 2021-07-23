# Developing echolocation: distinctive patterns in the ontogeny of the tympanoperiotic complex in toothed and baleen whales (Cetacea) 	:dolphin: :whale:
Analyses of growth allometry and shape change during ontogeny of the bulla and periotic of Cetacea 

Authors: [Agnese Lanzetti](mailto:agnese.lanzetti@gmail.com?subject=[GitHub]%20Earbones%20Paper%20Code)

To cite the paper:


Available at: URL repo

If using any of this code or data please cite the paper above and this repo

To cite this repo:


DOI ZENODO

## Data :floppy_disk: 

The data are provided in the Data folder.

- Data for size allometry analysis: *earbones_ossification_events.csv, growth_data.csv, measuraments.csv* 
Specimen codes are the same as the ones listed in table S1, where additonal details on the specimens are provided. For complete references for the growth data see the publication.

- Data for shape analysis of the tympanic bulla: *_gpsa_homologized_points_bulla_R.dat, _gpsa_ordination_projections_bulla_R.Rdata, _gpsa_ordination_values_bulla_R.RData, bulla_classifiers.csv*
Homologized points file contains the coordinates for the aligned points as produced by the GPSA. Ordination projections ar ethe PCOORD valeu also as produced by the GPSA. PCOORD oridination values were copied from the GPSA print out during analysis and the file produced was updated to include also the proportaional and cumulative values. Classifier files contains group, taxon, onotgentic category an measuraments of leght for the bulla and periotic for all specimens. For spciemns codes and numbers see the publiciation and table S1.

- Data for shape analysis of the periotic: *_gpsa_homologized_points_periotic_R.dat, _gpsa_ordination_projections_periotic_R.Rdata, _gpsa_ordination_values_periotic_R.RData, periotic_classifiers.csv*
Same as for the tympanic bulla.

- Silhouettes of taxa and groups for plots: *b.bona.png, b.physalus.png, phocoena.png, stenella.png*

## Analysis :computer:
In this repository you will find raw data (.csv and data files) and code for analyses (code supplied as .R files)

📁 Data

As described above

📁 Code for analyses

*ossification_size_allometry_analyses.R, tympanic_bulla_shape_analyses.R, periotic_shape_analyses.R*

Before running analyses, save Data folder in the same directory as the R project. This will allow to import the data as detailed in the code provided.

## License 📃
This project is licensed under the MIT License - see the LICENSE.md file for details

## Session Info 📋
For reproducibility purposes, here is the output of devtools::session_info() used to perform the analyses in the publication.

```
R version 4.1.0 (2021-05-18)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19043)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] Morpho_2.8          ggplotify_0.0.5     ggpubr_0.4.0        borealis_2021.02.04 ggphylomorpho_0.1.0 geiger_2.0.7        ape_5.5            
 [8] gridExtra_2.3       png_0.1-7           rphylopic_0.3.0     ggfortify_0.4.11    gginnards_0.1.0     broom_0.7.8         geomorph_4.0.0     
[15] Matrix_1.3-3        rgl_0.106.8         RRPP_1.0.0          ggthemes_4.2.4      RColorBrewer_1.1-2  polynom_1.4-0       ggrepel_0.9.1      
[22] AICcmodavg_2.3-1    forcats_0.5.1       stringr_1.4.0       purrr_0.3.4         readr_1.4.0         tibble_3.1.2        ggplot2_3.3.4      
[29] tidyverse_1.3.1     dplyr_1.0.7         tidyr_1.1.3        

loaded via a namespace (and not attached):
  [1] readxl_1.3.1            backports_1.2.1         fastmatch_1.1-0         VGAM_1.1-5              plyr_1.8.6              igraph_1.2.6           
  [7] sp_1.4-5                splines_4.1.0           crosstalk_1.1.1         unmarked_1.1.1          gridBase_0.4-7          digest_0.6.27          
 [13] foreach_1.5.1           htmltools_0.5.1.1       fansi_0.5.0             magrittr_2.0.1          phytools_0.7-70         doParallel_1.0.16      
 [19] openxlsx_4.2.3          modelr_0.1.8            jpeg_0.1-8.1            colorspace_2.0-2        rvest_1.0.0             haven_2.4.1            
 [25] xfun_0.24               crayon_1.4.1            jsonlite_1.7.2          iterators_1.0.13        survival_3.2-11         phangorn_2.6.2         
 [31] glue_1.4.2              gtable_0.3.0            car_3.0-10              maps_3.3.0              abind_1.4-5             scales_1.1.1           
 [37] mvtnorm_1.1-1           DBI_1.1.1               rstatix_0.7.0           miniUI_0.1.1.1          Rcpp_1.0.6              plotrix_3.8-1          
 [43] xtable_1.8-4            tmvnsim_1.0-2           gridGraphics_0.5-1      foreign_0.8-81          subplex_1.6             deSolve_1.28           
 [49] stats4_4.1.0            htmlwidgets_1.5.3       httr_1.4.2              ellipsis_0.3.2          pkgconfig_2.0.3         dbplyr_2.1.1           
 [55] utf8_1.2.1              crul_1.1.0              tidyselect_1.1.1        rlang_0.4.11            manipulateWidget_0.11.0 later_1.2.0            
 [61] munsell_0.5.0           cellranger_1.1.0        tools_4.1.0             cli_2.5.0               generics_0.1.0          fastmap_1.1.0          
 [67] knitr_1.33              fs_1.5.0                zip_2.1.1               nlme_3.1-152            mime_0.11               xml2_1.3.2             
 [73] compiler_4.1.0          rstudioapi_0.13         curl_4.3.1              ggsignif_0.6.1          reprex_2.0.0            clusterGeneration_1.3.7
 [79] stringi_1.6.1           lattice_0.20-44         vctrs_0.3.8             pillar_1.6.1            lifecycle_1.0.0         BiocManager_1.30.12    
 [85] combinat_0.0-8          Rvcg_0.19.2             data.table_1.14.0       raster_3.4-13           colorRamps_2.3          httpuv_1.6.1           
 [91] R6_2.5.0                promises_1.2.0.1        rio_0.5.26              codetools_0.2-18        MASS_7.3-54             gtools_3.8.2           
 [97] assertthat_0.2.1        withr_2.4.2             httpcode_0.3.0          mnormt_2.0.2            expm_0.999-6            parallel_4.1.0         
[103] hms_1.1.0               quadprog_1.5-8          grid_4.1.0              coda_0.19-4             rvcheck_0.1.8           carData_3.0-4          
[109] numDeriv_2016.8-1.1     scatterplot3d_0.3-41    shiny_1.6.0             lubridate_1.7.10      
```
