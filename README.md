# The healthy human airway microbial metagenome
<br>
Ajith Thavarasa,<sup>1*</sup> Marie-Madlen Pust,<sup>1*,#</sup> Ilona Rosenboom,<sup>1</sup> Lena Abeln,<sup>2</sup> Colin F. Davenport,<sup>2</sup> Lutz Wiehlmann,<sup>2</sup> Burkhard Tümmler<sup>1,3,#</sup> <br><br>
<sup>1</sup>Department for Pediatric Pneumology, Allergology and Neonatology, Hannover Medical School, Hannover, Germany<br>
<sup>2</sup>Research Core Unit Genomics, Hannover Medical School, Hannover, Germany<br>
<sup>3</sup>German Center for Lung Research, Biomedical Research in Endstage and Obstructive Lung Disease (BREATH), Hannover Medical School, Hannover, Germany<br>
<br><br>
<sup>*</sup>These authors contributed equally to this work
<br>
<sup>#</sup>Correspondence
<br><br>

# R environment 
<br>

```{r}
R version 4.1.1 (2021-08-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19044)

Matrix products: default

locale:
[1] LC_COLLATE=German_Germany.1252  LC_CTYPE=German_Germany.1252    LC_MONETARY=German_Germany.1252 LC_NUMERIC=C                   
[5] LC_TIME=German_Germany.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggVennDiagram_1.1.4 factoextra_1.0.7    ggplotify_0.1.0     pheatmap_1.0.12     compositions_2.0-4  forcats_0.5.1       tibble_3.1.3       
 [8] tidyverse_1.3.1     stringr_1.4.0       Boruta_7.0.0        randomForest_4.6-14 readxl_1.3.1        plyr_1.8.6          hrbrthemes_0.8.0   
[15] ggthemes_4.2.4      ggpubr_0.4.0        matrixStats_0.60.0  igraph_1.2.6        dplyr_1.0.7         Hmisc_4.5-0         Formula_1.2-4      
[22] survival_3.2-11     tidyr_1.1.3         viridis_0.6.1       viridisLite_0.4.0   ggrepel_0.9.1       ggplot2_3.3.5       vegan_2.5-7        
[29] lattice_0.20-45     permute_0.9-5       purrr_0.3.4         readr_2.0.1        

loaded via a namespace (and not attached):
  [1] backports_1.2.1     corrplot_0.90       systemfonts_1.0.2   splines_4.1.1       digest_0.6.27       yulab.utils_0.0.4   htmltools_0.5.1.1  
  [8] fansi_0.5.0         magrittr_2.0.1      checkmate_2.0.0     cluster_2.1.2       tzdb_0.1.2          openxlsx_4.2.4      modelr_0.1.8       
 [15] extrafont_0.17      bayesm_3.1-4        vroom_1.5.4         extrafontdb_1.0     jpeg_0.1-9          RVenn_1.1.0         colorspace_2.0-2   
 [22] rvest_1.0.2         textshaping_0.3.5   haven_2.4.3         xfun_0.25           crayon_1.4.1        jsonlite_1.7.2      glue_1.4.2         
 [29] gtable_0.3.0        car_3.0-11          Rttf2pt1_1.3.9      DEoptimR_1.0-9      abind_1.4-5         scales_1.1.1        DBI_1.1.1          
 [36] rstatix_0.7.0       Rcpp_1.0.7          htmlTable_2.2.1     units_0.7-2         gridGraphics_0.5-1  foreign_0.8-81      bit_4.0.4          
 [43] proxy_0.4-26        htmlwidgets_1.5.3   httr_1.4.2          RColorBrewer_1.1-2  ellipsis_0.3.2      pkgconfig_2.0.3     farver_2.1.0       
 [50] nnet_7.3-16         dbplyr_2.1.1        utf8_1.2.2          tidyselect_1.1.1    labeling_0.4.2      rlang_0.4.11        reshape2_1.4.4     
 [57] munsell_0.5.0       cellranger_1.1.0    tools_4.1.1         cli_3.0.1           generics_0.1.0      broom_0.7.9         evaluate_0.14      
 [64] yaml_2.2.1          ragg_1.1.3          knitr_1.33          bit64_4.0.5         fs_1.5.0            zip_2.2.0           robustbase_0.93-9  
 [71] nlme_3.1-152        xml2_1.3.2          compiler_4.1.1      rstudioapi_0.13     curl_4.3.2          png_0.1-7           e1071_1.7-8        
 [78] ggsignif_0.6.2      reprex_2.0.1        stringi_1.7.3       gdtools_0.2.4       Matrix_1.3-4        classInt_0.4-3      tensorA_0.36.2     
 [85] vctrs_0.3.8         pillar_1.6.2        lifecycle_1.0.0     data.table_1.14.0   cowplot_1.1.1       R6_2.5.1            latticeExtra_0.6-29
 [92] KernSmooth_2.23-20  gridExtra_2.3       rio_0.5.27          MASS_7.3-54         assertthat_0.2.1    withr_2.4.2         mgcv_1.8-36        
 [99] parallel_4.1.1      hms_1.1.0           grid_4.1.1          rpart_4.1-15        class_7.3-19        rmarkdown_2.10      carData_3.0-4      
[106] sf_1.0-2            lubridate_1.8.0     base64enc_0.1-3  
```
