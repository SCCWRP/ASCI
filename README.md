
## ASCI

#### *Marcus W. Beck (maintainer), <marcusb@sccwrp.org>, Susanna Theroux, <susannat@sccwrp.org>, Quynh-Thi Ho, <qthi.ho@gmail.com>, John Van Sickle*

[![Travis-CI Build
Status](https://travis-ci.org/SCCWRP/ASCI.svg?branch=master)](https://travis-ci.org/SCCWRP/ASCI)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/SCCWRP/ASCI?branch=master&svg=true)](https://ci.appveyor.com/project/SCCWRP/ASCI)
[![DOI](https://zenodo.org/badge/106055957.svg)](https://zenodo.org/badge/latestdoi/106055957)

R package materials to calculate the Algal Stream Condition Index (ASCI) using diatom, soft-bodied algae, or a hybrid appproach.

### Installation

Install the package as follows:

``` r
install.packages('devtools')
library(devtools)
install_github('SCCWRP/ASCI')
library(ASCI)
```

### Usage

The sample file `demo_algae_tax` is included to demonstrate the correct
format for input data. It is a `data.frame` of taxonomic data in long
format (one row per sample). See the help files for more information
(e.g., `?demo_algae_tax`)

``` r
demo_results <- ASCI(demo_algae_tax)
demo_results
```

    ## An object of class asci 
    ## Scores calculated for diatoms, sba, hybrid indices for 4 unique samples
    ## Use these functions for access: scores, Supp1_mmi

The output can be accessed using the defined methods.

``` r
.S4methods(class = 'asci')
```

    ## [1] scores    show      Supp1_mmi
    ## see '?methods' for accessing help and source code

``` r
show(demo_results)
```

    ## An object of class asci 
    ## Scores calculated for diatoms, sba, hybrid indices for 4 unique samples
    ## Use these functions for access: scores, Supp1_mmi

``` r
scores(demo_results)
```

    ## # A tibble: 12 x 4
    ##    taxa    SampleID              MMI MMI_Percentile
    ##    <chr>   <chr>               <dbl>          <dbl>
    ##  1 diatoms 000CAT148_8/10/10_1 1.27          0.620 
    ##  2 diatoms 102PS0139_8/9/10_1  0.920         0.0690
    ##  3 diatoms 105DLCDCC_5/19/09_1 1.31          0.696 
    ##  4 diatoms 105DLCDCC_6/23/09_1 1.34          0.747 
    ##  5 hybrid  000CAT148_8/10/10_1 1.19          0.836 
    ##  6 hybrid  102PS0139_8/9/10_1  0.934         0.0843
    ##  7 hybrid  105DLCDCC_5/19/09_1 1.09          0.507 
    ##  8 hybrid  105DLCDCC_6/23/09_1 1.13          0.648 
    ##  9 sba     000CAT148_8/10/10_1 1.10          0.845 
    ## 10 sba     102PS0139_8/9/10_1  0.966         0.298 
    ## 11 sba     105DLCDCC_5/19/09_1 0.915         0.129 
    ## 12 sba     105DLCDCC_6/23/09_1 1.07          0.740

``` r
Supp1_mmi(demo_results)
```

    ## # A tibble: 168 x 4
    ##    taxa    SampleID            Metric              Value
    ##    <chr>   <chr>               <chr>               <dbl>
    ##  1 diatoms 000CAT148_8/10/10_1 cnt.spp.BCG3       10    
    ##  2 diatoms 102PS0139_8/9/10_1  cnt.spp.BCG3        8    
    ##  3 diatoms 105DLCDCC_5/19/09_1 cnt.spp.BCG3        7    
    ##  4 diatoms 105DLCDCC_6/23/09_1 cnt.spp.BCG3       10    
    ##  5 diatoms 000CAT148_8/10/10_1 cnt.spp.BCG3_score  1.07 
    ##  6 diatoms 102PS0139_8/9/10_1  cnt.spp.BCG3_score  0.835
    ##  7 diatoms 105DLCDCC_5/19/09_1 cnt.spp.BCG3_score  0.716
    ##  8 diatoms 105DLCDCC_6/23/09_1 cnt.spp.BCG3_score  1.07 
    ##  9 diatoms 000CAT148_8/10/10_1 prop.Cyclotella     0    
    ## 10 diatoms 102PS0139_8/9/10_1  prop.Cyclotella     0    
    ## # ... with 158 more rows

<!-- Summary of ASCI peformance statewide: -->

<!-- ```{r} -->

<!-- perf(allscr) -->

<!-- ``` -->
