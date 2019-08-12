
## ASCI

#### *Marcus W. Beck (maintainer), <marcusb@sccwrp.org>, Susanna Theroux, <susannat@sccwrp.org>, Quynh-Thi Ho, <qthi.ho@gmail.com>, John Van Sickle*

[![Travis-CI Build
Status](https://travis-ci.org/SCCWRP/ASCI.svg?branch=master)](https://travis-ci.org/SCCWRP/ASCI)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/SCCWRP/ASCI?branch=master&svg=true)](https://ci.appveyor.com/project/SCCWRP/ASCI)
[![DOI](https://zenodo.org/badge/106055957.svg)](https://zenodo.org/badge/latestdoi/106055957)

R package materials to calculate the Algal Stream Condition Index (ASCI)
based on O/E and pMMI scores using diatom, soft-bodied algae, or a
hybrid appproach.

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

The output is in a wide format.

``` r
demo_results <- ASCI(demo_algae_tax)
demo_results
```

    ## # A tibble: 4 x 31
    ##   SampleID D_cnt.spp.BCG3 D_cnt.spp.BCG3_~ D_MMI D_MMI_Percentile
    ##   <chr>             <dbl>            <dbl> <dbl>            <dbl>
    ## 1 000CAT1~             10            1.07  1.09             0.798
    ## 2 102PS01~              8            0.835 1.03             0.309
    ## 3 105DLCD~              7            0.716 0.999            0.122
    ## 4 105DLCD~             10            1.07  1.09             0.798
    ## # ... with 26 more variables: D_prop.Cyclotella <dbl>,
    ## #   D_prop.Cyclotella_score <dbl>, D_prop.spp.OxyReq.DO_10 <dbl>,
    ## #   D_prop.spp.OxyReq.DO_10_score <dbl>, D_prop.Surirella <dbl>,
    ## #   D_prop.Surirella_score <dbl>, H_cnt.ind.most.tol <dbl>,
    ## #   H_cnt.ind.most.tol_score <dbl>,
    ## #   H_cnt.spp.IndicatorClass_Cu_high <dbl>,
    ## #   H_cnt.spp.IndicatorClass_Cu_high_score <dbl>, H_MMI <dbl>,
    ## #   H_MMI_Percentile <dbl>, H_prop.Cyclotella <dbl>,
    ## #   H_prop.Cyclotella_score <dbl>, H_prop.spp.OrgN.NHHONF <dbl>,
    ## #   H_prop.spp.OrgN.NHHONF_score <dbl>, S_cnt.spp.BCG5 <dbl>,
    ## #   S_cnt.spp.BCG5_score <dbl>, S_cnt.spp.IndicatorClass_Cu_high <dbl>,
    ## #   S_cnt.spp.IndicatorClass_Cu_high_score <dbl>,
    ## #   S_cnt.spp.IndicatorClass_DOC_high <dbl>,
    ## #   S_cnt.spp.IndicatorClass_DOC_high_score <dbl>,
    ## #   S_cnt.spp.IndicatorClass_TP_high <dbl>,
    ## #   S_cnt.spp.IndicatorClass_TP_high_score <dbl>, S_MMI <dbl>,
    ## #   S_MMI_Percentile <dbl>
