
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

``` r
head(demo_algae_tax)
```

    ##   StationCode SampleDate Replicate CollectionMethodCode SampleTypeCode
    ## 1   000CAT148    8/10/10         1          Algae_SWAMP     Integrated
    ## 2   000CAT148    8/10/10         1          Algae_SWAMP     Integrated
    ## 3   000CAT148    8/10/10         1          Algae_SWAMP     Integrated
    ## 4   000CAT148    8/10/10         1          Algae_SWAMP     Integrated
    ## 5   000CAT148    8/10/10         1          Algae_SWAMP     Integrated
    ## 6   000CAT148    8/10/10         1          Algae_SWAMP     Integrated
    ##   BAResult Result                                FinalID
    ## 1        6     NA Achnanthidium exiguum var heterovalvum
    ## 2      225     NA             Achnanthidium minutissimum
    ## 3        2     NA                      Adlafia minuscula
    ## 4        8     NA                    Adlafia suchlandtii
    ## 5        1     NA                       Amphora copulata
    ## 6        1     NA                     Amphora inariensis

The output is in a wide format.

``` r
demo_results <- ASCI(demo_algae_tax)
demo_results
```

    ## # A tibble: 4 x 51
    ##   SampleID StationCode SampleDate Replicate SampleType D_ValveCount
    ##   <chr>    <chr>       <chr>          <dbl> <chr>             <int>
    ## 1 000CAT1~ 000CAT148   8/10/10            1 Integrate~          600
    ## 2 102PS01~ 102PS0139   8/9/10             1 Integrate~          600
    ## 3 105DLCD~ 105DLCDCC   5/19/09            1 Integrate~          600
    ## 4 105DLCD~ 105DLCDCC   6/23/09            1 Integrate~          600
    ## # ... with 45 more variables: S_EntityCount <int>, S_Biovolume <dbl>,
    ## #   D_NumberTaxa <dbl>, S_NumberTaxa <dbl>, H_NumberTaxa <dbl>,
    ## #   UnrecognizedTaxa <chr>, D_ASCI <dbl>, S_ASCI <dbl>, H_ASCI <dbl>,
    ## #   D_cnt.spp.BCG3 <dbl>, D_cnt.spp.BCG3_score <dbl>,
    ## #   D_pcnt.attributed.BCG3 <dbl>, D_pcnt.attributed.Cyclotella <dbl>,
    ## #   D_pcnt.attributed.OxyReg.DO_10 <dbl>,
    ## #   D_pcnt.attributed.Surirella <dbl>, D_prop.Cyclotella <dbl>,
    ## #   D_prop.Cyclotella_score <dbl>, D_prop.spp.OxyReq.DO_10 <dbl>,
    ## #   D_prop.spp.OxyReq.DO_10_score <dbl>, D_prop.Surirella <dbl>,
    ## #   D_prop.Surirella_score <dbl>, H_cnt.spp.IndicatorClass_Cu_high <dbl>,
    ## #   H_cnt.spp.IndicatorClass_Cu_high_score <dbl>, S_cnt.spp.BCG5 <dbl>,
    ## #   S_cnt.spp.BCG5_score <dbl>, S_cnt.spp.IndicatorClass_Cu_high <dbl>,
    ## #   S_cnt.spp.IndicatorClass_Cu_high_score <dbl>,
    ## #   S_cnt.spp.IndicatorClass_DOC_high <dbl>,
    ## #   S_cnt.spp.IndicatorClass_DOC_high_score <dbl>,
    ## #   S_cnt.spp.IndicatorClass_TP_high <dbl>,
    ## #   S_cnt.spp.IndicatorClass_TP_high_score <dbl>,
    ## #   S_pcnt.attributed.BCG5 <dbl>, S_pcnt.attributed.HiCu <dbl>,
    ## #   S_pcnt.attributed.HiDOC <dbl>, S_pcnt.attributed.HiTP.DO_10 <dbl>,
    ## #   H_cnt.ind.most.tol <dbl>, H_cnt.ind.most.tol_score <dbl>,
    ## #   H_pcnt.attributed.Cyclotella <dbl>, H_pcnt.attributed.HiCu <dbl>,
    ## #   H_pcnt.attributed.HiTolerance <dbl>, H_pcnt.attributed.NHHONF <dbl>,
    ## #   H_prop.Cyclotella <dbl>, H_prop.Cyclotella_score <dbl>,
    ## #   H_prop.spp.OrgN.NHHONF <dbl>, H_prop.spp.OrgN.NHHONF_score <dbl>
