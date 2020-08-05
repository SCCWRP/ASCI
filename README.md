
## ASCI

#### *Marcus W. Beck, <marcusb@sccwrp.org>; Susanna Theroux, <susannat@sccwrp.org>; Robert Butler (maintainer), <robertb@sccwrp.org>; Quynh-Thi Ho, <qthi.ho@gmail.com>; Duy Nguyen, <duynguyen1993@csu.fullerton.edu> *

[![Travis-CI Build
Status](https://travis-ci.org/SCCWRP/ASCI.svg?branch=master)](https://travis-ci.org/SCCWRP/ASCI)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/SCCWRP/ASCI?branch=master&svg=true)](https://ci.appveyor.com/project/SCCWRP/ASCI)
[![DOI](https://zenodo.org/badge/106055957.svg)](https://zenodo.org/badge/latestdoi/106055957)

R package materials to calculate the Algal Stream Condition Index (ASCI)
based on predictive MMI scores using diatom, soft-bodied algae, or a
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

The sample files `demo_algae_tax` and `demo_station` are included to
demonstrate the correct formats for the input data. The `demo_algae_tax`
file is a `data.frame` of taxonomic data in long format (one row per
sample). The `demo_station` file is a `data.frame` of GIS predictors in
wide format, one row per station. See the help files for more
information (e.g., `?demo_algae_tax`). Also see the help file for
`chkinp()` and `calcgis()` for requirements of each file to work with
the ASCI.

``` r
head(demo_algae_tax)
```

    ## # A tibble: 6 x 7
    ##   StationCode SampleDate          Replicate SampleTypeCode BAResult  Result
    ##   <chr>       <dttm>                  <dbl> <chr>             <dbl>   <dbl>
    ## 1 909M24937   2016-06-22 00:00:00         1 Macroalgae           NA  1.21e9
    ## 2 909M24937   2016-06-22 00:00:00         1 Epiphyte             15 NA     
    ## 3 909M24937   2016-06-22 00:00:00         1 Epiphyte             83 NA     
    ## 4 909M24937   2016-06-22 00:00:00         1 Epiphyte              2 NA     
    ## 5 909M24937   2016-06-22 00:00:00         1 Integrated            5 NA     
    ## 6 909M24937   2016-06-22 00:00:00         1 Integrated           53 NA     
    ## # … with 1 more variable: FinalID <chr>

``` r
head(demo_station)
```

    ## # A tibble: 3 x 26
    ##   StationCode CondQR50 SITE_ELEV TEMP_00_09 KFCT_AVE  AtmCa PPT_00_09
    ##   <chr>          <int>     <int>      <int>    <dbl>  <dbl>     <dbl>
    ## 1 404M07357         NA       199       2456    0.278 0.0554    55570.
    ## 2 801M16916         NA       197       2685    0.185 0.0652    25406.
    ## 3 909M24937         NA       582       2442    0.202 0.106     37972.
    ## # … with 19 more variables: MAX_ELEV <int>, CaO_Mean <dbl>,
    ## #   MgO_Mean <dbl>, S_Mean <dbl>, UCS_Mean <dbl>, LPREM_mean <dbl>,
    ## #   AtmMg <dbl>, AtmSO4 <dbl>, MINP_WS <dbl>, MEANP_WS <dbl>,
    ## #   SumAve_P <dbl>, TMAX_WS <dbl>, XWD_WS <dbl>, MAXWD_WS <dbl>,
    ## #   LST32AVE <dbl>, BDH_AVE <dbl>, PRMH_AVE <dbl>, PSA6C <chr>,
    ## #   XerMtn <lgl>

The output is in a wide format.

``` r
demo_results <- ASCI(demo_algae_tax, demo_station)
demo_results
```

    ## # A tibble: 3 x 55
    ##   SampleID StationCode SampleDate          Replicate SampleType
    ##   <chr>    <chr>       <dttm>                  <int> <chr>     
    ## 1 404M073… 404M07357   2016-06-13 00:00:00         1 Integrate…
    ## 2 801M169… 801M16916   2016-05-25 00:00:00         1 Microalga…
    ## 3 909M249… 909M24937   2016-06-22 00:00:00         1 Macroalga…
    ## # … with 50 more variables: D_ValveCount <int>, S_EntityCount <int>,
    ## #   S_Biovolume <dbl>, D_NumberTaxa <dbl>, S_NumberTaxa <dbl>,
    ## #   H_NumberTaxa <dbl>, UnrecognizedTaxa <chr>, D_ASCI <dbl>,
    ## #   S_ASCI <dbl>, H_ASCI <dbl>, D_pct_att_prp_spp_BCG12 <dbl>,
    ## #   D_pct_att_prp_spp_OxRq_DO100_75 <dbl>,
    ## #   D_pct_att_prp_spp_Salinity_BF <dbl>,
    ## #   D_pct_att_prp_spp_Trophic_E <dbl>, D_prp_spp_BCG12_mod <dbl>,
    ## #   D_prp_spp_BCG12_mod_scr <dbl>, D_prp_spp_BCG12_pred <dbl>,
    ## #   D_prp_spp_OxRq_DO100_75_raw <dbl>,
    ## #   D_prp_spp_OxRq_DO100_75_raw_scr <dbl>,
    ## #   D_prp_spp_Salinity_BF_mod <dbl>, D_prp_spp_Salinity_BF_mod_scr <dbl>,
    ## #   D_prp_spp_Salinity_BF_pred <dbl>, D_prp_spp_Trophic_E_mod <dbl>,
    ## #   D_prp_spp_Trophic_E_mod_scr <dbl>, D_prp_spp_Trophic_E_pred <dbl>,
    ## #   H_OxRd_DO_30_richness_mod <dbl>, H_OxRd_DO_30_richness_mod_scr <dbl>,
    ## #   H_OxRd_DO_30_richness_pred <dbl>, H_pct_att_OxRd_DO_30_richness <dbl>,
    ## #   H_pct_att_prp_spp_BCG4 <dbl>, H_pct_att_prp_spp_IC_DOC_high <dbl>,
    ## #   H_pct_att_Salinity_BF_richness <dbl>, H_prp_spp_BCG4_mod <dbl>,
    ## #   H_prp_spp_BCG4_mod_scr <dbl>, H_prp_spp_BCG4_pred <dbl>,
    ## #   H_prp_spp_IC_DOC_high_raw <dbl>, H_prp_spp_IC_DOC_high_raw_scr <dbl>,
    ## #   H_Salinity_BF_richness_mod <dbl>,
    ## #   H_Salinity_BF_richness_mod_scr <dbl>,
    ## #   H_Salinity_BF_richness_pred <dbl>, S_cnt_spp_IC_DOC_high_raw <dbl>,
    ## #   S_cnt_spp_IC_DOC_high_raw_scr <dbl>,
    ## #   S_pct_att_cnt_spp_IC_DOC_high <dbl>, S_pct_att_prp_spp_BCG45 <dbl>,
    ## #   S_pct_att_prp_spp_Green <dbl>, S_prp_spp_BCG45_raw <dbl>,
    ## #   S_prp_spp_BCG45_raw_scr <dbl>, S_prp_spp_Green_raw <dbl>,
    ## #   S_prp_spp_Green_raw_scr <dbl>, Comments <chr>
