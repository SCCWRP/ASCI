
## ASCI

#### *Marcus W. Beck  <marcusb@sccwrp.org>, Robert Butler (maintainer) <robertb@sccwrp.org>, Susanna Theroux, <susannat@sccwrp.org>, Quynh-Thi Ho, <qthi.ho@gmail.com>

[![Travis-CI Build
Status](https://travis-ci.org/SCCWRP/ASCI.svg?branch=master)](https://travis-ci.org/SCCWRP/ASCI)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/SCCWRP/ASCI?branch=master&svg=true)](https://ci.appveyor.com/project/SCCWRP/ASCI)
[![DOI](https://zenodo.org/badge/106055957.svg)](https://zenodo.org/badge/latestdoi/106055957)

R package materials to calculate the Algal Stream Condition Index (ASCI)
based on pMMI scores using diatom, soft-bodied algae, or a
hybrid appproach. A link to the ASCI manuscript can be found here: https://www.sciencedirect.com/science/article/pii/S1470160X20303587

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

    ## # A tibble: 3 x 27
    ##   StationCode CondQR50 SITE_ELEV TEMP_00_09 KFCT_AVE  AtmCa PPT_00_09 MAX_ELEV
    ##   <chr>          <int>     <int>      <int>    <dbl>  <dbl>     <dbl>    <int>
    ## 1 404M07357         NA       199       2456    0.278 0.0554    55570.      783
    ## 2 801M16916         NA       197       2685    0.185 0.0652    25406.     3480
    ## 3 909M24937         NA       582       2442    0.202 0.106     37972.     1980
    ## # … with 19 more variables: CaO_Mean <dbl>, MgO_Mean <dbl>, S_Mean <dbl>,
    ## #   UCS_Mean <dbl>, LPREM_mean <dbl>, AtmMg <dbl>, AtmSO4 <dbl>, MINP_WS <dbl>,
    ## #   MEANP_WS <dbl>, SumAve_P <dbl>, TMAX_WS <dbl>, XWD_WS <dbl>,
    ## #   MAXWD_WS <dbl>, LST32AVE <dbl>, BDH_AVE <dbl>, PRMH_AVE <dbl>, PSA6C <chr>,
    ## #   XerMtn <lgl>, AREA_SQKM <dbl>

The output is in a wide format.

``` r
demo_results <- ASCI(demo_algae_tax, demo_station)
demo_results
```

    ## # A tibble: 3 x 84
    ##   SampleID StationCode SampleDate          Replicate SampleType D_ValveCount
    ##   <chr>    <chr>       <dttm>                  <dbl> <chr>             <int>
    ## 1 404M073… 404M07357   2016-06-13 00:00:00         1 Integrate…          600
    ## 2 801M169… 801M16916   2016-05-25 00:00:00         1 Microalga…          600
    ## 3 909M249… 909M24937   2016-06-22 00:00:00         1 Macroalga…          600
    ## # … with 78 more variables: S_EntityCount <int>, S_Biovolume <dbl>,
    ## #   D_NumberTaxa <int>, S_NumberTaxa <int>, H_NumberTaxa <int>,
    ## #   UnrecognizedTaxa <chr>, D_ASCI <dbl>, S_ASCI <dbl>, H_ASCI <dbl>,
    ## #   D_cnt.spp.most.tol_pct_att <dbl>, D_cnt.spp.most.tol_pred <dbl>,
    ## #   D_cnt.spp.most.tol_raw <int>, D_cnt.spp.most.tol_scr <dbl>,
    ## #   D_EpiRho.richness_pct_att <dbl>, D_EpiRho.richness_pred <dbl>,
    ## #   D_EpiRho.richness_raw <int>, D_EpiRho.richness_scr <dbl>,
    ## #   D_prop.spp.IndicatorClass_TN_low_pct_att <dbl>,
    ## #   D_prop.spp.IndicatorClass_TN_low_pred <dbl>,
    ## #   D_prop.spp.IndicatorClass_TN_low_raw <dbl>,
    ## #   D_prop.spp.IndicatorClass_TN_low_scr <dbl>,
    ## #   D_prop.spp.Planktonic_pct_att <dbl>, D_prop.spp.Planktonic_pred <dbl>,
    ## #   D_prop.spp.Planktonic_raw <dbl>, D_prop.spp.Planktonic_scr <dbl>,
    ## #   D_prop.spp.Trophic.E_pct_att <dbl>, D_prop.spp.Trophic.E_pred <dbl>,
    ## #   D_prop.spp.Trophic.E_raw <dbl>, D_prop.spp.Trophic.E_scr <dbl>,
    ## #   D_Salinity.BF.richness_pct_att <dbl>, D_Salinity.BF.richness_pred <dbl>,
    ## #   D_Salinity.BF.richness_raw <int>, D_Salinity.BF.richness_scr <dbl>,
    ## #   S_prop.spp.IndicatorClass_DOC_high_pct_att <dbl>,
    ## #   S_prop.spp.IndicatorClass_DOC_high_raw <dbl>,
    ## #   S_prop.spp.IndicatorClass_DOC_high_scr <dbl>,
    ## #   S_prop.spp.IndicatorClass_NonRef_pct_att <dbl>,
    ## #   S_prop.spp.IndicatorClass_NonRef_raw <dbl>,
    ## #   S_prop.spp.IndicatorClass_NonRef_scr <dbl>,
    ## #   S_prop.spp.IndicatorClass_TP_high_pct_att <dbl>,
    ## #   S_prop.spp.IndicatorClass_TP_high_raw <dbl>,
    ## #   S_prop.spp.IndicatorClass_TP_high_scr <dbl>, S_prop.spp.ZHR_pct_att <dbl>,
    ## #   S_prop.spp.ZHR_raw <dbl>, S_prop.spp.ZHR_scr <dbl>,
    ## #   H_cnt.spp.IndicatorClass_TP_high_pct_att <dbl>,
    ## #   H_cnt.spp.IndicatorClass_TP_high_pred <dbl>,
    ## #   H_cnt.spp.IndicatorClass_TP_high_raw <int>,
    ## #   H_cnt.spp.IndicatorClass_TP_high_scr <dbl>,
    ## #   H_cnt.spp.most.tol_pct_att <dbl>, H_cnt.spp.most.tol_pred <dbl>,
    ## #   H_cnt.spp.most.tol_raw <int>, H_cnt.spp.most.tol_scr <dbl>,
    ## #   H_EpiRho.richness_pct_att <dbl>, H_EpiRho.richness_pred <dbl>,
    ## #   H_EpiRho.richness_raw <int>, H_EpiRho.richness_scr <dbl>,
    ## #   H_OxyRed.DO_30.richness_pct_att <dbl>, H_OxyRed.DO_30.richness_pred <dbl>,
    ## #   H_OxyRed.DO_30.richness_raw <int>, H_OxyRed.DO_30.richness_scr <dbl>,
    ## #   H_prop.spp.Planktonic_pct_att <dbl>, H_prop.spp.Planktonic_pred <dbl>,
    ## #   H_prop.spp.Planktonic_raw <dbl>, H_prop.spp.Planktonic_scr <dbl>,
    ## #   H_prop.spp.Trophic.E_pct_att <dbl>, H_prop.spp.Trophic.E_pred <dbl>,
    ## #   H_prop.spp.Trophic.E_raw <dbl>, H_prop.spp.Trophic.E_scr <dbl>,
    ## #   H_prop.spp.ZHR_pct_att <dbl>, H_prop.spp.ZHR_raw <dbl>,
    ## #   H_prop.spp.ZHR_scr <dbl>, H_Salinity.BF.richness_pct_att <dbl>,
    ## #   H_Salinity.BF.richness_pred <dbl>, H_Salinity.BF.richness_raw <int>,
    ## #   H_Salinity.BF.richness_scr <dbl>, Comments <chr>, version_number <chr>


### FAQ

#### Missing data 
If a single algal assemblage is submitted (e.g. no soft algae taxa submitted), then the corresponding metrics and indices to the missing assemblage will return NA. However, if an assemblage is submitted but no taxa are attributed for the corresponding metrics, then the metrics will score the worst possible score. Samples with single assemblages submitted will result in a warning message. With few exceptions, missing values in stations data are not allowed. 

#### Bad or Missing Field Names
All required field names must be present in input files. Please be sure to match the field names provided above. Although we have implemented scripts to make the inputs case-insensitive, we recommend conforming to the capitalizations shown above.

#### Stations with Catchments that Include Parts in Mexico
Portions of some streams include areas in Mexico. Because the geodatabases used to calculate ASCI predictors do not currently include this area, the ASCI cannot be calculated properly for these sites. The geodatabases will be updated within the next few months. In the interim, we make the following recommendations: If more than 90% of the area of a watershed is within California, treat the state boundary as the edge of the watershed and calculate the predictors accordingly. However, you should interpret these results with caution, particularly if the portion within Mexico contains substantially different natural features. For watersheds that are less than 90% within California, we recommend using the Southern California Algal Indices of Biotic Integrity (Fetscher et al. 2014) as a substitute index.

#### Unrecognized Taxa
Novel or misspelled species names will not be recognized by the calculator and will be output as unrecognized taxa. Users should modify these species in agreement with the SWAMP species lists and re-run the calculator. The calculator STE currently reflects the 2019 SWAMP lookup lists and is available to view [HERE](https://github.com/SCCWRP/ASCI/blob/master/data/STE.RData).

### Metadata
Resources: <a href="https://github.com/SCCWRP/ASCIsop">SOP</a><br>
Contact: <a href="https://www.sccwrp.org/about/staff/susanna-theroux/">Susanna Theroux</a><br>


