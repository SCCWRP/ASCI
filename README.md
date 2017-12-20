---
output:
  html_document:
    keep_md: yes
    toc: no
    self_contained: no
---

## ASCI

#### *Marcus W. Beck, marcusb@sccwrp.org, Susanna Theroux, susannat@sccwrp.org, John Van Sickle, VanSickle.John@epa.gov*

[![Travis-CI Build Status](https://travis-ci.org/SCCWRP/ASCI.svg?branch=master)](https://travis-ci.org/SCCWRP/ASCI)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/SCCWRP/ASCI?branch=master&svg=true)](https://ci.appveyor.com/project/SCCWRP/ASCI)
[![DOI](https://zenodo.org/badge/106055957.svg)](https://zenodo.org/badge/latestdoi/106055957)

R package materials to calculate the Algal Stream Condition Index (ASCI) based on O/E and pMMI scores using diatom, soft-bodied algae, or a hybrid appproach.

### Installation

Install the package as follows:


```r
install.packages('devtools')
library(devtools)
install_github('SCCWRP/ASCI')
library(ASCI)
```

### Citation

Please cite the current release as follows:

Marcus W Beck, Susanna Theroux, John Van Sickle (2017, October 11). SCCWRP/ASCI: v0.0.1. Zenodo. http://doi.org/10.5281/zenodo.1008839

### Usage

Two sample files are included to demonstrate the correct format for input data. The `demo_algae_tax` file is a `data.frame` of taxonomic data in long format (one row per sample).  The `demo_algae_sitedata` file is a `data.frame` of site data.  Site names must match between the files. See the help files for more information (e.g., `?demo_algae_tax`)

The core function is `ASCI` that estimates O/E and pMMI scores.  This provides a wrapper to the separate `oefun` and `pmmifun` functions. All functions require taxonomic and site data as inputs. 


```r
demo_results <- ASCI(demo_algae_tax, demo_algae_sitedata)
demo_results
```

```
## An object of class asci 
## Scores calculated for diatoms, sba, hybrid indices for 4 unique samples
## Use these functions for access: perf, scores, Supp1_mmi, Supp1_OE, Supp2_OE
```

The output can be accessed using the defined methods.

```r
.S4methods(class = 'asci')
```

```
## [1] perf      scores    show      Supp1_mmi Supp1_OE  Supp2_OE 
## see '?methods' for accessing help and source code
```

```r
scores(demo_results)
```

```
## # A tibble: 12 x 8
##       taxa            SampleID       MMI MMI_Percentile     O         E
##  *   <chr>               <chr>     <dbl>          <dbl> <dbl>     <dbl>
##  1 diatoms 000CAT148_8/10/10_1 0.8078584     0.86431742     8  6.666486
##  2 diatoms  102PS0139_8/9/10_1 0.7534252     0.31653108     8  7.052908
##  3 diatoms 105DLCDCC_5/19/09_1 0.7877054     0.69705073     8  6.411555
##  4 diatoms 105DLCDCC_6/23/09_1 0.7306138     0.12746846     8  6.411555
##  5  hybrid 000CAT148_8/10/10_1 0.6909954     0.91676192     9  9.639016
##  6  hybrid  102PS0139_8/9/10_1 0.5397984     0.20206578    10 10.349072
##  7  hybrid 105DLCDCC_5/19/09_1 0.5543945     0.26757698     9  8.757375
##  8  hybrid 105DLCDCC_6/23/09_1 0.6014989     0.52822442    10  8.407181
##  9     sba 000CAT148_8/10/10_1 0.7329157     0.70136266     1  2.332670
## 10     sba  102PS0139_8/9/10_1 0.6364676     0.46639972     2  2.822729
## 11     sba 105DLCDCC_5/19/09_1 0.4347732     0.08604756     3  5.394556
## 12     sba 105DLCDCC_6/23/09_1 0.7948133     0.82160584     4  5.430622
## # ... with 2 more variables: OoverE <dbl>, OoverE_Percentile <dbl>
```

Summary of ASCI peformance statewide:

```r
perf(allscr)
```

```
## # A tibble: 18 x 9
##        grp   ind        typ   cls       ave        fst    prc_amg
##      <chr> <chr>      <chr> <chr>     <dbl>      <dbl>      <dbl>
##  1 diatoms    oe Predictive    rc 1.0331498  0.9040594 0.15001039
##  2 diatoms    oe Predictive    rv 0.9610270  0.7859889 0.22526853
##  3  hybrid    oe Predictive    rc 1.0396242  1.5954080 0.15977061
##  4  hybrid    oe Predictive    rv 1.0303972  2.2887879 0.18613214
##  5     sba    oe Predictive    rc 1.0055350  2.2819198 0.35127413
##  6     sba    oe Predictive    rv 0.7367404  1.6572882 0.36675319
##  7 diatoms   mmi Predictive    rc 0.7376118 30.2288213 0.09417096
##  8 diatoms   mmi Predictive    rv 0.7487854  9.2282463 0.07980960
##  9  hybrid   mmi Predictive    rc 0.5733118  2.5709406 0.09956694
## 10  hybrid   mmi Predictive    rv 0.5869762  1.3315131 0.09867553
## 11     sba   mmi Predictive    rc 0.6732414  1.6250796 0.12878218
## 12     sba   mmi Predictive    rv 0.6665803  1.0523804 0.13849540
## 13 diatoms    oe       Null    rc 0.9828557 11.4249511 0.19914499
## 14 diatoms    oe       Null    rv 0.9630963  2.4226541 0.23680292
## 15  hybrid    oe       Null    rc 1.0613450  1.2757346 0.18416892
## 16  hybrid    oe       Null    rv 1.0279032  1.0932668 0.20920504
## 17     sba    oe       Null    rc 0.8893550  5.7128060 0.38847344
## 18     sba    oe       Null    rv 0.7677371  1.7063871 0.39954320
## # ... with 2 more variables: prc_wth <dbl>, res_tst <dbl>
```




