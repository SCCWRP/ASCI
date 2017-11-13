
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
results <- ASCI(demo_algae_tax, demo_algae_sitedata)
results
```

```
## An object of class asci 
## Scores calculated for diatoms, sba, hybrid 
## Use these functions for access: scores, Supp1_mmi, Supp1_OE, Supp2_OE, Supp3_OE
```

Output from one index type (diatoms, soft-bodied algae, or hybrid) can be returned with the `tax` argument.

```r
ASCI(demo_algae_tax, demo_algae_sitedata, tax = 'diatoms')
```

```
## An object of class asci 
## Scores calculated for diatoms 
## Use these functions for access: scores, Supp1_mmi, Supp1_OE, Supp2_OE, Supp3_OE
```

The output can be accessed using the defined methods.

```r
.S4methods(class = 'asci')
```

```
## [1] show      scores    Supp1_mmi Supp1_OE  Supp2_OE  Supp3_OE 
## see '?methods' for accessing help and source code
```

```r
scores(results)
```

```
## # A tibble: 12 x 8
##       taxa            SampleID       MMI MMI_Percentile     O         E
##  *   <chr>               <chr>     <dbl>          <dbl> <dbl>     <dbl>
##  1 diatoms 000CAT148_8/10/10_1 0.8134297     0.87518058     8  6.666486
##  2 diatoms  102PS0139_8/9/10_1 0.7534252     0.31270715     8  7.052908
##  3 diatoms 105DLCDCC_5/19/09_1 0.7877054     0.67306700     7  6.411555
##  4 diatoms 105DLCDCC_6/23/09_1 0.7306138     0.13319087     8  6.411555
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



