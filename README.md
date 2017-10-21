
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

Two core functions are used to estimate scores.  The `oefun` function estimates O/E scores and the `pmmifun` funtion estimates pMMI scores. Both functions require taxonomic and site data as inputs. 



```r
results <- ASCI(demo_algae_tax, demo_algae_sitedata)
results
```

```
##       taxa            SampleID       MMI MMI_Percentile  O         E
## 1  diatoms 000CAT148_8/10/10_1 0.8078584     0.84992625  8  6.666486
## 2  diatoms  102PS0139_8/9/10_1 0.7690314     0.44229715  8  7.052908
## 3  diatoms 105DLCDCC_5/19/09_1 0.7877054     0.66384698  8  6.411555
## 4  diatoms 105DLCDCC_6/23/09_1 0.7306138     0.09443072  8  6.411555
## 5   hybrid 000CAT148_8/10/10_1 0.6909954     0.91676192  9  9.639016
## 6   hybrid  102PS0139_8/9/10_1 0.5397984     0.20206578 10 10.349072
## 7   hybrid 105DLCDCC_5/19/09_1 0.5543945     0.26757698  9  8.757375
## 8   hybrid 105DLCDCC_6/23/09_1 0.6014989     0.52822442 10  8.407181
## 9      sba 000CAT148_8/10/10_1 0.7329157     0.70136266  1  2.332670
## 10     sba  102PS0139_8/9/10_1 0.6364676     0.46639972  2  2.822729
## 11     sba 105DLCDCC_5/19/09_1 0.4347732     0.08604756  3  5.394556
## 12     sba 105DLCDCC_6/23/09_1 0.7948133     0.82160584  4  5.430622
##       OoverE OoverE_Percentile
## 1  1.2000325        0.44506519
## 2  1.1342840        0.08657412
## 3  1.2477472        0.77341730
## 4  1.2477472        0.77341730
## 5  0.9337053        0.20022150
## 6  0.9662702        0.28966901
## 7  1.0277053        0.49445596
## 8  1.1894593        0.92059418
## 9  0.4286934        0.10584953
## 10 0.7085342        0.75988848
## 11 0.5561162        0.35987830
## 12 0.7365638        0.81640445
```

Output from one index type (diatoms, soft-bodied algae, or hybrid) can be returned with the `tax` argument.

```r
ASCI(demo_algae_tax, demo_algae_sitedata, tax = 'diatoms')
```

```
##      taxa            SampleID       MMI MMI_Percentile O        E   OoverE
## 1 diatoms 000CAT148_8/10/10_1 0.8134297     0.84538563 8 6.666486 1.200033
## 2 diatoms  102PS0139_8/9/10_1 0.7809986     0.53232892 8 7.052908 1.134284
## 3 diatoms 105DLCDCC_5/19/09_1 0.7877054     0.60820055 8 6.411555 1.247747
## 4 diatoms 105DLCDCC_6/23/09_1 0.7306138     0.08493848 8 6.411555 1.247747
##   OoverE_Percentile
## 1        0.44506519
## 2        0.08657412
## 3        0.77341730
## 4        0.77341730
```

The output contains supplementary data as attributes.

```r
names(attributes(results))
```

```
## [1] "class"     "row.names" "names"     "Supp1_mmi" "Supp1_OE"  "Supp2_OE"
```

```r
attr(results, 'Supp1_mmi')
```

```
## # A tibble: 96 x 4
##       taxa            SampleID                  Metric      Value
##      <chr>               <chr>                   <chr>      <dbl>
##  1 diatoms 000CAT148_8/10/10_1       prop.ind.most.tol 0.04545455
##  2 diatoms  102PS0139_8/9/10_1       prop.ind.most.tol 0.15384615
##  3 diatoms 105DLCDCC_5/19/09_1       prop.ind.most.tol 0.06666667
##  4 diatoms 105DLCDCC_6/23/09_1       prop.ind.most.tol 0.10526316
##  5 diatoms 000CAT148_8/10/10_1 prop.ind.most.tol_score 0.89017919
##  6 diatoms  102PS0139_8/9/10_1 prop.ind.most.tol_score 0.66804689
##  7 diatoms 105DLCDCC_5/19/09_1 prop.ind.most.tol_score 0.84670814
##  8 diatoms 105DLCDCC_6/23/09_1 prop.ind.most.tol_score 0.76761043
##  9 diatoms 000CAT148_8/10/10_1           prop.spp.BCG3 0.45454545
## 10 diatoms  102PS0139_8/9/10_1           prop.spp.BCG3 0.46153846
## # ... with 86 more rows
```



