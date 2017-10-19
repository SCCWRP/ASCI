
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
# O/E scores
oescr <- oefun(demo_algae_tax, demo_algae_sitedata)

# pMMI scores
pmmiscr <- pmmifun(demo_algae_tax, demo_algae_sitedata)
```

The output from each function is a `list` with the relevant score information.

```r
names(oescr)
```

```
## [1] "Scores.diatoms" "Scores.sba"     "Scores.hybrid"
```

```r
names(pmmiscr)
```

```
## [1] "taxa"     "met"      "SampleID" "val"
```

Scores for each index can be viewed as follows.

```r
lapply(oescr, function(x) head(x$OE.scores))
```

```
## $Scores.diatoms
##              SampleID O        E   OoverE OoverE_Percentile
## 1 000CAT148_8/10/10_1 8 6.666486 1.200033        0.44506519
## 2  102PS0139_8/9/10_1 8 7.052908 1.134284        0.08657412
## 3 105DLCDCC_5/19/09_1 8 6.411555 1.247747        0.77341730
## 4 105DLCDCC_6/23/09_1 8 6.411555 1.247747        0.77341730
## 
## $Scores.sba
##              SampleID O        E    OoverE OoverE_Percentile
## 1 000CAT148_8/10/10_1 1 2.332670 0.4286934         0.1058495
## 2  102PS0139_8/9/10_1 2 2.822729 0.7085342         0.7598885
## 3 105DLCDCC_5/19/09_1 3 5.394556 0.5561162         0.3598783
## 4 105DLCDCC_6/23/09_1 4 5.430622 0.7365638         0.8164044
## 
## $Scores.hybrid
##              SampleID  O         E    OoverE OoverE_Percentile
## 1 000CAT148_8/10/10_1  9  9.639016 0.9337053         0.2002215
## 2  102PS0139_8/9/10_1 10 10.349072 0.9662702         0.2896690
## 3 105DLCDCC_5/19/09_1  9  8.757375 1.0277053         0.4944560
## 4 105DLCDCC_6/23/09_1 10  8.407181 1.1894593         0.9205942
```

```r
head(pmmiscr)
```

```
## # A tibble: 6 x 4
##      taxa                   met            SampleID        val
##     <chr>                 <chr>               <chr>      <dbl>
## 1 diatoms  prop.spp.Salinity.BF 000CAT148_8/10/10_1 0.07444986
## 2 diatoms  prop.spp.Salinity.BF  102PS0139_8/9/10_1 0.01584088
## 3 diatoms  prop.spp.Salinity.BF 105DLCDCC_5/19/09_1 0.01231552
## 4 diatoms  prop.spp.Salinity.BF 105DLCDCC_6/23/09_1 0.09878168
## 5 diatoms prop.spp.HighMotility 000CAT148_8/10/10_1 0.17391304
## 6 diatoms prop.spp.HighMotility  102PS0139_8/9/10_1 0.20000000
```


