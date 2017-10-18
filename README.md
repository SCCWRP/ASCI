
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
## [1] "d.results.scored"      "sba.results.scored"    "hybrid.results.scored"
```

Scores for each index can be viewed as follows.

```r
lapply(oescr, function(x) head(x$OE.scores))
```

```
## $Scores.diatoms
##                     O        E   OoverE
## 000CAT148_8/10/10_1 8 6.666486 1.200033
## 102PS0139_8/9/10_1  8 7.052908 1.134284
## 105DLCDCC_5/19/09_1 7 6.411555 1.091779
## 105DLCDCC_6/23/09_1 8 6.411555 1.247747
## 
## $Scores.sba
##                     O        E    OoverE
## 000CAT148_8/10/10_1 1 2.332670 0.4286934
## 102PS0139_8/9/10_1  2 2.822729 0.7085342
## 105DLCDCC_5/19/09_1 3 5.394556 0.5561162
## 105DLCDCC_6/23/09_1 4 5.430622 0.7365638
## 
## $Scores.hybrid
##                      O         E    OoverE
## 000CAT148_8/10/10_1  9  9.639016 0.9337053
## 102PS0139_8/9/10_1  10 10.349072 0.9662702
## 105DLCDCC_5/19/09_1  9  8.757375 1.0277053
## 105DLCDCC_6/23/09_1 10  8.407181 1.1894593
```

```r
lapply(pmmiscr, head)
```

```
## $d.results.scored
##   prop.spp.Salinity.BF prop.spp.HighMotility prop.ind.most.tol
## 1            0.7712797             0.5882098         0.8942293
## 2            0.7745123             0.4962933         0.6905672
## 3            0.7746800             0.4587652         0.8256892
## 4            0.7698189             0.3850259         0.7676104
##   prop.spp.BCG3 diatom.pMMI
## 1             1   0.8134297
## 2             1   0.7403432
## 3             1   0.7647836
## 4             1   0.7306138
## 
## $sba.results.scored
##   cnt.spp.IndicatorClass_TP_high prop.spp.IndicatorClass_DOC_high
## 1                      0.7080179                        0.3161360
## 2                      0.7080179                        0.7287864
## 3                      0.7080179                        0.1235659
## 4                      0.7080179                        0.5637263
##   prop.spp.Green prop.spp.BCG3  sba.pMMI
## 1      0.9075090     1.0000000 0.7329157
## 2      0.7223031     0.3867628 0.6364676
## 3      0.9075090     0.0000000 0.4347732
## 4      0.9075090     1.0000000 0.7948133
## 
## $hybrid.results.scored
##   prop.spp.IndicatorClass_DOC_high prop.spp.Trophic.I prop.spp.ZHR
## 1                        0.6000610          1.0000000    0.1639208
## 2                        0.5763431          0.5828506    0.0000000
## 3                        0.4360691          0.9159080    0.0000000
## 4                        0.5283145          0.6010646    0.2766163
##   prop.spp.BCG3 hybri.pMMI
## 1     1.0000000  0.6909954
## 2     1.0000000  0.5397984
## 3     0.8656009  0.5543945
## 4     1.0000000  0.6014989
```


