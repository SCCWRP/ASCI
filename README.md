
## ASCI

#### *Marcus W. Beck (maintainer), <marcusb@sccwrp.org>, Susanna Theroux, <susannat@sccwrp.org>*

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

### Citation

Please cite the current release as follows:

Marcus W Beck, Susanna Theroux, John Van Sickle (2017, October 11).
SCCWRP/ASCI: v0.0.1. Zenodo. <http://doi.org/10.5281/zenodo.1008839>

### Usage

Two sample files are included to demonstrate the correct format for
input data. The `demo_algae_tax` file is a `data.frame` of taxonomic
data in long format (one row per sample). The `demo_algae_sitedata` file
is a `data.frame` of site data. Site names must match between the files.
See the help files for more information (e.g., `?demo_algae_tax`)

The core function is `ASCI` that estimates O/E and pMMI scores. This
provides a wrapper to the separate `oefun` and `pmmifun` functions. All
functions require taxonomic and site data as inputs.

``` r
demo_results <- ASCI(demo_algae_tax, demo_algae_sitedata)
demo_results
```

    ## An object of class asci 
    ## Scores calculated for diatoms, sba, hybrid indices for 4 unique samples
    ## Use these functions for access: perf, scores, Supp1_mmi, Supp1_OE, Supp2_OE

The output can be accessed using the defined methods.

``` r
.S4methods(class = 'asci')
```

    ## [1] perf      scores    show      Supp1_mmi Supp1_OE  Supp2_OE 
    ## see '?methods' for accessing help and source code

``` r
scores(demo_results)
```

    ## # A tibble: 12 x 8
    ##    taxa   SampleID   MMI MMI_Percentile     O     E OoverE OoverE_Percenti~
    ##    <chr>  <chr>    <dbl>          <dbl> <dbl> <dbl>  <dbl>            <dbl>
    ##  1 diato~ 000CAT1~ 0.808         0.828      8  6.67  1.20            0.445 
    ##  2 diato~ 102PS01~ 0.781         0.551      8  7.05  1.13            0.0866
    ##  3 diato~ 105DLCD~ 0.788         0.630      8  6.41  1.25            0.773 
    ##  4 diato~ 105DLCD~ 0.731         0.0798     8  6.41  1.25            0.773 
    ##  5 hybrid 000CAT1~ 0.691         0.917      9  9.64  0.934           0.200 
    ##  6 hybrid 102PS01~ 0.540         0.202     10 10.3   0.966           0.290 
    ##  7 hybrid 105DLCD~ 0.554         0.268      9  8.76  1.03            0.494 
    ##  8 hybrid 105DLCD~ 0.601         0.528     10  8.41  1.19            0.921 
    ##  9 sba    000CAT1~ 0.733         0.701      1  2.33  0.429           0.106 
    ## 10 sba    102PS01~ 0.636         0.466      2  2.82  0.709           0.760 
    ## 11 sba    105DLCD~ 0.435         0.0860     3  5.39  0.556           0.360 
    ## 12 sba    105DLCD~ 0.795         0.822      4  5.43  0.737           0.816

Summary of ASCI peformance statewide:

``` r
perf(allscr)
```

    ## # A tibble: 18 x 9
    ##    grp     ind   typ        cls     ave    fst prc_amg prc_wth res_tst
    ##    <chr>   <chr> <chr>      <chr> <dbl>  <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 diatoms oe    Predictive rc    1.03   0.904  0.150   0.0963   11.6 
    ##  2 diatoms oe    Predictive rv    0.961  0.786  0.225   0.0822    2.40
    ##  3 hybrid  oe    Predictive rc    1.04   1.60   0.160   0.0800   10.9 
    ##  4 hybrid  oe    Predictive rv    1.03   2.29   0.186   0.109     5.49
    ##  5 sba     oe    Predictive rc    1.01   2.28   0.351   0.193    12.8 
    ##  6 sba     oe    Predictive rv    0.737  1.66   0.367   0.178     3.60
    ##  7 diatoms mmi   Predictive rc    0.738 30.2    0.0942  0.0467   33.1 
    ##  8 diatoms mmi   Predictive rv    0.749  9.23   0.0798  0.0412   24.4 
    ##  9 hybrid  mmi   Predictive rc    0.573  2.57   0.0996  0.0539   29.2 
    ## 10 hybrid  mmi   Predictive rv    0.587  1.33   0.0987  0.0399   19.7 
    ## 11 sba     mmi   Predictive rc    0.673  1.63   0.129   0.0953   17.0 
    ## 12 sba     mmi   Predictive rv    0.667  1.05   0.138   0.0743    9.72
    ## 13 diatoms oe    Null       rc    0.983 11.4    0.199   0.0853    4.29
    ## 14 diatoms oe    Null       rv    0.963  2.42   0.237   0.0915    1.36
    ## 15 hybrid  oe    Null       rc    1.06   1.28   0.184   0.0870    8.47
    ## 16 hybrid  oe    Null       rv    1.03   1.09   0.209   0.0852    3.05
    ## 17 sba     oe    Null       rc    0.889  5.71   0.388   0.174    10.0 
    ## 18 sba     oe    Null       rv    0.768  1.71   0.400   0.183     4.10
