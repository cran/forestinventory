<!-- README.md is generated from README.Rmd. Please edit that file -->

forestinventory
===============

The R-package `forestinventory` addresses the current interest of
combining existing forest inventory data, which are derived by field
surveys, with additional information sources such as remote sensing
data. The major benefit of these so-called *multisource inventory
methods* is the potential increase of estimation precision without an
increase in the number of expensive field surveys. Additionally, it also
allows for deriving estimates of sufficient accuracy for spatial units
where terrestrial information is scarcely available if not absent.

The aim of `forestinventory` is to facilitate the application of
multiphase forest inventories by providing an extensive set of functions
for global and small-area estimation procedures. The implementation
includes all estimators for simple and cluster sampling published by
Daniel Mandallaz between 2007 and 2014, providing point estimates, their
external- and design-based variances as well as confidence intervals.
The procedures have also been optimized for the use of remote sensing
data as auxiliary information.

Quick demo
==========

We look at the example dataset `grisons` which comes with our package:

``` r
library(forestinventory)
?grisons
```

As the help tells us, `grisons` contains the data of a twophase
inventory: We are provided with LiDAR canopy height metrics at 306
inventory locations, and at 67 subsamples we have the terrestrially
measured timber volume values. We now want to estimate the timber volume
in m<sup>3</sup>/ha within four subdomains A, B, C and D (*small
areas*).

If we only use the terrestrial information within the small areas, we
call the `onephase`-function:

``` r
op <- onephase(formula = tvol~1,
               data = grisons,
               phase_id = list(phase.col = "phase_id_2p", terrgrid.id = 2),
                 area = list(sa.col = "smallarea", areas = c("A", "B", "C", "D")))

summary(op)
#> 
#> One-phase estimation
#>  
#> Call: 
#> onephase(formula = tvol ~ 1, data = grisons, phase_id = list(phase.col = "phase_id_2p", 
#>     terrgrid.id = 2), area = list(sa.col = "smallarea", areas = c("A", 
#>     "B", "C", "D")))
#> 
#> Method used:
#> One-phase estimator
#>  
#> Estimation results:
#>  area estimate variance n2
#>     A 410.4047 1987.117 19
#>     B 461.4429 3175.068 17
#>     C 318.0091 1180.853 15
#>     D 396.8496 2290.652 16
```

We now try to increase the precision of our estimates by applying a
*twophase estimation method*, where we use the large sample of
LiDAR-metrics and a linear regression model to specify the relationship
between the remote sensing derived predictor variables and the
terrestrial timber volume:

``` r
sae.2p.uv<- twophase(formula = tvol ~ mean + stddev + max + q75, data = grisons,
                     phase_id = list(phase.col = "phase_id_2p", terrgrid.id = 2),
                     small_area = list(sa.col = "smallarea", areas = c("A", "B","C", "D"),
                                       unbiased = TRUE))

summary(sae.2p.uv)
#> 
#> Two-phase small area estimation
#>  
#> Call: 
#> twophase(formula = tvol ~ mean + stddev + max + q75, data = grisons, 
#>     phase_id = list(phase.col = "phase_id_2p", terrgrid.id = 2), 
#>     small_area = list(sa.col = "smallarea", areas = c("A", "B", 
#>         "C", "D"), unbiased = TRUE))
#> 
#> Method used:
#> Extended pseudosynthetic small area estimator
#>  
#> Regression Model:
#> tvol ~ mean + stddev + max + q75 + smallarea
#> 
#> Estimation results:
#>  area estimate ext_variance g_variance  n1 n2 n1G n2G r.squared
#>     A 391.1605     995.5602   1016.956 306 67  94  19 0.6526503
#>     B 419.6746    1214.6053   1019.270 306 67  81  17 0.6428854
#>     C 328.0117     916.2266   1035.091 306 67  66  15 0.6430018
#>     D 371.0596    1272.7056   1112.735 306 67  65  16 0.6556178
```

We now want to compare the results and performances of the onephase and
twophase method. For such issues, the package provides the `estTable()`
function that concatenates the results from the different methods in one
`list`:

``` r
sae.table<- estTable(est.list = list(op, sae.2p.uv), sae = TRUE)

data.frame(sae.table[c(1:6,9)])
#>    area    domain   method       estimator      vartype estimate error
#> 1     A    global onephase        onephase     variance 410.4047 10.86
#> 2     A smallarea twophase psynth extended ext_variance 391.1605  8.07
#> 3     A smallarea twophase psynth extended   g_variance 391.1605  8.15
#> 4     B    global onephase        onephase     variance 461.4429 12.21
#> 5     B smallarea twophase psynth extended ext_variance 419.6746  8.30
#> 6     B smallarea twophase psynth extended   g_variance 419.6746  7.61
#> 7     C    global onephase        onephase     variance 318.0091 10.81
#> 8     C smallarea twophase psynth extended ext_variance 328.0117  9.23
#> 9     C smallarea twophase psynth extended   g_variance 328.0117  9.81
#> 10    D    global onephase        onephase     variance 396.8496 12.06
#> 11    D smallarea twophase psynth extended ext_variance 371.0596  9.61
#> 12    D smallarea twophase psynth extended   g_variance 371.0596  8.99
```

We can already see that the estimation errors of the twophase estimation
are up to 5% smaller than the onephase errors.

The function `mphase.gain()` can now be used to further compare the
performance of the methods:

``` r
mphase.gain(sae.table)
#>   area var_onephase var_multiphase   method       estimator gain  rel.eff
#> 1    A     1987.117       1016.956 twophase psynth extended 48.8 1.953986
#> 2    B     3175.068       1019.270 twophase psynth extended 67.9 3.115041
#> 3    C     1180.853       1035.091 twophase psynth extended 12.3 1.140821
#> 4    D     2290.652       1112.735 twophase psynth extended 51.4 2.058579
```

The column `gain` tells us that the twophase estimation procedure here
leads to a 67.9 % reduction in variance compared to the one- phase
procedure“. The column `rel.eff` speciﬁes the relative efﬁciency that
can be interpreted as the relative sample size of the one-phase
estimator needed to achieve the variance of the multi-phase (here
twophase) estimator. For small area”B" we can thus see that we would
have to increase the terrestrial sample size by factor 3 in the
one-phase approach in order to get the same estimation precision as the
twophase extended psynth estimator.

So in our short example, we were able to considerably improve the
estimation precision when combining the terrestrial data with the remote
sensing data. The package `forestinventory` offers further estimators
that can be applied to a wide range of multiphase inventory scenarios.

Installation
============

The package can be installed from CRAN:

``` r
install.packages("forestinventory")
```
