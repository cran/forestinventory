<!-- README.md is generated from README.Rmd. Please edit that file -->
forestinventory
===============

The R-package `forestinventory` adresses the current interest of combining existing forest inventory data, which are derived by field surveys, with additional information sources such as remote sensing data. The major benefit of these so-called *multisource inventory methods* is the potential increase of estimation precision without an increase in the number of expensive field surveys. Additionally, it also allows for deriving estimates of sufficient accuracy for spatial units where terrestrial information is scarcely available if not absent.

The aim of `forestinventory` is to facilitate the application of multiphase forest inventories by providing an extensive set of functions for global and small-area estimation procedures. The implementation includes all estimators for simple and cluster sampling published by Daniel Mandallaz between 2007 and 2014, providing point estimates, their external- and design-based variances as well as confidence intervals. The procedures have also been optimized for the use of remote sensing data as auxiliary information.

Quick demo
==========

We look at the example dataset `grisons` which comes with our package:

``` r
library(forestinventory)
?grisons
```

As the help tells us, `grisons` contains the data of a twophase inventory: We are provided with LiDAR canopy height metrics at 306 inventory locations, and at 67 subsamples we have the terrestrially measured timber volume values. We now want to estimate the timber volume in m<sup>3</sup>/ha within four subdomains A, B, C and D (*small areas*).

If we only use the terrestrial information within the small areas, we call the `onephase`-function:

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

We compute the estimation error in % as the *Standard Deviation / Point Estimate*:

``` r
ee.o<- 100*(sqrt(op$estimation$variance) / op$estimation$estimate)
names(ee.o)<- op$estimation$area
ee.o
#>        A        B        C        D 
#> 10.86174 12.21120 10.80583 12.06018
```

We now try to increase the precision of our estimates by applying a *twophase estimation method*, where we use the large sample of LiDAR-metrics and a linear regression model to specify the relationship between the remote sensing derived predictor variables and the terrestrial timber volume:

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

We again compute the estimation error in % and see how the efficiency has been increased from the onephase to the twophase approach:

``` r
ee.uv<- 100*(sqrt(sae.2p.uv$estimation$g_variance) / sae.2p.uv$estimation$estimate)
names(ee.uv)<- op$estimation$area
ee.uv
#>        A        B        C        D 
#> 8.152598 7.607323 9.808440 8.989844

efficiency<- 100* ((ee.o - ee.uv) / ee.o)
names(efficiency)<- op$estimation$area
efficiency
#>         A         B         C         D 
#> 24.942062 37.702103  9.230142 25.458441
```

In our short example, we were able to reduce the estimation error by 9 % up to 37 % when combining the terrestrial data with the remote sensing data. The package `forestinventory` offers further estimators that can be applied to a wide range of multiphase inventory scenarios.

Installation
============

The package will soon be available from CRAN:

``` r
install.packages("forestinventory")
```

Coming soon ...
===============

A vignette will be released as soon as possible, giving detailed information about the implemented estimation methods and estimators and how to apply them to inventory scenarios in practice.
