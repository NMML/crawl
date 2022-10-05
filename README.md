<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![crawl status
badge](https://dsjohnson.r-universe.dev/badges/crawl)](https://dsjohnson.r-universe.dev)
[![R-CMD-check](https://github.com/NMML/crawl/workflows/R-CMD-check/badge.svg)](https://github.com/NMML/crawl/actions)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

## Correlated RAndom Walk Library of R functions

The Correlated RAndom Walk Library of R functions was designed for
fitting continuous-time correlated random walk (CTCRW) models with time
indexed covariates. The model is fit using the Kalman-Filter on a state
space version of the continuous-time stochastic movement process. The
use case in the estimation of animal movement paths where the observed
locations are determined from Argos or FastLoc enabled bio-loggers. In
addition to the continuous-time component, `{crawl}` was specifically
developed to incorporate the measurement error often associated with
these observed locations. Lastly, `{crawl}` provides a framework for
multiple imputation workflows for incorporation of model uncertainty.

## The Future of crawl

The original code base and concepts for `{crawl}` were developed almost
15 years ago. Much has changed in the world of movement ecology, spatial
statistics, R, bio-logging, and many other fields. In some cases, we’ve
done a fairly good job keeping pace; in other cases, we’ve fallen
behind. We feel it is time for a new approach and will, from now on, be
focusing our development efforts on `{crawl2}`. We will continue to
maintain `{crawl}`, improve the documentation, and ensure compatibility
with dependent packages.

## Installation

## Install via CRAN

`{crawl}` is currently available on CRAN and R \>= 4.0 is highly
recommended.

``` r
# install latest version of crawl from CRAN
install.packages("crawl")
```

### Install via R-Universe

The latest version of `{crawl}` is also available via R-Universe.

``` r
# Install crawl from my R-Universe repository
# Enable repository from dsjohnson

options(repos = c(
  dsjohnson = 'https://dsjohnson.r-universe.dev',
  CRAN = 'https://cloud.r-project.org'))

# Download and install crawl in R
install.packages('crawl')

# Browse the crawl manual pages
help(package = 'crawl')
```

You can also add the repository to your local list of repositories in
your *.Rprofile* and this will ensure `update.packages()` pulls any new
releases of `{crawl}` from R-Universe

``` r
#install.packages("usethis")
usethis::edit_r_profile()

# add the following text or replace existing repos option

options(repos = c(dsjohnson = 'https://dsjohnson.r-universe.dev',
                  CRAN = 'https://cloud.r-project.org'))
```

### Install via Github

A development version of `{pathroutr}` is also available from
[GitHub](https://github.com/NMML/crawl). This version should be used
with caution and only after consulting with package authors.

``` r
# install.packages("remotes")
remotes::install_github("NMML/crawl@devel")
```

### Disclaimer

This repository is a scientific product and is not official
communication of the National Oceanic and Atmospheric Administration, or
the United States Department of Commerce. All NOAA GitHub project code
is provided on an ‘as is’ basis and the user assumes responsibility for
its use. NOAA and DOC have relinquished control of the information and
no longer has responsibility to protect the integrity, confidentiality,
or availability of the information. Any claims against the Department of
Commerce or Department of Commerce bureaus stemming from the use of this
GitHub project will be governed by all applicable Federal law. Any
reference to specific commercial products, processes, or services by
service mark, trademark, manufacturer, or otherwise, does not constitute
or imply their endorsement, recommendation or favoring by the Department
of Commerce. The Department of Commerce seal and logo, or the seal and
logo of a DOC bureau, shall not be used in any manner to imply
endorsement of any commercial product or activity by DOC or the United
States Government.
