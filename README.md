
<!-- README.md is generated from README.Rmd. Please edit that file -->

# exactLTRE

<!-- badges: start -->
<!-- badges: end -->

The goal of exactLTRE is to provide a set of friendly tools to increase
use of exact LTRE methods. Life Table Response Experiments are a
valuable tool for advancing our knowledge about how stage-specific vital
rates influence population dynamics. The results of LTRE’s also inform
conservation goals, by identifying the life stage and vital rate (for
example, juvenile survival or adult reproductive success) that have the
largest impact on population growth rate.

The standard methods of LTRE, presented in the Caswell (2001) textbook,
are based on approximations and are vulnerable to mistakes both due to
the approximation of specific terms and the fact that higher-order terms
are neglected from the analyses. Exact LTRE avoids errors of
interpretation by directly calculating the impact on population growth
rate of differences or variation in vital rates. With the improved
computing power that we have today, we no longer need to rely on
approximation methods that can be faulty.

This package includes functions for the approximate methods as well as
the exact methods that we introduce. This will allow researchers to use
both the methods that are standard in the literature up until now
(2021), and the new methods that we encourage as a new standard.

## Installation

You can install the development version of exactLTRE from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("chrissy3815/exactLTRE")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(exactLTRE)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
