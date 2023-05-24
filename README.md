
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Welcome to CHROMO

<!-- badges: start -->

[![Build
Status](https://travis-ci.com/saezlab/reheat.svg?branch=master)](https://travis-ci.com/saezlab/reheat)
<!-- badges: end -->

> the Ch

### About

### How to access

There are basically three ways of how to access ReHeaT:

  - You can access a live version running
    [here](https://saezlab.shinyapps.io/reheat/) on the server from
    `shinyapps.io`.

  - You can run the app locally in an interactive R session. Before make
    sure you have all packages installed listed in
    [`sub/global.R`](https://github.com/saezlab/reheat/blob/master/sub/global.R).

<!-- end list -->

``` r
shiny::runGitHub("reheat", "saezlab")
```

  - You can run the app locally by cloning this repository. All required
    packages can be easily installed using the
    [`renv`](https://rstudio.github.io/renv/index.html) package.

<!-- end list -->

``` r
# install all required packages using the renv package
renv::restore()
shiny::runApp()
```
