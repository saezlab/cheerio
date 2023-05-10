
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Welcome to ReHeaT

<!-- badges: start -->

[![Build
Status](https://travis-ci.com/saezlab/reheat.svg?branch=master)](https://travis-ci.com/saezlab/reheat)
<!-- badges: end -->

> the **Re**ference for the **Hea**rt failure **T**ranscriptome

### About

ReHeaT is a user friendly portal to explore, analyze and download the
transcriptional consensus signature of human end-stage heart failure
(HF-CS).

The HF-CS is based on a
[meta-analysis](https://www.medrxiv.org/content/10.1101/2020.05.23.20110858v1)
of 16 independent studies. Corresponding code is available
[here](https://github.com/saezlab/HF_meta-analysis).

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
