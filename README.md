
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Welcome to CHEERIO

<!-- badges: start -->

[![Build
Status](https://travis-ci.com/saezlab/reheat.svg?branch=master)](https://travis-ci.com/saezlab/reheat)
<!-- badges: end -->

> Cardiac Hypertrophy Transcriptome Appication

### About

CHEERIO (Cardiac HypErtrophy gEne expRessIOn) database is a comprehensive resource of cardiac gene expression during hypertrophic growth, supporting researchers in elaborating, testing, and refining their hypotheses on the molecular nature of cardiac hypertrophy assisted through a free and easy-to-use web application. 

This is a collaborative project with the Mirko Völkers AG, Heidelberg University. 


### How to access

There are basically three ways of how to access Cheerio:

  - You can access a live version running
    [here](https://jlanzer.shinyapps.io/shiny_hypertophy/) on the server from
    `shinyapps.io`.

  - You can run the app locally in an interactive R session. Before make
    sure you have all packages installed listed in
    [`sub/global.R`](https://github.com/saezlab/cheerio/sub/global.R).

<!-- end list -->

``` r
shiny::runGitHub("cheerio", "saezlab")
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
