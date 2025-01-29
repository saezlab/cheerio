
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Welcome to CHEERIO

> Cardiac HypErtrophy gEne expRessIOn database

### About

CHEERIO (Cardiac HypErtrophy gEne expRessIOn) is an integrated resource designed to facilitate the study of gene expression in cardiac hypertrophy. It provides a comprehensive collection of data from rodent and human models, enabling researchers to explore gene regulation during cardiac hypertrophy and heart failure. The platform supports hypothesis generation and testing, offering an accessible and user-friendly interface for data analysis.

By combining various datasets, CHEERIO provides insights into the regulation of the transcriptome, translatome and proteome in both physiological and pathological hypertrophy. It highlights key differences between early and advanced stages of cardiac hypertrophy, and offers a detailed comparison across species. Additionally, CHEERIO enables the identification of gene expression changes specifically associated with cardiac failure in the context of pathological hypertrophy.

### How to cite

A cross-species multimodal database of gene expression in cardiac hypertrophy
Christoph Sandmann, Jan Lanzer, Frank Stein, Ellen Malovrh, Nicholas Rudinger, Norbert Frey, Julio Saez Rodriguez, Mirko Volkers
bioRxiv 2025.01.27.635074; doi: https://doi.org/10.1101/2025.01.27.635074 

### How to access

There are three ways to access CHEERIO:

  - You can access a live version running
    [here](https://voelkerslab.shinyapps.io/cheerio/) on the server from
    `shinyapps.io`.

  - You can run the application locally in an interactive R session.
    Ensure beforehand that you have installed all packages listed in
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
