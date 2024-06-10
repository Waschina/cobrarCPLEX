
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cobrarCPLEX

<!-- badges: start -->
<!-- badges: end -->

An extension of the R-package
[cobrar](https://github.com/Waschina/cobrar), which makes IBM’s ILOG
CPLEX solver accessible.

## Requirements

- Installed R-package [cobrar](https://github.com/Waschina/cobrar)
- Working CPLEX Installation version 20.1.0 or higher

## Installation

You can install the development version of cobrarCPLEX like so:

``` r
remotes::install_github("Waschina/cobrarCPLEX", configure.args="--with-cplex-dir=/path/to/cplex")
```

<mark>Note:</mark> Please replace `/path/to/cplex` with the path to your
cplex installation.

E.g., in unix systems, cplex is by default installed to something like
`/opt/ibm/ILOG/CPLEX_Studio2211/cplex`. Thus the command to install
cobrarCPLEX would be:

``` r
remotes::install_github("Waschina/cobrarCPLEX", configure.args="--with-cplex-dir=/opt/ibm/ILOG/CPLEX_Studio2211/cplex")
```

## Example

This is a basic example which shows how to use cobrarCPLEX and compares
the computation time for flux balance analysis:

``` r
library(cobrarCPLEX)
#> Loading required package: cobrar
#> Loading required package: Matrix
#> cobrar uses...
#>  - libSBML (v. 5.19.0)
#>  - glpk (v. 5.0)
#> Using cplex version 22010100

# Download the E. coli model "iML1515" from http://bigg.ucsd.edu/
modfile <- tempfile(fileext = ".xml.gz")
download.file("http://bigg.ucsd.edu/static/models/iML1515.xml.gz", modfile)
mod <- readSBMLmod(modfile)

# FBA with GLPK
COBRAR_SETTINGS("SOLVER","glpk")
print(system.time(fba(mod)))
#>    user  system elapsed 
#>   0.345   0.000   0.346

# FBA with GLPK
COBRAR_SETTINGS("SOLVER","cplex")
system.time(fba(mod))
#>    user  system elapsed 
#>   0.075   0.000   0.075
```
