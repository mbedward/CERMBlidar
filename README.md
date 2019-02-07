# CERMBlidar
A collection of R functions to process aerial LIDAR data

This R package provides functions that extend those in the [lidR package](https://github.com/Jean-Romain/lidR).
It is developed mainly for use by the Centre for Environmental Risk Management of Bushfires, University of Wollongong, 
Australia.

To install the package on your system, first install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) if you do not already have it (required to compile C++ functions included in the package).

Next enter the following commands from the R console:

```
# install devtools if not already present
install.packages("devtools")

# install this package
devtools::install_github("mbedward/CERMBlidar")
```
