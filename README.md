# CERMBlidar
A collection of R functions to process aerial LIDAR data

This R package provides functions that extend those in the [lidR package](https://github.com/Jean-Romain/lidR).
It is developed mainly for use by the Centre for Environmental Risk Management of Bushfires, University of Wollongong, 
Australia.

To install the package on your system, first install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) 
if you do not already have it (required to compile C++ functions included in the package).

Next enter the following commands from the R console:

```
# install devtools if not already present
install.packages("devtools")

# install this package
devtools::install_github("mbedward/CERMBlidar")
```
The package may be freely used by anyone, but please note that: (1) Many functions are particular to data 
conventions used by UOW and New South Wales State Government and might not be suitable for general use; 
(2) The package is at an early stage of development and features are likely to change without notice between 
successive 0.x versions.
