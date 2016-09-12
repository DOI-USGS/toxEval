# toxEval

Initial code for studying ToxCast data in relation to measured concentrations.

To pull the upstream changes in from github in RStudio. In the Git tab, choose the More dropdown, and choose Shell. In that window, enter:

```
git pull upstream master
```

Then, rebuild the package...go to the "Build" tab in RStudio, and choose Build and Reload.

Alternatively, use the `devtools` package. In R:

```r
devtools::install_github("USGS-R/toxEval")
```



Disclaimer
----------
This software is in the public domain because it contains materials that originally came from the U.S. Geological Survey  (USGS), an agency of the United States Department of Interior. For more information, see the official USGS copyright policy at [http://www.usgs.gov/visual-id/credit_usgs.html#copyright](http://www.usgs.gov/visual-id/credit_usgs.html#copyright)

Although this software program has been used by the USGS, no warranty, expressed or implied, is made by the USGS or the U.S. Government as to the accuracy and functioning of the program and related program material nor shall the fact of distribution constitute any such warranty, and no responsibility is assumed by the USGS in connection therewith.

This software is provided "AS IS."

Package Installation
---------------------------------

To install the `toxEval` package you need to be using R 3.0 or greater. Then use the following command:

```R
install.packages("toxEval", repos=c("http://owi.usgs.gov/R",
        "http://cran.us.r-project.org"))

library(toxEval)
explore_endpoints()

```

Linux: [![travis](https://api.travis-ci.org/USGS-R/toxEval.svg?branch=master)](https://travis-ci.org/USGS-R/toxEval)


 [
   ![CC0](http://i.creativecommons.org/p/zero/1.0/88x31.png)
 ](http://creativecommons.org/publicdomain/zero/1.0/)
