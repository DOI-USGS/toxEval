# toxEval

Initial code for studying ToxCast data in relation to measured concentrations.

## Installation of R and RStudio

This section should only need to be done once per computer.

The following link walks you through an installation of R and RStudio:
[Installation Instructions](https://owi.usgs.gov/R/training-curriculum/intro-curriculum/Before/)

If you follow those instructions exactly, you should have the USGS R repository (GRAN) added to your R profile. If that step doesn't ring a bell, paste the following into your R console:

```r
rprofile_path = file.path(Sys.getenv("HOME"), ".Rprofile")
write('\noptions(repos=c(getOption(\'repos\'),
    CRAN=\'https://cloud.r-project.org\',
    USGS=\'https://owi.usgs.gov/R\'))\n',
      rprofile_path, 
      append =  TRUE)

cat('Your Rprofile has been updated to include GRAN.
    Please restart R for changes to take effect.')
```

*RESTART RSTUDIO!*

Useful links:
[Download R Windows](https://cran.r-project.org/bin/windows/base/)
[Download R Mac](https://cran.r-project.org/bin/macosx/)
[Download RStudio](https://www.rstudio.com/products/rstudio/download/)


## Installation of toxEval

This section should also only have to be done once. It assumes the USGS R repository (GRAN) was added to your R profile as described above.

```r
install.packages("toxEval")
```

Regularly, it is a good idea to update *ALL* your packages in R. If using RStudio, this is quite easy, there's an Update button in the "Packages" tab. This checks CRAN and GRAN for updates. It is a good idea to click this update regularly.

![update](images/update.png)

## Run toxEval

To run the toxEval app:

1. Open RStudio
2. In the Console (lower-left window of RStudio) paste the following:

```r
library(toxEval)
explore_endpoints()

```

## At-your-own-risk, highly developmental, probably buggy code installation instructions

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
