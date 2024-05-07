FROM code.usgs.gov:5001/water/dockerr:latest

# Necessary R libraries
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
  r-cran-dataretrieval \
  r-cran-tidyverse \
  r-cran-checkmate \
  r-cran-openxlsx \
  r-cran-here \
  r-cran-shinydashboard \
  r-cran-shiny \
  r-cran-readxl \
  r-cran-rcolorbrewer \
  r-cran-shinyace \
  r-cran-rmarkdown \
  r-cran-knitr \
  r-cran-connectapi \
  r-cran-dt \
  r-cran-leaflet \
  r-cran-rmarkdown \
  && rm -rf /var/lib/apt/lists/*



