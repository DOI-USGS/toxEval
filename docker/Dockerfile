FROM code.chs.usgs.gov:5001/ctek/docker/r-lang/r-base:4.4

# Necessary R libraries
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
  r-cran-devtools \
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



