library(DT)
library(dplyr)
library(toxEval)
endPointInfo <- endPointInfo
wData <- wData
pCodeInfo <- pCodeInfo

packagePath <- system.file("extdata", package="toxEval")
filePath <- file.path(packagePath, "stationINFO.RData")
load(file=filePath)

endPointInfo <- endPointInfo
summary <- readRDS("data/summary.rds")
endPoint <- readRDS("data/endPoint.rds")
chemicalSummary <- readRDS("data/chemicalSummary.rds")

shinyUI(fluidPage(
  
  titlePanel("toxEval"),
  
  
  sidebarLayout(
    sidebarPanel(
      selectInput("groupCol", label = "Column to group", 
                  choices = names(endPointInfo),
                  selected = names(endPointInfo)[11], multiple = FALSE),
      uiOutput("groupControl")
      # submitButton("Submit")
    ),
    mainPanel(
#       dygraphOutput("dygraph1"),
#       dygraphOutput("dygraph2"),
#       dygraphOutput("dygraph3")
    )
  )
))