library(DT)
library(dplyr)
library(toxEval)
endPointInfo <- endPointInfo
wData <- wData
pCodeInfo <- pCodeInfo

packagePath <- system.file("extdata", package="toxEval")
filePath <- file.path(packagePath, "stationINFO.RData")
load(file=filePath)
siteKey <- setNames(stationINFO$shortName, stationINFO$fullSiteID)

endPointInfo <- endPointInfo

# For debugging!
pathToApp <- "D:/LADData/RCode/toxEval/toxEvalUI"

summary <- readRDS(file.path(pathToApp,"data/summary.rds"))
endPoint <- readRDS(file.path(pathToApp,"data/endPoint.rds"))
chemicalSummary <- readRDS(file.path(pathToApp,"data/chemicalSummary.rds"))

shinyUI(fluidPage(
  
  titlePanel("toxEval"),
  
  
  sidebarLayout(
    sidebarPanel(
      selectInput("groupCol", label = "Column to group", 
                  choices = names(endPointInfo),
                  selected = names(endPointInfo)[10], multiple = FALSE),
      uiOutput("groupControl")
    ),
    mainPanel(
      dataTableOutput('table')
    )
  )
))