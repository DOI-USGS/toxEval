library(toxEval)
endPointInfo <- endPointInfo
choicesPerGroup <- apply(endPointInfo[,-3], 2, function(x) length(unique(x)))
groupChoices <- paste0(names(choicesPerGroup)," (",choicesPerGroup,")")
pathToApp <- system.file("extdata", package="toxEval")
summaryFile <- readRDS(file.path(pathToApp,"summary.rds"))

shinyUI(
  fluidPage(
  
    fluidRow(
      column(1),
      column(4,
        h1("toxEval"),
        selectInput("data", label = "Data", 
                         choices = c("Water Sample",
                                     "Passive Samples"),
                         selected = "Water Sample + ToxCast", multiple = FALSE),
        selectInput("sites", label = "Site", 
                                      choices = c("All","Potential 2016",summaryFile$site),
                                      selected = "All", multiple = FALSE),
        selectInput("groupCol", label = "Annotation (# choices)", 
                                      choices = setNames(names(endPointInfo)[-3],groupChoices),
                                      selected = names(endPointInfo)[20], multiple = FALSE),
        tags$div(class="header", checked=NA,
                 tags$p("For annotation information, see: "),
                 tags$a(href="http://www.epa.gov/ncct/toxcast/files/ToxCast%20Assays/ToxCast_Assay_Annotation_Data_Users_Guide_20141021.pdf", 
                        "ToxCast"))
        ),
      column(6, 
             leaflet::leafletOutput("mymap"),
             h5("Size range represents number of collected samples from 1-64")
      ),
      column(1)),
        
    fluidRow(
      column(1),
      column(10,
        tabsetPanel(
            tabPanel("Annotation Summary",
                htmlOutput("TableHeaderColumns"),
                radioButtons("radioMaxGroup", label = "",inline = TRUE,
                                  choices = list("Choice" = 1, "Chemical" = 2, "Class" = 3), 
                                  selected = 1),
                tabsetPanel(
                  tabPanel("Visualizations",
                           plotOutput("stackBarGroup"),
                           h4("All EARs"),
                           plotOutput("graphGroup")),
                  tabPanel("EAR Summary",
                           h5("maxEAR = Maximum summation of EARs per sample"),
                           h5("freq = Fraction of samples with hits"),
                          DT::dataTableOutput('tableSumm')),
                  tabPanel("Chemical Summary",
                           htmlOutput("nGroup"),
                           DT::dataTableOutput('tableGroupSumm'))
                  )
            ),
            tabPanel("Group Summary",
                 fluidRow(
                   column(5,uiOutput("groupControl")),
                   column(5,radioButtons("radio", label = "", inline = TRUE,
                                       choices = list("Chemical" = 1, "Class" = 2), 
                                       selected = 1))
                 ),
                 tabsetPanel(
                   tabPanel("Visualizations",
                            htmlOutput("BoxHeader"),
                            h5("Only shading EARs with hits (> 0.1)"),
                            plotOutput("stackBar"),
                            h5("All EARs"),
                            plotOutput("graph")
                            
                   ),
                   tabPanel("Table",
                            htmlOutput("TableHeader"),
                            DT::dataTableOutput('table'))

                   )
                 )
            
          )
        ),
      column(1)
    )
    
  )
)