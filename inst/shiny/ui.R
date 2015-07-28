library(toxEval)
endPointInfo <- endPointInfo
choicesPerGroup <- apply(endPointInfo[,-3], 2, function(x) length(unique(x)))
groupChoices <- paste0(names(choicesPerGroup)," (",choicesPerGroup,")")

shinyUI(
  
  navbarPage("toxEval",
    tabPanel("Water Sample + ToxCast",
      fluidRow(
        column(4,
          selectInput("sites", label = "Site", 
                                        choices = c("All",summary$site),
                                        selected = "All", multiple = FALSE)),
        column(4,
          selectInput("groupCol", label = "Annotation (# choices)", 
                                        choices = setNames(names(endPointInfo)[-3],groupChoices),
                                        selected = names(endPointInfo)[20], multiple = FALSE)),
        column(4, 
          tags$div(class="header", checked=NA,
                   tags$p("See: "),
                   tags$a(href="http://www.epa.gov/ncct/toxcast/files/ToxCast%20Assays/ToxCast_Assay_Annotation_Data_Users_Guide_20141021.pdf", "ToxCast"),
                   tags$p("for annotation information")
          ))
    
        ),
        
      fluidRow(
          tabsetPanel(
              tabPanel("Annotation Summary", 
#                        tabsetPanel(
#                          tabPanel("Table",
                                  tabsetPanel(
                                      tabPanel("EAR Summary",
                                               htmlOutput("TableHeaderColumns"),
                                               DT::dataTableOutput('tableSumm')),
                                      tabPanel("Chemical Summary",
                                               htmlOutput("TableHeaderColumns2"),
                                               DT::dataTableOutput('tableGroupSumm'))
                                    )
                                  # )
#                          tabPanel("Box Plot",
#                                   htmlOutput("BoxHeaderColumns"),
#                                   h3("Maybe grid")),
#                          tabPanel("Map",
#                                   h3("Hmm...")
#                                   )
#                       )
              ),
              tabPanel("Group Summary",
                       uiOutput("groupControl"),
                       tabsetPanel(
                         tabPanel("Table",
                                  htmlOutput("TableHeader"),
                                  DT::dataTableOutput('table')),
                         tabPanel("Box Plots",
                                    htmlOutput("BoxHeader"),
                                    plotOutput("stackBar"),
                                    plotOutput("graph")
                                  ),
                         tabPanel("Map",
                                  h3(""),
                                  leaflet::leafletOutput("mymap"))
                       )
              )
            )
          )
        
    ),
    tabPanel("Passive Samples + ToxCast"),
    tabPanel("Water Samples + Benchmarks")
  )
)