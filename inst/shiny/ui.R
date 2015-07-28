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
                       tabsetPanel(
                         tabPanel("Table",
                                  tabsetPanel(
                                    tabPanel("Chemical Summary",
                                             htmlOutput("TableHeaderColumns"),
                                             DT::dataTableOutput('tableSumm')),
                                    tabPanel("Group Summary",
                                             htmlOutput("TableHeaderColumns2"),
                                             DT::dataTableOutput('tableGroupSumm')))
                                  ),
                         tabPanel("Box Plot",
                                  htmlOutput("BoxHeaderColumns")),
                         tabPanel("Map",
                                  h3("Hmm...")
                                  )
                       )
              ),
              tabPanel("Group Summary",
                       uiOutput("groupControl"),
                       tabsetPanel(
                         tabPanel("Table",
                                  htmlOutput("TableHeader"),
                                  DT::dataTableOutput('table')),
                         tabPanel("Box Plot",
                                    htmlOutput("BoxHeader"),
                                    plotOutput("graph")),
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