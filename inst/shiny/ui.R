shinyUI(
  
  navbarPage("toxEval",
    tabPanel("Water Sample + ToxCast",
    sidebarLayout(
      sidebarPanel(
        selectInput("sites", label = "Site", 
                                      choices = c("All",summary$site),
                                      selected = "All", multiple = FALSE),
        selectInput("groupCol", label = "Annotation (# choices)", 
                                      choices = setNames(names(endPointInfo)[-3],groupChoices),
                                      selected = names(endPointInfo)[20], multiple = FALSE),
        
        tags$div(class="header", checked=NA,
                 tags$p("See: "),
                 tags$a(href="http://www.epa.gov/ncct/toxcast/files/ToxCast%20Assays/ToxCast_Assay_Annotation_Data_Users_Guide_20141021.pdf", "ToxCast"),
                 tags$p("for annotation information")
        )
  
      ),
      
      mainPanel(
        tabsetPanel(
            tabPanel("Annotation Summary", 
                     tabsetPanel(
                       tabPanel("Table",
                                htmlOutput("TableHeaderColumns"),
                                DT::dataTableOutput('tableSumm')),
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
                                leaflet::leafletOutput("mymap"))
                     )
            )
          )
        )
      )
    ),
    tabPanel("Passive Samples + ToxCast"),
    tabPanel("Water Samples + Benchmarks")
  )
)