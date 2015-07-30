library(toxEval)
endPointInfo <- endPointInfo
choicesPerGroup <- apply(endPointInfo[,-3], 2, function(x) length(unique(x)))
groupChoices <- paste0(names(choicesPerGroup)," (",choicesPerGroup,")")
pathToApp <- system.file("extdata", package="toxEval")
summaryFile <- readRDS(file.path(pathToApp,"summary.rds"))

shinyUI(
  fluidPage(
  
    fluidRow(
      column(4,
        h1("toxEval"),
        selectInput("data", label = "Data", 
                         choices = c("Water Sample",
                                     "Passive Samples",
                                     "Water Samples + Passive"),
                         selected = "Water Sample + ToxCast", multiple = FALSE),
        selectInput("sites", label = "Site", 
                                      choices = c("All",summaryFile$site),
                                      selected = "All", multiple = FALSE),
        selectInput("groupCol", label = "Annotation (# choices)", 
                                      choices = setNames(names(endPointInfo)[-3],groupChoices),
                                      selected = names(endPointInfo)[20], multiple = FALSE),
        tags$div(class="header", checked=NA,
                 tags$p("For annotation information, see: "),
                 tags$a(href="http://www.epa.gov/ncct/toxcast/files/ToxCast%20Assays/ToxCast_Assay_Annotation_Data_Users_Guide_20141021.pdf", "ToxCast")),
        
        uiOutput("groupControl")),
      column(8, 
             leaflet::leafletOutput("mymap")
      )),
        
    fluidRow(
        tabsetPanel(
            tabPanel("Annotation Summary",
                h2("Information about Annotations"),
                tabsetPanel(
                    tabPanel("EAR Tally Summary",
                             htmlOutput("TableHeaderColumns"),
                             DT::dataTableOutput('tableSumm')),
                    tabPanel("Chemical Tally Summary",
                             htmlOutput("TableHeaderColumns2"),
                             DT::dataTableOutput('tableGroupSumm'))
                  )
  
            ),
            
            tabPanel("Group Summary",
                 h2("Information about Groups"),
                 tabsetPanel(
                   tabPanel("Table",
                            htmlOutput("TableHeader"),
                            DT::dataTableOutput('table')),
                   tabPanel("Box Plots",
                              htmlOutput("BoxHeader"),
                              h3("Only shading EARs with hits (> 0.1)"),
                              plotOutput("stackBar"),
                              h3("All EARs"),
                              plotOutput("graph")
                            )
                 )
            )
          )
        )
  )
)