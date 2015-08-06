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
                                      choices = c("All",summaryFile$site),
                                      selected = "All", multiple = FALSE),
        selectInput("groupCol", label = "Annotation (# choices)", 
                                      choices = setNames(names(endPointInfo)[-3],groupChoices),
                                      selected = names(endPointInfo)[20], multiple = FALSE),
        tags$div(class="header", checked=NA,
                 tags$p("For annotation information, see: "),
                 tags$a(href="http://www.epa.gov/ncct/toxcast/files/ToxCast%20Assays/ToxCast_Assay_Annotation_Data_Users_Guide_20141021.pdf", 
                        "ToxCast")),
        radioButtons("radio", label = "",
                     choices = list("Chemical" = 1, "Class" = 2), 
                     selected = 1)
        ),
      column(6, 
             leaflet::leafletOutput("mymap")
      ),
      column(1)),
        
    fluidRow(
      column(1),
      column(10,
        tabsetPanel(
            tabPanel("Annotation Summary",
                h2("Information about Annotations"),
                tabsetPanel(
                  tabPanel("Visualizations",
                           htmlOutput("BoxHeader2"),
                           h4("Only shading EARs with hits (> 0.1)"),
                           plotOutput("stackBarGroup"),
                           h4("All EARs"),
                           plotOutput("graphGroup")),
                  tabPanel("EAR Tally Summary",
                             htmlOutput("TableHeaderColumns"),
                             h4("maxEAR = Maximum (summation of EARs per sample)"),
                             h4("freq = Fraction of samples with hits"),
                             DT::dataTableOutput('tableSumm')),
                    tabPanel("Chemical Tally Summary",
                             htmlOutput("TableHeaderColumns2"),
                             htmlOutput("maxGroup"),
                             htmlOutput("meanGroup"),
                             htmlOutput("nGroup"),
                             DT::dataTableOutput('tableGroupSumm'))
                  )
  
            ),
            
            tabPanel("Group Summary",
                 h2("Information about Groups"),
                 uiOutput("groupControl"),
                 tabsetPanel(
                   tabPanel("Visualizations",
                            htmlOutput("BoxHeader"),
                            h4("Only shading EARs with hits (> 0.1)"),
                            plotOutput("stackBar"),
                            h4("All EARs"),
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