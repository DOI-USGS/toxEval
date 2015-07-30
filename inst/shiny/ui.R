library(toxEval)
endPointInfo <- endPointInfo
choicesPerGroup <- apply(endPointInfo[,-3], 2, function(x) length(unique(x)))
groupChoices <- paste0(names(choicesPerGroup)," (",choicesPerGroup,")")

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
                                      choices = c("All",summary$site),
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
                tabsetPanel(
                    tabPanel("EAR Summary",
                             htmlOutput("TableHeaderColumns"),
                             DT::dataTableOutput('tableSumm')),
                    tabPanel("Chemical Summary",
                             htmlOutput("TableHeaderColumns2"),
                             DT::dataTableOutput('tableGroupSumm'))
                  )
  
            ),
            tabPanel("Group Summary",
                 tabsetPanel(
                   tabPanel("Table",
                            htmlOutput("TableHeader"),
                            DT::dataTableOutput('table')),
                   tabPanel("Box Plots",
                              htmlOutput("BoxHeader"),
                              plotOutput("stackBar"),
                              plotOutput("graph")
                            )
                 )
            )
          )
        )
  )
)