shinyUI(fluidPage(
  
  titlePanel("toxEval"),
  fluidRow(column(width = 3, selectInput("sites", label = "Site", 
                                         choices = c("All",summary$site),
                                         selected = "All", multiple = FALSE)),
           column(width = 4, selectInput("groupCol", label = "Column to group (# choices)", 
                                choices = setNames(names(endPointInfo)[-3],groupChoices),
                                selected = names(endPointInfo)[20], multiple = FALSE)),
           column(width = 4, uiOutput("groupControl"))
           # column(width = 4, actionButton("calculate","Calculate"))
           # column(width = 4, dataTableOutput('groupSumTable'))
           ),

    fluidRow(
      DT::dataTableOutput('table'),
      plotOutput("graph"),
      leaflet::leafletOutput("mymap")
    )
  )
)