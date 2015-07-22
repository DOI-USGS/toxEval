shinyUI(fluidPage(
  
  titlePanel("toxEval"),
  fluidRow(column(width = 4, selectInput("groupCol", label = "Column to group", 
                                choices = names(endPointInfo),
                                selected = names(endPointInfo)[10], multiple = FALSE)),
           column(width = 4, uiOutput("groupControl")),
           column(width = 4, actionButton("calculate","Calculate"))
           # column(width = 4, dataTableOutput('groupSumTable'))
           ),

    fluidRow(
      dataTableOutput('table'),
      plotOutput("graph"),
      leaflet::leafletOutput("mymap")
    )
  )
)