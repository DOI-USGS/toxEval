shinyUI(fluidPage(
  
  titlePanel("toxEval"),
  fluidRow(column(6, selectInput("groupCol", label = "Column to group", 
                                choices = names(endPointInfo),
                                selected = names(endPointInfo)[10], multiple = FALSE)),
           column(6, uiOutput("groupControl"))),

    fluidRow(
      dataTableOutput('table'),
      plotOutput("graph")
    )
  )
)