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
           ),

    tabsetPanel(
      tabPanel("Group Summary",
               fluidRow(
                 fluidRow(h4("Table of summations summaries:", style  = "text-align:center")),
                 DT::dataTableOutput('table'),
                 fluidRow(h4("Boxplot of summations:", style  = "text-align:center")),
                 plotOutput("graph")
               )
      ),
      tabPanel("Column Summary", DT::dataTableOutput('tableSumm'))
    ),
    fluidRow(
      leaflet::leafletOutput("mymap")
    )
  )
)