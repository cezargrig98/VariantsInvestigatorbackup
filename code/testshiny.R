library(shiny)

ui <- fluidPage(
  checkboxGroupButtons(
    inputId = "Id002",
    label = "Choices", 
    choices = c("Choice 1", "Choice 2", "Choice 3")
    
), textOutput("txt")
)
server <- function(input, output, session) {
  output$txt <- renderText({paste(input$Id002)})
}

shinyApp(ui, server)