library(tidyverse)
library(shiny)
library(bslib)
library(readxl)
library(plotly)


#### Import map database 
data_map <- import(here("data base and clean", "map_brazil.rda"))  # suppose it loads 'data_map'

data_map <- data_map %>% 
  st_as_sf()

######## Dashboard
ui <- page_fluid(
  theme = bs_theme(bootswatch = "minty"),
  
  sliderInput(
    "year", "Year", 
    min = 2007, max = 2023, 
    value = 2007, 
    animate = animationOptions(interval = 2000, loop = TRUE)
  ),
  plotlyOutput("map", height = "800px", width = "100%")  # Show plot here
)

server <- function(input, output) {
  output$map <- renderPlotly({
    filtered_data <- data_map %>% 
      filter(year == input$year)
    
#breaks <- c(10, 50, 100, 150, 200, 250, 300)
    
    p1 <- ggplot(filtered_data) +
      geom_sf(aes(fill = n, color = name_muni), lwd = 0.1) + 
      labs(fill = "Cases") +
      scale_fill_viridis_c(limits = c(1, 300), option = "plasma", na.value = "white", 
                           breaks = breaks) +
      theme_minimal(base_size = 16)
    
    ggplotly(p1)
  })
}

shinyApp(ui = ui, server = server)




