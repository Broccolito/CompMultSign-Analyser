library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyWidgets)
library(scales)
library(dplyr)
library(pracma)
library(ggplot2)
library(plotly)

source("peak_detection.R")

ui = dashboardPage(skin = "black",title = "CompMultSign Analyser",
  header = dashboardHeaderPlus(title = "CompMultSign Analyser",
                               titleWidth = 400),
  sidebar = dashboardSidebar(width = 400,div(
    sidebarMenu(
      menuItem("Chromosome",
               tabName = "Chromosome",
               startExpanded = TRUE,
               div(
                 numericInput(
                   inputId = "chr",
                   label = "Enter Chromosome Number: ",
                   value = 1,
                   min = 1,
                   max = 22,
                   step = 1,
                   width = NULL
                 ),
                 br()
               )),
      menuItem("Signal Processing",
               tabName = "Signal Processing",
               startExpanded = FALSE,
               div(
                 pickerInput(
                   inputId = "exponential",
                   label = "Exponential Enhancer: ",
                   choices = c(1.1,2,exp(1),10),
                   selected = exp(1)
                 ),
                 numericInput(
                   inputId = "kernal_size",
                   label = "Kernal Size: ",
                   value = 3,
                   min = 1,
                   max = 30,
                   step = 1,
                   width = NULL
                 ),
                 pickerInput(
                   inputId = "smooth_method",
                   label = "Smooth Method", 
                   choices = c("3RS3R", "3RSS", "3RSR", "3R", "3", "S"),
                   selected = "3RS3R",
                   options = list(
                     `live-search` = TRUE)
                 ),
                 prettySwitch(
                   inputId = "threshold",
                   label = "Select Peak by Thresholding: ",
                   value = FALSE
                 ),
                 numericInput(
                   inputId = "threshold_percentage",
                   label = "Threshold Percentage: ",
                   value = 0.95,
                   min = 0,
                   max = 1,
                   step = 0.01,
                   width = NULL
                 ),
                 prettySwitch(
                   inputId = "ranking",
                   label = "Select Peak by Ranking: ",
                   value = TRUE
                 ),
                 numericInput(
                   inputId = "top_ranking",
                   label = "Only Select Top: ",
                   value = 50,
                   min = 1,
                   max = 150,
                   step = 1,
                   width = NULL
                 )
               )),
      menuItem("Visualization", 
               tabName = "Visualization",
               startExpanded = FALSE,
               div(
                 textInput(
                   inputId = "lower_limit",
                   label = "Lower Limit",
                   value = "0"
                 ),
                 textInput(
                   inputId = "upper_limit",
                   label = "Upper Limit",
                   value = "1e20"
                 ),
                 pickerInput(
                   inputId = "peak_color",
                   label = "Peak Color", 
                   choices = colors(),
                   selected = "coral",
                   options = list(
                     `live-search` = TRUE)
                 ),
                 pickerInput(
                   inputId = "peak_shape",
                   label = "Peak Shape", 
                   choices = 0:25,
                   selected = 21,
                   options = list(
                     `live-search` = TRUE)
                 ),
                 sliderInput(
                   inputId = "peak_size",
                   label = "Peak Size",
                   min = 0.1,
                   max = 3,
                   step = 0.1,
                   value = 1
                 )
               ))
    ),
    actionButton(inputId = "generate_plot",
                 label = "Update Plot")
  )),
  body = dashboardBody(
    plotlyOutput(outputId = "generated_plot"),
    br(),
    verbatimTextOutput(outputId = "status")
  )
)

server = function(input, output, session) {
  
  observeEvent(input$generate_plot,{
    output$generated_plot = renderPlotly({
      generate_processed_df(chromo = isolate(as.numeric(input$chr)),
                            exponential = isolate(as.numeric(input$exponential)),
                            kernal_size = isolate(as.numeric(input$kernal_size)),
                            smooth_method = isolate(as.character(input$smooth_method)),
                            threshold = isolate(input$threshold),
                            ranking = isolate(input$ranking),
                            threshold_percentage = isolate(as.numeric(input$threshold_percentage)),
                            top_ranking = isolate(input$top_ranking)) %>%
        plot_df(lower_limit = isolate(as.numeric(input$lower_limit)),
                upper_limit = isolate(as.numeric(input$upper_limit)),
                peak_color = isolate(as.character(input$peak_color)),
                peak_shape = isolate(as.numeric(input$peak_shape)),
                peak_size = isolate(as.numeric(input$peak_size)))
    })
    output$status = renderText({
      paste(
        "[[Signal Processing Parameters]]",
        paste0("Chromosome Numebr: ", isolate(input$chr)),
        paste0("Exponential: ", isolate(input$exponential)),
        paste0("Kernal Size: ", isolate(input$kernal_size)),
        paste0("Smooth Method: ", isolate(input$smooth_method)),
        paste0("Select Peaks by: ", 
               ifelse(isolate(input$ranking),"Ranking ",""),
               ifelse(isolate(input$threshold), "Threshold","")),
        "",
        "[[Visualization Status]]",
        paste0("Lower Limit: ", isolate(input$lower_limit)),
        paste0("Upper Limit: ", isolate(input$upper_limit)),
        paste0("Peak Color: ", isolate(input$peak_color)),
        paste0("Peak Shape: ", isolate(input$peak_shape)),
        paste0("Peak Size: ", isolate(input$peak_size)),
      sep = "\n")
    })
  })
  
  onSessionEnded(fun = function(){stopApp()})
}

shinyApp(ui, server)