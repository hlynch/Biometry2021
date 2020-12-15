#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#library(rsconnect)
#deployApp()

library(shiny)

# Define UI for application that draws a histogram
server <- function(input,output){
  
  ####################################### 
  # These are the reactive calculations # These do not get displayed 
  ####################################### 
  
  rand_samp <- reactive({
    rnorm(input$obs, input$mean, input$sd)
  })
  
  # calculate some stuff needed to define the quantiles 
#  min_quant <- reactive({
#    input$mean - 3*input$sd # set the minimum of the quantiles 
#  })
  
#  max_quant <- reactive({
#    input$mean + 3*input$sd # set the maximum of the quantiles 
#  })  
  
  quantiles <- reactive({
    seq(-10, 10, 0.01) 
  })
  
  # probability density 
  prob_dens <- reactive({
    dnorm(quantiles(), input$mean, input$sd) 
  })
  
  # cumulative probability 
  cuml_prob <- reactive({
    pnorm(quantiles(), input$mean, input$sd) 
  })
  
  # quantiles plot 
  quant_plot <- reactive({
    qnorm(seq(0,1,0.01), input$mean, input$sd)
  })
  
  
  ################################## 
  # These are the reactive outputs # These get displayed
  ################################## 
  
  # Plot the random sample 
  output$Plot_sample <-renderPlot({
    plot(rand_samp(),main="",xlab="Sample index",ylab="Value",pch=16,col=rgb(225/255,74/255,141/255,0.5),ylim=c(-10,10))
  },
  width=400,height=400
  )
  
  # Mean 
  output$summary <- renderPrint({
    summary(rand_samp())
  })
  
  # Standard deviation 
  output$stdev <- renderPrint({
    sd(rand_samp())
  })
  
  # Plot the probability density
  output$Plot_prob_dens <-renderPlot({
    hist(rand_samp(),main="",xlab="Value",ylab="Probability density",col=rgb(225/255,74/255,141/255,0.5),freq=F)
    lines(quantiles(), prob_dens(),lwd=2)
  },
  width=400,height=400
  )
  
  # Plot the cumulative probability 
  output$Plot_cuml_prob <-renderPlot({
    plot(quantiles(), cuml_prob(), type="l", ylab="Cumulative Probability", xlab="Values")
    lines(ecdf(rand_samp()),col=rgb(225/255,74/255,141/255,0.5),lwd=2)
  },
  width=400,height=400
  )
  
  # Plot the quantiles
  output$Plot_quantiles <-renderPlot({
    plot(seq(0,1,0.01),quant_plot(),type="l",ylab="Values", xlab="Probability")
    points(seq(0,1,0.01),quantile(rand_samp(),probs=seq(0,1,0.01)),col=rgb(225/255,74/255,141/255,0.7),pch=16)
  },
  width=400,height=400
  )
  
}

ui <- shinyUI(pageWithSidebar(
  
  # Title
  headerPanel(""),
  
  sidebarPanel(
    sliderInput("obs", 
                label = "Number of observations:", 
                min = 1, 
                max = 1000,
                value = 500,
                step=1),
    
    sliderInput("mean",
                "Mean:",
                min=-10,
                max=10,
                value=0,
                step=0.1),
    
    sliderInput("sd",
                "Standard deviation:",
                value=1,
                min=0.1,
                max=10,
                step=0.1)
  ),
  
  # GGPLOT
  
  mainPanel(
    tabsetPanel(
      tabPanel("Random Sample",
               plotOutput("Plot_sample"),
               h6("Summary:"),
               verbatimTextOutput("summary"),
               h5("Standard deviation:"),
               verbatimTextOutput("stdev")
      ),
      tabPanel("Probability Density",
               plotOutput("Plot_prob_dens"),
               h5("dnorm(quantiles, mean, sd)")
      ),
      tabPanel("Cumulative Probability",
               plotOutput("Plot_cuml_prob"),
               h5("pnorm(quantiles, mean, sd)")
      ),
      tabPanel("Quantile Plot",
               plotOutput("Plot_quantiles"),
               h5("qnorm(seq(0,1,0.01), mean, sd)")
      )  
    )
  )
  
))

# Run the application 
shinyApp(ui = ui, server = server)

