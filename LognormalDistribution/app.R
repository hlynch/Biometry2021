#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#to deploy app for the first time, set folder with correct app.R file as the working directory then run the two lines below
#library(rsconnect)
#deployApp()

library(shiny)

# Define UI for application that draws a histogram
server <- function(input,output){
  
  ####################################### 
  # These are the reactive calculations # These do not get displayed 
  ####################################### 
  
  rand_samp <- reactive({
    rlnorm(input$obs, input$mu, input$sigma)
  })
  
  # calculate some stuff needed to define the quantiles 
#  min_quant <- reactive({
#    input$mean - 3*input$sd # set the minimum of the quantiles 
#  })
  
#  max_quant <- reactive({
#    input$mean + 3*input$sd # set the maximum of the quantiles 
#  })  
  
  quantiles <- reactive({
    seq(0, 20, 0.01) 
  })
  
  prob_dens <- reactive({
    dlnorm(quantiles(), input$mu, input$sigma) 
  })
  
  cuml_prob <- reactive({
    plnorm(quantiles(), input$mu, input$sigma) 
  })
  
  quant_plot <- reactive({
    qlnorm(seq(0,1,0.01), input$mu, input$sigma)
  })
  
  
  
  ################################## 
  # These are the reactive outputs # These get displayed
  ################################## 
  
  # Plot the random sample 
  output$Plot_sample <-renderPlot({
    plot(rand_samp(),main="",xlab="Sample index",ylab="Value",pch=16,col=rgb(225/255,74/255,141/255,0.5),ylim=c(0,40))
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
    
    sliderInput("mu",
                "mu:",
                min=-5,
                max=5,
                value=0,
                step=0.1),
    
    sliderInput("sigma",
                "sigma:",
                value=1,
                min=0.1,
                max=3,
                step=0.1)
    
  ),  # GGPLOT
  
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
               h5("dlnorm(quantiles, mu, sigma)")
      ),
      tabPanel("Cumulative Probability",
               plotOutput("Plot_cuml_prob"),
               h5("plnorm(quantiles, mu, sigma)")
      ),
      tabPanel("Quantile Plot",
               plotOutput("Plot_quantiles"),
               h5("qlnorm(seq(0,1,0.01), mu, sigma)")
      )  
    )
  )
  
))

# Run the application 
shinyApp(ui = ui, server = server)

