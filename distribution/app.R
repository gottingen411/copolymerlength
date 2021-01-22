#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
source('global.R')
# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel( h1(
        tags$b("Distribution of End-to-End Length for Co-Polymer A/B"),
        align = "center",
        style = "font-family: 'Times', serif"
    )),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            numericInput(
                "na",
                label = h4("Number of monomers A"),
                value = 5,
                min = 1,
                max = 100000,
                step = 0.1,
                width = NULL
            ),
            numericInput(
                "thetaa",
                label = h4("Bond angle for monomer A (deg.)"),
                value = 180,
                min = 0,
                max = 180,
                step = 1,
                width = NULL
            ),
            
            numericInput(
                "nb",
                label = h4("Number of monomers B"),
                value = 5,
                min = 1,
                max = 100000,
                step = 1,
                width = NULL
            ),
            numericInput(
                "thetab",
                label = h4("Bond angle for monomer B (deg.)"),
                value = 180,
                min = 0,
                max = 180,
                step = 0.1,
                width = NULL
            ),
            numericInput(
                "nsamples",
                label = h4("Number of polymer molecules to sample"),
                value = 100,
                min = 1,
                max = 10000,
                step = 1,
                width = NULL
            ),
            numericInput(
                "rigidflag",
                label = h4("Stiff bond? (Yes=1, No=0)"),
                value = 0,
                min = 0,
                max = 1,
                step = 1,
                width = NULL
            ),
            actionButton("button1", "Generate Distribution"),
            style = "font-family: 'Times', serif",
            align = "center"
),
        # Show a plot of the generated distribution
        mainPanel(
        fluidRow(    
           plotOutput("distPlot")
        ),  
        fluidRow(
           textOutput("predictedmean")
        ),
        fluidRow(
           textOutput("theoreticalmean")
        ),
          
           style = "font-family: 'Times', serif",
           align = "center",
           size=18
        )

    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    inputList <- eventReactive(input$button1, {
       
        nA = input$na
        nB = input$nb
        tA = input$thetaa
        tB = input$thetab
        n_samples = input$nsamples
        flag = input$rigidflag
      
        return(c(nA,nB,tA,tB,n_samples,flag))
    })
    nA = reactive({
        inputList()[1]
    })
    nB = reactive({
        inputList()[2]
    })
    tA = reactive({
        inputList()[3]
    })
    tB = reactive({
        inputList()[4]
    })
    lengthsq = reactive({
        nA()*(sin(tA()*pi/360)*2)**2+nB()*(sin(tB()*pi/360)*2)**2
    })
    n_samples = reactive({
        inputList()[5]
    })
    iflag = reactive({
        inputList()[6]
    })
    lengths=reactive({
        computedist(nA(),nB(),tA(),tB(),n_samples(),flag=iflag())
    })
    output$distPlot <- renderPlot({
        
        
        #bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        if(!iflag()){
        ggplot(data=as.data.frame(lengths()),aes(x=lengths(),colour="simulation"))+geom_histogram(stat="density")+
        stat_function(fun=function(x) {2*x/lengthsq()*exp(-x^2/lengthsq())}, aes(colour = "theory"))+
        scale_colour_manual("Legend: ", values = c("blue","red"))+
        xlab("polymer length") + ylab("probability density")+
        theme(text=element_text(size=20,  family="Times")) 
        }else{
        ggplot(data=as.data.frame(lengths()),aes(x=lengths(),colour="simulation"))+geom_histogram(stat="density")+
        xlab("polymer length") + ylab("probability density")+
        theme(text=element_text(size=20,  family="Times"))     
        }
           
        
    })
    output$predictedmean <- renderPrint({
        
        c("Predicted mean: ", mean(lengths()))
        
            
        
    }) 
    output$theoreticalmean <- renderPrint({
    if(!iflag()){
    c("Theoretical mean: " , sqrt(pi)/2*sqrt(nA()*(sin(tA()*pi/360)*2)**2+nB()*(sin(tB()*pi/360)*2)**2))
    }else{
        'No available theoretical model'
    }    
})
}
    

# Run the application 
shinyApp(ui = ui, server = server)
