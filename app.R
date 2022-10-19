library(shiny)
library(readxl)
if(!require(grdR)) {
  devtools::install_github("mauromiguelm/grdR")
}

# Define UI for application
ui <- 
  fluidPage(
    tags$head(
      tags$link(rel = "stylesheet", href = "styles.css")
    ),
  navbarPage("grdR",
             tags$link(rel = "stylesheet", type="text/css", href="www/style.css"),
             tabPanel("Introduction",
                      
                      fluidRow(h2("Introduction to grdR")),
                      fluidRow(p("Drug sensitivity metrics are important tools that aid basic research both drug discovery and pharmacology, but also clinical drug therapy assessment. The widely used metrics IC50, GI50 and GR50 incorporate the effect of concentration at endpoint, including or not corrections for baseline (treatment naïve) and growth rate. These metrics assume a monotonic time response, meaning the assumption is that the growth effect occurs at an early time point after drug treatment. We show that this assumption does not hold true in many cases, leading to biases in the evaluation of drug sensitivity assessment. We solved this problem by developing the metric growth response delay (GRD), an estimation of a drug’s growth effect lag time, together with a model that can simulate both growth response delay and drug potency effects. The growth response delay metric can be used to accurately estimate and study time-dependent drug effects, and to correct the conventional drug sensitivity metrics by the growth effect lag time. To demonstrate, we show a large growth response delay diversity both within and between 31 anticancer drugs, and that basal enzyme levels involved in asparaginase drug sensitivity response have a strong growth response delay pattern.")),
                      sidebarLayout(
                        sidebarPanel(
                          sliderInput("GRD", label = h3("GRD"), min = 0, 
                                      max = 3, value = 1.6,step = 0.1),
                          sliderInput("halfeff", label = h3("Drug potency"), min = 0.1, 
                                      max = 50, value = 2, step = 0.1),
                          sliderInput("k", label = h3("Growth rate"), min = 0, 
                                      max = 10, value = 0.6, step = 0.1),
                          sliderInput("tz", label = h3("Initial seed"), min = 1, 
                                      max = 1000, value = 5, step = 10),
                        ),
                        mainPanel(
                          plotOutput("plotExample")
                        )
                      )
               
               
             ),
             
             tabPanel("Calculator",
    
    #import data

    fluidRow(
      
      column(4,
             fileInput(inputId = "inputData", 
                       label = "Choose Input Data", 
                       multiple = FALSE, 
                       accept = ".csv", 
                       placeholder = "Input data in GRD format")
             ),
      column(1,h3("OR")),
      
      column(3,
             br(),
             actionButton(inputId = "giveSampleData", 
                          label="Use sample data as input", 
                          icon = NULL),
             
             ),
      column(1,
             br(),
             actionButton(inputId = "clearData", 
                          label="Clear input data", 
                          icon = NULL)
             
      )
             
             
             ),
    
  
      fluidRow(column(4,selectInput("columns", "Select Columns", choices = NULL,multiple = TRUE)),
               column(3,p("Columns that define treatment groups. For sample data, only agent should be selected."))
      ),
    
    fluidRow(tableOutput("checkFunction")),
    
    conditionalPanel(
      condition = "output.existsData",
    
    fluidRow(actionButton(inputId = "crunchData", 
                          label="Calculate GRD", 
                          icon = NULL)),
    fluidRow(tableOutput("GRD")),
    
    )),
    tabPanel("Figures",
             h3("Select data in order to display figures."),
             conditionalPanel(
               condition = "output.existsData",
             fluidRow(plotOutput("plot_GRD"),width = "50%"),{
               fluidRow(column(4,selectInput("chooseGroup", "Select Columns", choices = NULL,multiple = FALSE)),
                        column(6,plotOutput("plot_growth")))
               
             }
             
             
             ))
             
    
    )
  )
             
             
    
# Define server logic required to draw a histogram
server <- function(input, output) {
  
  
  # validate input data 
  
  storeData <- reactiveValues(data = NULL)
  
  storeData <- reactiveValues(GRD = NULL)
  
  observeEvent(input$inputData,{
    file <- input$inputData
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "csv", "Please upload a csv file"))
    
    storeData$data = read.csv(file$datapath,check.names = FALSE)
    
  })
  
  observeEvent(input$giveSampleData,{
    storeData$data = grdR::sample_data
  })
  
  observeEvent(input$clearData,{
    storeData$data = NULL
    output$data = NULL
    output$plot = NULL
    output$GRD = NULL
    output$plot_GRD = NULL
    output$plot_growth = NULL
  })
  

  output$checkFunction <- renderTable({
    if(is.null(storeData$data)){return()}
    
    vars = colnames(storeData$data)
    
    groups = unique(storeData$data$perturbation)
    
    updateSelectInput(inputId =  "columns",label = "Select Columns", choices = vars)
    
    updateSelectInput(inputId =  "chooseGroup",label = "Select group to plot", choices = groups)
  
    return(head(storeData$data))
  })
  
  
  observeEvent(input$crunchData,{
    
    storeData$GRD <- get_fit(storeData$data, groupingVariables = input$columns,
                          smoothData = FALSE, upperLimitThreshold = 1,
                          timeTreatment = 1,upperLimit = 1.0,orderConc = TRUE,saveModel = FALSE)
    
    insertUI("#crunchData",where = "afterEnd",
             ui = h3("Results:"))
    
  })
  
  
  output[["existsData"]] <- reactive({
    class(storeData$data) == "data.frame"

  })

  outputOptions(output, "existsData", suspendWhenHidden = FALSE)
  
  
  
  output$GRD <- renderTable({
    if(is.null(storeData$GRD)){return()}
    
    return(storeData$GRD)
    
  })
  
  output$plot_GRD <- renderPlot({
    # hist(x = storeData$GRD$grd)
    
    my_bar <- grdR::plot_GRD(storeData$GRD, groupingVariables = input$columns)
    
  })
  
  finalInput <- reactive({
    list(subset(storeData$data, perturbation == input$chooseGroup),
         subset(storeData$GRD, perturbation == input$chooseGroup,"GRD"))
  })
  
  
  output$plot_growth <- renderPlot({
    # hist(x = storeData$GRD$grd)
    
    data = finalInput()
    
    my_curve <- plot_growth(inputData =data[[1]], groupingVariables = "perturbation", GRD = data[[2]])
    
    return(my_curve)
    
    
  })
  
  exampleInput <- reactive({
    nreps <- 1
    time <- seq(0,3,0.1)
    
    conc <- matrix(base::rep(c(0,seq(0.5,8,length.out = 6)),  each = nreps), dimnames = list(NULL, "concentration"))
    
    params <- matrix(data = c(input$halfeff, input$GRD, input$k), dimnames = list(NULL, c("halfeff", "t_onset","k")), nrow = nreps, ncol = 3)
    params <- cbind(do.call("rbind", rep(list(params), 7)), conc)
    params <- do.call(rbind, replicate(length(time), params, simplify=FALSE)) #replicate matrix to acommodate time
    params <- cbind(params, time = rep(time, times = nreps))
    
    #generate simulated data
    
    sapply(1:dim(params)[1], function(x)
      
      model_growth(cell_tz = input$tz, t = params[x,'time'], t_onset = params[x,"t_onset"], k =  params[x,"k"], maxeff = 1.0, halfeff = params[x,"halfeff"], conc = params[x,"concentration"], hill = 1.6)
      
    ) -> sample_data
    
    sample_data <- cbind(params, output = sample_data)
    
    sample_data <- reshape2::dcast(data.frame(sample_data),halfeff + k +t_onset +concentration ~time, value.var="output")
    
    # prepare data for input in GRDR
    
    tidyr::gather(data.frame(sample_data, check.names = F), key = "time", value = "cell_count", -c("halfeff","k", "t_onset", "concentration")) -> sample_data
    
    sample_data <- cbind(data.frame(perturbation = paste(sample_data$halfeff,sample_data$k, sample_data$t_onset, sep = "_")), sample_data)
    
    sample_data$perturbation <- factor(sample_data$perturbation,labels = LETTERS[1:length(unique(sample_data$perturbation))])
    
    sample_data$time <- as.numeric(sample_data$time)
    
    sample_data <- subset(sample_data,perturbation == "A")
    
  })
  
  
  output$plotExample <- renderPlot({
    grdR::plot_growth(inputData = exampleInput(),groupingVariables = "perturbation",GRD = NULL,addLegend = TRUE)
  })
  
  
  }

# Run the application 
shinyApp(ui = ui, server = server)




