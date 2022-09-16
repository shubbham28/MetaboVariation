#' MetaboVariation Shiny app launcher
#'
#' Launches the Shiny application for the MetaboVariation functionality. The
#' interactive application allows to analyse data and flag individuals that have variations in their metabolite levels.
#'
#' @export
#'
Shiny.MetaboVariation <- function(){

  ui <- shiny::navbarPage(
    shiny::tags$style( shiny::HTML("#green {
      background-color: #53C1BE;
      }
      .other-sortable .rank-list-container{
      height: 400px;
      }
      ")),
    theme = shinythemes::shinytheme("cerulean"),
    # Application title
    # "Exploring individual variation in metabolites levels",

    # Sidebar with a slider input for number of bins

    shiny::tabPanel("About",shiny::uiOutput("about")),
    shiny::tabPanel("Analysis",shiny::uiOutput("analysis"))
  )

  # Define server logic required to draw a histogram
  server <- function(input, output,session) {
    # Creating data set for analysis
    main = shiny::reactiveValues(data = data.frame())
    # REading file or example dataset to analyse
    shiny::observeEvent(input$data_type,{
      switch (input$data_type,
              "example" = main$data <- MetaboVariation::metabol.data,
              "file" = shiny::observeEvent(input$data_file,{
                main$data = readxl::read_excel(input$data_file$datapath)
              })
      )
    })
    #About tab
    output$about <- renderUI({
      shiny::fluidPage(
        shiny::fluidRow(
          shiny::tags$h1("MetaboVariation"),
          shiny::tags$h4("Metabovariation flags individuals that have significant variation in their
    metabolite levels based on assessment of repeated measures of metabolites. A
    Bayesian generalised linear model (BGLM) is used to flag individuals with
    variation in their metabolite levels and the model can include covariates such
    as sex or age. Individuals are flagged when their metabolite levels lie outside
    of specified posterior predictive intervals."),
    shiny::tags$hr(),
    shiny::tags$hr(),
    shiny::HTML('<center><img src="www/Flowchart_shiny.png" width="720" height="360"></center>')#tags$img(src = "Flowchart_shiny.png", height = 540, width = 1080, align = "center")
        )
      )
    })

    # Analysis tab
    output$analysis <- renderUI({
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::radioButtons("data_type","Type of data",
                       c("Example Dataset" = "example",
                         "Upload your own dataset" = "file"),
                       selected = "example"),
          shiny::fileInput("data_file","Upload the data file"),
          shiny::uiOutput("sortable")
        ),

        # Show a plot of the generated distribution
        shiny::mainPanel(
          shiny::uiOutput("select_metabolites"),
          shiny::uiOutput("results")

        )
      )
    })
    # Actual modelling done here
    model = shiny::reactiveValues(result = list())
    ## Add modeling function here
    shiny::observeEvent(input$start_,{
      main[['data']][,covariates_continous()] = round(main[['data']][,covariates_continous()],2)
      main[['data']][,covariates_categorical()] = factor(main[['data']][,covariates_categorical()])
      main[['data']][,individual_id()] = factor(main[['data']][,individual_id()])
      #model$result = readRDS("~/Data/UCD/Ph.D/Main Research/R codes/result_model.rds")
      #model$result = model$result[final_metabolites()]
      #class(model$result) = c("MetaboVariation","meta.multi_model")
      force(main[['data']][,c(input$individual_id,input$covariates_categorical,input$covariates_continous,input$metabolites)])
      force(individual_id())
      force(c(covariates_continous(),covariates_categorical()))
      force(as.numeric(input$iter))
      force(as.numeric(input$warmup))
      force(as.numeric(input$thinning))
      force(as.numeric(input$cutoff))
      force(final_metabolites())
      model$result = MetaboVariation::MetaboVariation(data = main[['data']][,c(input$individual_id,input$covariates_categorical,input$covariates_continous,input$metabolites)], individual_ids = individual_id(),
                                                      metabolite = final_metabolites(),covariates = c(covariates_continous(),covariates_categorical()),
                                                      iter=as.numeric(input$iter),thin=as.numeric(input$thinning),
                                                      warmup=as.numeric(input$warmup),cutoff=as.numeric(input$cutoff))
    })
    m_single = shiny::reactiveValues(result = list())
    shiny::observeEvent(input$single_metabolite,{
      #m_single$result = readRDS("~/Data/UCD/Ph.D/Main Research/R codes/result_model.rds")
      #m_single$result = m_single$result[[input$single_metabolite]]
      #class(m_single$result) = c("MetaboVariation","meta.single_model")
      m_single$result = model$result[[input$single_metabolite]]
    })

    # Selecting iter, burnin and thinning
    shiny::observeEvent(input$select,{
      output$analysis_start <-shiny::renderUI({
        shiny::fluidPage(shiny::fluidRow(shiny::column(3,shiny::textInput("iter", shiny::HTML("Number of iterations per chain <br/>(Recommended: 3000 - 10000)"), value = "5000", width = '300px', placeholder = "5000")),
                                         shiny::column(3,shiny::textInput("warmup", shiny::HTML("Number of burn-in iterations <br/>(< no. of iterations)"), value = "2500", width = '300px', placeholder = "2500")),
                                         shiny::column(3,shiny::textInput("thinning", shiny::HTML("Thinning rate (k) <br/>(Recommended: 1 - 10)"), value = "2", width = '300px', placeholder = "2")),
                                         shiny::column(3,shiny::textInput("cutoff", shiny::HTML("HPD interval <br/>(Recommended: 0.95)"), value = "0.95", width = '300px', placeholder = "0.95"))),
                  shiny::tags$h6("Number of iterations represents the number of iterations per chain."),
                  shiny::tags$h6("Number of burn-in iterations specifies the number of iterations used for stepsize adaptation of the chains."),
                  shiny::tags$h6("Thinning means storing the k-th sample of the chain to achieve independent sampling."),
                  shiny::tags$h6("Highest Posterior Distribution (HPD) interval is an interval within which an unobserved parameter value falls with a particular probability."),
                  shiny::actionButton("start_","Start Analysis"))
      })
    })
    # Sorting buckets of columns
    output$sortable <- shiny::renderUI({
      if(length(main$data)>0){
        shiny::fluidRow(
          shiny::textInput("divider","Enter the character that seperates the metabolite name from timepoint number in metabolite columns",value = "_"),
          shiny::tags$h5("For column name such as 'MetaboliteA_1', divider is '_'"),
          shiny::tags$b("List of columns"),
          width = 12,
          shiny::actionButton("sort_done","Start when covariates are seperated."),
          sortable::bucket_list(
            header = "Drag the covariates in the desired sections",
            group_name = "bucket_list_group",
            orientation = "vertical",
            sortable::add_rank_list(
              text = "Individual id column here",
              labels = NULL,
              class = c("default-sortable", "other-sortable"),
              input_id = "individual_id"
            ),sortable::add_rank_list(
              text = "Categorical covariates here",
              labels = NULL,
              class = c("default-sortable", "other-sortable"),
              input_id = "covariates_categorical"
            ),
            sortable::add_rank_list(
              text = "Continous covariates here",
              labels = NULL,
              class = c("default-sortable", "other-sortable"),
              input_id = "covariates_continous"
            ),
            sortable::add_rank_list(
              text = "Metabolites columns here",
              labels = colnames(main$data),
              input_id = "metabolites"
            )
          )

        )

      }

    })
    # Data output
    covariates_categorical <- shiny::reactiveVal()
    covariates_continous <- shiny::reactiveVal()
    individual_id <- shiny::reactiveVal()
    metabolites <- shiny::reactiveVal()
    shiny::observeEvent(input$sort_done,{
      covariates_categorical(input$covariates_categorical)
      covariates_continous(input$covariates_continous)
      individual_id(input$individual_id)
      metabolites(MetaboVariation::get.metabolites(input$metabolites,divider = input$divider))
    },once = TRUE)

    output$table <- shiny::renderTable({
      if(!identical(input$covariates_categorical, character(0)) | !identical(input$covariates_continous, character(0)) | !identical(input$individual_id, character(0))){
        head(main$data[,c(input$individual_id,input$covariates_categorical,input$covariates_continous)])
      }
      else{
        head(main$data[,c(individual_id(),covariates_categorical(),covariates_continous())])
      }
    })
    # Final metabolites selected for analysis
    final_metabolites = shiny::eventReactive(input$select,{input$view_metabolites})
    shiny::observeEvent(input$select,{
      output$results_metabolites <-shiny::renderUI({
        shiny::HTML(paste(paste("<br/>These following metabolites will be analysed"),paste(final_metabolites(),collapse = ", "),sep='<br/><br/>'))
      })
    })

    # Analysis start section
    output$select_metabolites <- shiny::renderUI({
      if(length(metabolites()) >1){
        shiny::fluidRow(            shiny::tags$h2("Dataset"),
                             shiny::tableOutput("table"),
                             shiny::downloadButton('download_initial', 'Download dataset'),
                             shiny::tags$br(),
                             shiny::tags$br(),
                             shiny::selectizeInput("view_metabolites", "Select metabolites to analyse (>1).", choices  = metabolites(),
                                            selected = metabolites()[1], multiple = TRUE,
                                            options = list(plugins = list("remove_button",  'drag_drop'))),
                             shiny::tags$h6("The analysis will be performed on the selected metabolites:"),
                             shiny::actionButton("select","Select"),
                             shiny::htmlOutput("results_metabolites"),
                             shiny::tags$br(),
                             shiny::uiOutput("analysis_start")
        )

      }

    })
    # Results section
    shiny::observeEvent(input$start_,{
      output$results <- shiny::renderUI({
        shiny::fluidPage(
          shiny::titlePanel("Analysis for selected individuals"),
          shiny::mainPanel(
            #          verbatimTextOutput("data_type"),
            shiny::tabsetPanel(type = "tabs",
                               shiny::tabPanel("Individual Metabolite",shiny::uiOutput("indimetabolite")),
                               shiny::tabPanel("Summary",shiny::uiOutput("summary")))
          )
        )
      })
    })
    #  output$data_type = renderPrint({input$individual_id})
    # Summary tab of results
    output$summary = shiny::renderUI({
      shiny::fluidRow(
        shiny::tags$hr(),
        shiny::tags$h3("Number of times individuals have been flagged in each metabolite"),
        shiny::dataTableOutput("multiple_table"),
        shiny::downloadButton('download_group_data', 'Download results'),
        shiny::tags$br(),
        shiny::tags$h3("Number of individuals detected in each timepoint for every metabolites"),
        plotly::plotlyOutput("plot_metabolite_count_multiple"),
        shiny::tags$hr()

      )
    })
    # Individual Metabolite tab of results
    output$indimetabolite = shiny::renderUI({
      shiny::fluidRow(
        shiny::selectizeInput("single_metabolite", "Select metabolite to show", choices  = final_metabolites(),
                       selected = final_metabolites()[1]),
        shiny::tags$h3("Circos Plot"),
        shiny::plotOutput("circos",width = 720,height = 540),
        shiny::tags$hr(),
        shiny::tags$h3("Number of individuals flagged in each timepoint"),
        plotly::plotlyOutput("plot_metabolite_count_single"),
        shiny::tags$hr(),
        shiny::tags$h3("Flagged Individuals"),
        shiny::dataTableOutput("single_table"),
        shiny::downloadButton('download_individual_data', 'Download results'),
        shiny::tags$br(),
        shiny::tags$h4("Model summary"),
        shiny::htmlOutput("model_details")

      )
    })
    # Plots
    output$circos <- shiny::renderPlot(plot(m_single$result,type = "circos"))
    output$plot_metabolite_count_single <- plotly::renderPlotly(plot(m_single$result,type = "metabolites_count"))
    output$plot_metabolite_count_multiple <- plotly::renderPlotly(plot(model$result,type = "metabolites_count"))
    output$single_table <- shiny::renderDataTable(MetaboVariation::flagged.Individuals(m_single$result,data = main[['data']][,c(input$individual_id,input$covariates_categorical,input$covariates_continous,input$metabolites)], individual_id = individual_id()))
    output$multiple_table <- shiny::renderDataTable(MetaboVariation::flagged.Individuals(model$result,data = main[['data']][,c(input$individual_id,input$covariates_categorical,input$covariates_continous,input$metabolites)], individual_id = individual_id()))
    output$model_details <- shiny::renderUI({
      covar<-paste0("<br/>Significant covariates for the metabolite: ",paste(rownames(m_single$result$significant_covariates[which(m_single$result$significant_covariates[,5]==1),]),collapse = ", "))
      iter<-paste0("<br/>Number of iterations per chain: ",m_single$result$iter)
      warmup<-paste0("<br/>Number of burn-in(warmup) iterations per chain: ",m_single$result$warmup)
      thinning<-paste0("<br/>Thinning rate: ",input$thinning)
      cutoff<-paste0("<br/>HPD interval probability of the model: ",input$cutoff)
      rhat<-paste0("<br/>Convergence of chains(Rhat): ",round(m_single$result$Rhat,2))
      note<-paste("<br/>Rhat value should lies between 0.9 and 1.05. If not, increase the number of iterations")
      shiny::HTML(paste(covar,iter,warmup,thinning,cutoff,rhat,note,sep = '<br/>'))
    })

    # Download handlers
    output$download_initial <- shiny::downloadHandler(
      filename = function() {
        paste('data_', format(Sys.time()), '.xlsx', sep='')
      },
      content = function(file) {
        writexl::write_xlsx(main$data[,c(input$individual_id,input$covariates_categorical,input$covariates_continous,input$metabolites)], file)
      },
      contentType = ".xlsx"
    )
    output$download_group_data <- shiny::downloadHandler(
      filename = function() {
        paste('results_', format(Sys.time()), '.xlsx', sep='')
      },
      content = function(file) {
        writexl::write_xlsx(MetaboVariation::flagged.Individuals(model$result,data = main$data, individual_id = individual_id()), file)
      },
      contentType = ".xlsx"
    )
    output$download_individual_data <- shiny::downloadHandler(
      filename = function() {
        paste('results_', format(Sys.time()), '.xlsx', sep='')
      },
      content = function(file) {
        writexl::write_xlsx(MetaboVariation::flagged.Individuals(m_single$result,data = main$data, individual_id = individual_id()), file)
      },
      contentType = ".xlsx"
    )


    # session$allowReconnect(TRUE)
  }

  # Run the application
  shinyApp(ui = ui, server = server)

}

.onAttach <- function(libname, pkgname) {
  shiny::addResourcePath('www',
                         system.file('www',
                                     package = 'MetaboVariation'))
}
