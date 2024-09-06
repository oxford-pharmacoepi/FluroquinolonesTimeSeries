library(shinydashboard)
library(dplyr)
library(readr)
library(here)
library(stringr)
library(DT)
library(shinycssloaders)
library(shinyWidgets)
library(gt)
library(scales)
library(kableExtra)
library(tidyr)
library(stringr)
library(ggplot2)
library(fresh)
library(plotly)
library(IncidencePrevalence)
# library(DiagrammeR)
# library(DiagrammeRsvg)
# library(lubridate)
# library(orca)
load("mergeData.Rdata")

source(here("functions.R"))



# theme -----
mytheme <- create_theme(
  adminlte_color(
    light_blue = "#62AFD6"
  ),
  adminlte_sidebar(
    # width = "400px",
    dark_bg = "#9AD1ED",
    dark_hover_bg = "#3B9AB2",
    dark_color = "black"
  ), 
  adminlte_global(
    content_bg = "#eaebea" 
  ),
  adminlte_vars(
    border_color = "#112446",
    active_link_hover_bg = "#FFF",
    active_link_hover_color = "#112446",
    active_link_hover_border_color = "#112446",
    link_hover_border_color = "#112446"
  )
)

# ui -----
ui <- dashboardPage(
  dashboardHeader(title = "Fluroquinolones studyathon - interrupted time series modelling"),
  ## menu ----
  dashboardSidebar(
    sidebarMenu(
      menuItem(
        text = "ARIMA",
        tabName = "arima_tab",
        menuSubItem(
          text = "ARIMA estimations",
          tabName = "arima_pred_tab"
        ),
        menuSubItem(
          text = "ARIMA parameters",
          tabName = "arima_para"
        ),
        menuSubItem(
          text = "ARIMA coefficients",
          tabName = "arima_coef"
        ),
        menuSubItem(
          text = "ARIMA diagnostics",
          tabName = "arima_diag"
        )
      ),
      menuItem(
        text = "Segemented regression",
        tabName = "segreg_tab",
        menuSubItem(
          text = "Segmented regression estimations",
          tabName = "segreg_pred"
        ),
        menuSubItem(
          text = "Segmented regression coefficients",
          tabName = "segreg_coef"
        ),
        menuSubItem(
          text = "Segmented regression diagnostics",
          tabName = "segreg_diag"
        )
      )
    )
  ),
  
  ## body ----
  dashboardBody(
    use_theme(mytheme),
    tabItems(
      # prediction ----
      tabItem(
        tabName = "arima_pred_tab",
        h4("Interrupted time series using ARIMA"),
        selectors(arima_pred, "pred", c("Model", "cdm_name", "age_group"), multiple = TRUE),
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "Raw table",
            downloadButton("arima_raw_download", "Download csv"),
            DT::dataTableOutput("arima_pred") %>% withSpinner()
          ),
          tabPanel(
            "Plot",
            downloadButton("arima_plot_download_png", "Download png"),
            plotlyOutput("arima_plot") %>% withSpinner()
          )
        )
      ),
      # diagnostics ---
      tabItem(
        tabName = "arima_para",
        h4("Parameters of ARIMA models"),
        selectors(arima_para, "arima_para", c("cdm_name", "age_group"), multiple = FALSE),
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "Raw table",
            downloadButton("diagnostics_raw_download", "Download csv"),
            DT::dataTableOutput("arima_para_table") %>% withSpinner()
          )
        )
      ),
      #--ARIMA coefficients
      tabItem(
        tabName = "arima_coef",
        h4("Coefficients of ARIMA models"),
        #selectors(newuserinc, "newuserinc", c("cdm_name", "outcome_cohort_name", "denominator_age_group", "denominator_sex"), multiple = TRUE),
        #selectors(newuserinc, "newuserinc", c("analysis_interval"), multiple = FALSE),
        selectors(arima_coef, "arima_coef", c("Model", "cdm_name", "age_group"), multiple = FALSE),
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "Raw table",
            downloadButton("coef_raw_download", "Download csv"),
            DT::dataTableOutput("arima_coef_table") %>% withSpinner()
          )
        )
      ),
      tabItem(
        tabName = "segreg_pred",
        h4("Interrupted time series using segmented regression"),
        selectors(segreg_pred, "pred_segreg", c("Model", "cdm_name", "age_group"), multiple = TRUE),
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "Raw table",
            downloadButton("segreg_raw_download", "Download csv"),
            DT::dataTableOutput("segreg_pred_table") %>% withSpinner()
          ),
          tabPanel(
            "Plot",
            downloadButton("segreg_plot_download_png", "Download png"),
            plotlyOutput("segreg_plot") %>% withSpinner()
          )
        )
      ),
      #--Segmented regression coefficients
      tabItem(
        tabName = "segreg_coef",
        h4("Coefficients of segmented regression models"),
        #selectors(newuserinc, "newuserinc", c("cdm_name", "outcome_cohort_name", "denominator_age_group", "denominator_sex"), multiple = TRUE),
        #selectors(newuserinc, "newuserinc", c("analysis_interval"), multiple = FALSE),
        selectors(segreg_coef, "segreg_coef", c("Model", "cdm_name", "age_group"), multiple = FALSE),
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "Raw table",
            downloadButton("coef_raw_download_segreg", "Download csv"),
            DT::dataTableOutput("segreg_coef_table") %>% withSpinner()
          )
        )
      ),
      
      # New tab item for Residual Diagnostics
      tabItem(
        tabName = "arima_diag",
        h4("ACF PACF Plots for ARIMA model residuals"),
        tabPanel(title = "ACF PACF Plots",
                 fluidRow(
                   column(4,
                          pickerInput(
                            inputId = "cdm_name_arima",
                            label = "cdm_name:",
                            choices = unique(arima_pred$cdm_name),  
                            selected = "CPRD GOLD",
                            multiple = FALSE
                          )
                   ),
                   column(4,
                          pickerInput(
                            inputId = "age_group_arima",
                            label = "age_group:",
                            choices = unique(arima_pred$age_group),  
                            selected = "19 - 59",
                            multiple = FALSE
                          )
                   )
                 ),
                 imageOutput("acf_pacf_plot_arima"))),
      # New tab item for Residual Diagnostics
      tabItem(
        tabName = "segreg_diag",
        h4("ACF PACF Plots for segmented regression model residuals"),
        tabPanel(title = "ACF PACF Plots",
                 fluidRow(
                   column(4,
                          pickerInput(
                            inputId = "cdm_name_segreg",
                            label = "cdm_name:",
                            choices = unique(segreg_pred$cdm_name),  
                            selected = "CPRD GOLD",
                            multiple = FALSE
                          )
                   ),
                   column(4,
                          pickerInput(
                            inputId = "age_group_segreg",
                            label = "age_group:",
                            choices = unique(segreg_pred$age_group),  
                            selected = "19 - 59",
                            multiple = FALSE
                          )
                   )
                 ),
                 imageOutput("acf_pacf_plot_segreg")))
      
    ),
    
    # end ----
  ),
  
  tags$head(
    tags$style(HTML("
            .skin-blue .sidebar-menu .treeview-menu li a {
                color: black !important;
            }
            .skin-blue .main-header .logo {
                color: black !important;
            }
        ")))
)


server <- function(input, output, session) {
  getForecastData_segreg <- reactive({
    filterData(segreg_pred, "pred_segreg", input)
  })
  
  output$segreg_pred_table <- renderDataTable({
    datatable(getForecastData_segreg(), options = list(scrollX = TRUE))
  })
  
  
  output$segreg_plot <- renderPlotly({
    prev <- getForecastData_segreg()
    prev <- prev %>% mutate(group = paste(Model,
                                                   cdm_name,
                                                   age_group,
                                                   sep = " "),
                            group_ori = paste(age_group,
                                                       cdm_name,
                                                       sep = " ")) %>%
      dplyr::arrange(Date) %>% mutate(Date = as.Date(Date))
    
    p <- ggplot(prev, aes(x = Date, y = Forecast, color = group)) +
      ylab("Incidence Rate") +
      geom_line(linetype = "dotted") +
      geom_line(aes(y = incidence_100000_pys, color = group_ori)) +
      geom_point(aes(y = incidence_100000_pys, color = group_ori)) +
      geom_errorbar(aes(y = incidence_100000_pys, 
                        ymin = incidence_100000_pys_95CI_lower, 
                        ymax = incidence_100000_pys_95CI_upper, 
                        color = group_ori), width = 0.2) +
      geom_point(aes(y = Forecast, color = group)) +
      geom_errorbar(aes(y = Forecast, 
                        ymin = pmax(0, Lower_95_CI, na.rm = TRUE), 
                        ymax = pmax(0, Upper_95_CI, na.rm = TRUE), 
                        color = group), width = 0.2) +
      geom_vline(xintercept = as.numeric(as.Date("2019-04-01")), linetype="dashed", color = "red") +
      theme_minimal()+
      theme(text = element_text(size = 12),  
            axis.title = element_text(size = 14),  
            axis.text = element_text(size = 12),  
            legend.title = element_text(size = 14),  
            legend.text = element_text(size = 12)) 
    
    
    
    p_plotly <- ggplotly(p) %>%
      layout(
        legend = list(
          orientation = "h", 
          x = 0.5, 
          y = -0.3,  
          xanchor = "center", 
          yanchor = "top"
        ),
        margin = list(
          b = 100  
        )
      )
    
    
    p_plotly
  })
  
  output$segreg_plot_download_png <- downloadHandler(
    filename = "segmented_regression_plot.png",
    content = function(file) {
      prev <- getForecastData_segreg()
      prev <- prev %>% mutate(group = paste(Model,
                                                     cdm_name,
                                                     age_group,
                                                     sep = " "),
                              group_ori = paste(age_group,
                                                         cdm_name,
                                                         sep = " ")) %>%
        dplyr::arrange(Date) %>% mutate(Date = as.Date(Date))
      
      p <- ggplot(prev, aes(x = Date, y = Forecast, color = group)) +
        ylab("Incidence Rate") +
        geom_line(linetype = "dotted") +
        geom_line(aes(y = incidence_100000_pys, color = group_ori)) +
        geom_point(aes(y = incidence_100000_pys, color = group_ori)) +
        geom_errorbar(aes(y = incidence_100000_pys, 
                          ymin = incidence_100000_pys_95CI_lower, 
                          ymax = incidence_100000_pys_95CI_upper, 
                          color = group_ori), width = 0.2) +
        geom_point(aes(y = Forecast, color = group)) +
        geom_errorbar(aes(y = Forecast, 
                          ymin = pmax(0, Lower_95_CI, na.rm = TRUE), 
                          ymax = pmax(0, Upper_95_CI, na.rm = TRUE), 
                          color = group), width = 0.2) +
        geom_vline(xintercept = as.numeric(as.Date("2019-04-01")), linetype="dashed", color = "red") +
        theme_minimal()
      # Save the plotly object as an HTML file temporarily
      temp_html <- tempfile(fileext = ".html")
      
      plotly_obj <- ggplotly(p) # Convert to plotly object

      saveWidget(plotly_obj, temp_html, selfcontained = TRUE)
      
      # Use webshot to capture this HTML file as a PNG image
      webshot::webshot(temp_html, file = file, delay = 1) 
    },
    contentType = "image/png"
  )

  
  getForecastData <- reactive({
    filterData(arima_pred, "pred", input)
  })
  
  output$arima_pred <- renderDataTable({
    datatable(getForecastData(), options = list(scrollX = TRUE))
  })
  
  output$arima_plot <- renderPlotly({
    prev <- getForecastData()
    prev <- prev %>% mutate(group = paste(Model,
                                                   cdm_name,
                                                   age_group,
                                                   sep = " "),
                            group_ori = paste(age_group,
                                                       cdm_name,
                                                       sep = " ")) %>%
      dplyr::arrange(Date) %>% mutate(Date = as.Date(Date))
    
    p <- ggplot(prev, aes(x = Date, y = Forecast, color = group)) +
      ylab("Incidence Rate") +
      geom_line(linetype = "dotted") +
      geom_line(aes(y = incidence_100000_pys, color = group_ori)) +
      geom_point(aes(y = incidence_100000_pys, color = group_ori)) +
      geom_errorbar(aes(y = incidence_100000_pys, 
                        ymin = incidence_100000_pys_95CI_lower, 
                        ymax = incidence_100000_pys_95CI_upper, 
                        color = group_ori), width = 0.2) +
      geom_point(aes(y = Forecast, color = group)) +
      geom_errorbar(aes(y = Forecast, 
                        ymin = pmax(0, Lower_95_CI, na.rm = TRUE), 
                        ymax = pmax(0, Upper_95_CI, na.rm = TRUE), 
                        color = group), width = 0.2) +
      geom_vline(xintercept = as.numeric(as.Date("2019-04-01")), linetype="dashed", color = "red") +
      theme_minimal()+
      theme(text = element_text(size = 12),  
            axis.title = element_text(size = 14),  
            axis.text = element_text(size = 12),  
            legend.title = element_text(size = 14),  
            legend.text = element_text(size = 12)) 

    

    p_plotly <- ggplotly(p) %>%
      layout(
        legend = list(
          orientation = "h", 
          x = 0.5, 
          y = -0.3,  
          xanchor = "center", 
          yanchor = "top"
        ),
        margin = list(
          b = 100  
        )
      )
    
    p_plotly
    # Return the plot

    
  })

  output$arima_plot_download_png <- downloadHandler(
    filename = function() {
      paste("arima-plot-", Sys.Date(), ".png", sep="")
    },
    content = function(file) {
      prev <- getForecastData()
      prev <- prev %>% mutate(group = paste(Model,
                                                     cdm_name,
                                                     age_group,
                                                     sep = " "),
                              group_ori = paste(age_group,
                                                         cdm_name,
                                                         sep = " ")) %>%
        dplyr::arrange(Date) %>% mutate(Date = as.Date(Date))
      
      p <- ggplot(prev, aes(x = Date, y = Forecast, color = group)) +
        ylab("Incidence Rate") +
        geom_line(linetype = "dotted") +
        geom_line(aes(y = incidence_100000_pys, color = group_ori)) +
        geom_point(aes(y = incidence_100000_pys, color = group_ori)) +
        geom_errorbar(aes(y = incidence_100000_pys, 
                          ymin = incidence_100000_pys_95CI_lower, 
                          ymax = incidence_100000_pys_95CI_upper, 
                          color = group_ori), width = 0.2) +
        geom_point(aes(y = Forecast, color = group)) +
        geom_errorbar(aes(y = Forecast, 
                          ymin = pmax(0, Lower_95_CI, na.rm = TRUE), 
                          ymax = pmax(0, Upper_95_CI, na.rm = TRUE), 
                          color = group), width = 0.2) +
        geom_vline(xintercept = as.numeric(as.Date("2019-04-01")), linetype="dashed", color = "red") +
        theme_minimal()+
        theme(text = element_text(size = 12),  
              axis.title = element_text(size = 14),  
              axis.text = element_text(size = 12),  
              legend.title = element_text(size = 14),  
              legend.text = element_text(size = 12))  

      
      
      plotly_obj <- ggplotly(p)
      
      plotly_obj <- plotly_obj %>% layout(legend = list(y = -0.2, yanchor = "top")) 
      
      # Save the plotly object as an HTML file temporarily
      temp_html <- tempfile(fileext = ".html")
      saveWidget(plotly_obj, temp_html, selfcontained = TRUE)
      
      # Use webshot to capture this HTML file as a PNG image
      webshot::webshot(temp_html, file = file, delay = 1) 
    },
    contentType = "image/png"
  )
  getDiagData <- reactive({
    filterData(arima_para, "arima_para", input)
  })
  
  output$arima_para_table <- renderDataTable({
    datatable(getDiagData(), options = list(scrollX = TRUE))
  })
  
  getCoefData <- reactive({
    filterData(arima_coef, "arima_coef", input)
  })
  
  output$arima_coef_table <- renderDataTable({
    datatable(getCoefData(), options = list(scrollX = TRUE))
  })
  

  getCoefData_seg <- reactive({
    filterData(segreg_coef, "segreg_coef", input)
  })
  
  output$segreg_coef_table <- renderDataTable({
    datatable(getCoefData_seg(), options = list(scrollX = TRUE))
  })
  
  # Other reactive expressions and output renderings...
  
  # Correctly define reactive expressions for image paths
  
  
  
  reactiveARIMAPlotPath <- reactive({

    
    file_name <- if(input$age_group_arima == "19 - 59") {
        paste0("ARIMA_diag_", input$cdm_name_arima, "_young.png")
    } else{
      paste0("ARIMA_diag_", input$cdm_name_arima, "_old.png")
      }
        
    file_path <- here::here("www_gold/arima/", file_name)
    file_path # This line returns the path from the reactive expression
  })
  
  reactiveSegregPlotPath <- reactive({
    file_name <- if(input$age_group_segreg == "19 - 59") {
      paste0("SegReg_diag_", input$cdm_name_segreg, "_young.png")
    } else{
      paste0("SegReg_diag_", input$cdm_name_segreg, "_old.png")
      
    }
    file_path <- here::here("www_gold/segreg/", file_name)
    file_path # This line returns the path from the reactive expression
  })
  
  
  # Use reactive expressions in renderImage for scatterPlot
  output$acf_pacf_plot_arima <- renderImage({
    src <- reactiveARIMAPlotPath()
    style <- paste("display: block; margin: 0 auto;")
    list(src = src,
         width = "auto",
         height = "auto",
         style = style)
  }, deleteFile = FALSE)
  
  # Use reactive expressions in renderImage for qqPlot
  output$acf_pacf_plot_segreg <- renderImage({
    src <- reactiveSegregPlotPath()
    style <- paste("display: block; margin: 0 auto;")
    list(src = src,
         contentType = 'image/png',
         width = "auto",
         height = "auto",
         style = style)
  }, deleteFile = FALSE)
  
  # end ----
}

shinyApp(ui = ui, server = server)