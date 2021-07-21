library(shiny)
library(base64enc)
library(png) 
library(reticulate)
#library(EBImage)
#library(xlsx)
library(ggplot2)
library(readxl)
#use_python("C:/Users/10158/Anaconda3")
# Sys.setenv(RETICULATE_PYTHON = "C:\\Users\\10158\\Anaconda3\\python.exe")
# C:\Users\10158\Anaconda3\envs\CTC_web
#use_python("C:\\Users\\10158\\Anaconda3\\python.exe", required = FALSE)
#use_condaenv(condaenv =  "CTC_web", conda = "auto", required = FALSE)
#C:\Users\10158\Anaconda3\pythonw.exe
#py_config()
#source_python("preprocessing_function.py")
#source_python("equality_preprocess.py")
source_python("CTC_detection_function.py")
#options(shiny.maxRequestSize = 30*1024^10)
data1 <- read_excel("estimate_result.xls")
if (interactive()) {
  ui <- tagList(
    shinythemes::themeSelector(),
    navbarPage(
      theme = "flatly",  # <--- To use a theme, uncomment this
      "CityU CTC Detection Website",
      tabPanel("Cell Segmentation",
               sidebarPanel(
                 titlePanel("Upload Data"),
                 fileInput("upload_BF", "Upload BF image", accept = c('image/png', 'image/jpeg','image/tiff')),
                 fileInput("upload_DAPI", "Upload DAPI image", accept = c('image/png', 'image/jpeg','image/tiff')),
                 fileInput("upload_TRITC", "Upload TRITC image", accept = c('image/png', 'image/jpeg','image/tiff')),
                 actionButton("previwButton", "Preview your images"),
                 actionButton("analysisButton", "Analysis your images")
               ),
               
               mainPanel(
                 tabsetPanel(
                   tabPanel("Original Images", fluidRow(
                     column(10, imageOutput("image1.1", height = 500
                     ))),
                     fluidRow(
                       column(10, imageOutput("image1.2", height = 500
                       ))),
                     fluidRow(
                       column(10, imageOutput("image1.3", height = 500
                       )))),
                   tabPanel("Cell detection result",
                            fluidRow(
                              column(4, imageOutput("image2", height = 300
                              ))))
                   
                   #tabPanel("Tab 3", "This panel is intentionally left blank")
                 )
               )
      ),
      tabPanel("CTC Detection",
               mainPanel(
                 tabsetPanel(
                   tabPanel("CTC Result Image",
                            fluidRow(
                              column(4, imageOutput("image3", height = 500
                              )))),
                   tabPanel("CTC Result Table",
                            fluidPage(
                              titlePanel("Cell Information"),
                              
                              # Create a new Row in the UI for selectInputs
                              fluidRow(
                                column(4,
                                       selectInput("Cell",
                                                   "Statistics:",
                                                   c("All","P_value","Bayes_Factor"
                                                     ))
                                ),
                                column(4,
                                       selectInput("Cell_diameter",
                                                   "Cell diameter:",
                                                   c("All",
                                                     unique(as.character(data1$Cell_diameter))))
                                ),
                                column(4,
                                       selectInput("All",
                                                   "Integrated significance:",
                                                   c("All",
                                                     unique(as.character(data1$All))))
                                )
                              ),
                              # Create a new row for the table.
                              DT::dataTableOutput("table")
                            ),
                            downloadButton('downloadData', 'Download'),) 
                   ))),
      tabPanel("Contact Us", includeMarkdown("About.md"))
    )
  )
  
  server <- function(input, output){
    
    observeEvent(input[["upload_BF"]], {
      inFile <- input[["upload_BF"]]
      if (is.null(inFile))
        return()
      file.copy(inFile$datapath, file.path("C:/Users/10158/OneDrive - City University of Hong Kong/ProgramCode/CTC_Image/02_shiny/CTC/CTC_shiny", inFile$name))
    })
    observeEvent(input[["upload_DAPI"]], {
      inFile <- input[["upload_DAPI"]]
      if (is.null(inFile))
        return()
      file.copy(inFile$datapath, file.path("C:/Users/10158/OneDrive - City University of Hong Kong/ProgramCode/CTC_Image/02_shiny/CTC/CTC_shiny", inFile$name))
    })
    observeEvent(input[["upload_TRITC"]], {
      inFile <- input[["upload_TRITC"]]
      if (is.null(inFile))
        return()
      file.copy(inFile$datapath, file.path("C:/Users/10158/OneDrive - City University of Hong Kong/ProgramCode/CTC_Image/02_shiny/CTC/CTC_shiny", inFile$name))
    })
    num<-'1'
    name<-'1_BF.tif' 
    imageP<-"C:/Users/10158/OneDrive - City University of Hong Kong/ProgramCode/CTC_Image/02_shiny/CTC/CTC_shiny"
    observeEvent(input$previwButton, {
      v1<-original(num,name)
      names(v1)<-c("img","img_nucle","img_yellow")
      output$image1.1 <- renderImage({
        outfile <- tempfile(fileext = ".png")
        #cat(v1$img)
        #browser()
        #cat(table(is.na(v1$img)))
        writePNG(v1$img/255, target = outfile)
        list(src = outfile,
             contentType = "image/png/jpg/jpeg",
             width = 400,
             height = 400,
             alt = "This is alternate text")
      })
      output$image1.2 <- renderImage({
        outfile <- tempfile(fileext = ".png")
        writePNG(v1$img_nucle/255, target = outfile)
        list(src = outfile,
             contentType = "image/png/jpg",
             width = 400,
             height = 400,
             alt = "This is alternate text")
      })
      output$image1.3 <- renderImage({
        outfile <- tempfile(fileext = ".png")
        writePNG(v1$img_yellow/255, target = outfile)
        list(src = outfile,
             contentType = "image/png/jpg",
             width = 400,
             height = 400,
             alt = "This is alternate text")
      })
      
    })
    observeEvent(input$analysisButton, {
      v2<-cell_detection(num,imageP)
      names(v2)<-c("seg_all","All_num","CTC_num","seg_CTC","Added_image")
      
      output$image2 <- renderImage({
        outfile <- tempfile(fileext = ".png")
        writePNG(v2$seg_all, target = outfile)
        # Return a list containing information about the image
        list(src = outfile,
             contentType = "image/png",
             width = 400,
             height = 400,
             alt = "This is alternate text")
      })
      
      output$image3 <- renderImage({
        outfile <- tempfile(fileext = ".png")
        writePNG(v2$seg_CTC/255, target = outfile)
        # Return a list containing information about the image
        list(src = outfile,
             contentType = "image/png",
             width = 400,
             height = 400,
             alt = "This is alternate text")
      })
      output$table <- DT::renderDataTable(DT::datatable({
        data1 <- read_excel("estimate_result.xls")
        if (input$Cell != "All") {
          data1 <- data1[data1$Cell == input$Cell,]
        }
        if (input$Cell_diameter != "All") {
          data1 <- data1[data1$Cell_diameter == input$Cell_diameter,]
        }
        if (input$All != "All") {
          data1 <- data1[data1$All == input$All,]
        }
        data1
      }))
      output$downloadData <- downloadHandler(
        filename = 'CTC_result_table.csv',
        content = function(file) {
          write.csv(data1, file)
        }
      )
    })
    
  }
  shinyApp(ui, server)
}