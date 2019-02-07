library(shiny)
library(dplyr)
#library(shinyFiles) #this library is used, but does not need to be loaded since we only call it via shinyFiles::[function]

ui = fluidPage(
  sidebarLayout(
    sidebarPanel(
      #folder input. Creates the input button with some informative texts. 
      p(style = "font-size: 105%; font-weight: bold;", "Select folder:"),
      shinyFiles::shinyDirButton("dir", "Browse...", "Select folder"), #this function is coupled to shinyDirChoose in the server section.
      textOutput("selected.folder", inline = T), #states the currently selected folder
      
      #lists the names of the .bam files found in the x-week-y folder.
      uiOutput("choose_sample", style = "margin-top: 20px;")
    ),
    
    mainPanel(
      tableOutput("table1"),
      #htmlOutput("table2")
      tableOutput("table2")
    )
  )
)

server = function(input, output, session) {
  #########################################################
  # folder input
  #########################################################
  
  #1) create a named list of places where users can start searching
  roots = c(unique(substr(list.dirs(paste0(LETTERS, ":/"), recursive = F), 1, 3)), getwd(), "O:/Transfer/Siebren/17-CHPv4c_week5run3-NORMALStest2")
  names(roots) <- c(roots[1:(length(roots)-2)], getwd(), "demo_location")
  
  #2) use the named list to allow users to select a folder. This function is coupled to shinyDirButton in the ui section.
  shinyFiles::shinyDirChoose(input, 'dir',  updateFreq = 0, roots = roots, defaultRoot = "demo_location")
  
  #3) harvest data from the user's selected folder
  folder.info <- reactiveValues() #all data will be stored in this reactiveValues object.
  observeEvent(input$dir, #priority = 2,
    #input$dir only returns a second variable if a folder is selected
    if(!is.na(input$dir[2])){
      folder.info$path <- shinyFiles::parseDirPath(roots, input$dir) #this location needs to be specified alot of times, so we save it here
      
      folder.info$content = list.files(folder.info$path)
      #collect the various sample names from the x-week-y folder
      week.folder = paste0(folder.info$path, "/", folder.info$content[grep("week", folder.info$content)])
      folder.info$vcf = sub(".vcf", "",
        list.files(week.folder, pattern = ".vcf$")
      )
      folder.info$bam = sub(".bam", "",
        list.files(week.folder, pattern = ".bam$")
      )
      
      #aggregate all unique names in this folder, regardless of file type
      folder.info$samples = unique(
        gsub("-unaligned", "", 
             gsub("\\..*","",
                  list.files(c(week.folder, paste0(folder.info$path, "/originals")))
      )))
      
      #collect the sample names from the /originals folder
      folder.info$unaligned.samples = sub("-unaligned.bam", "",
                                          list.files(paste0(folder.info$path, "/originals"), pattern = ".bam$")
      )
      #check presence of standard files. Returns TRUE or FALSE
      folder.info$QCsample = file.exists(paste0(folder.info$path, "/QCsample.html"))
      folder.info$coverage = file.exists(paste0(folder.info$path, "/Coverage.txt"))
      folder.info$gender   = file.exists(paste0(folder.info$path, "/gender.txt"))
      folder.info$NGSE     = file.exists(paste0(folder.info$path, "/NGSE.html"))
      folder.info$igv      = file.exists(paste0(folder.info$path, "/igv_session.xml"))
    }
  )
  
  
  
  #########################################################
  # Sidebar
  #########################################################
  
  #states the currently selected folder
  output$selected.folder <- renderPrint(
    if(!is.na(input$dir[2]))cat("Selected folder: ", folder.info$path) else cat("No folder selected")
  )
  
  #lists the names of the .bam files found in the x-week-y folder.
  output$choose_sample <- renderUI(if(!is.na(input$dir[2])){
    radioButtons(
      inputId = "sampleselect",   #This does not do anything (yet).
      label = "Samples",
      choices = folder.info$samples
    )
  })
  
  
  
  #########################################################
  # Mainpanel 
  #########################################################
  
  output$table1 <- renderTable(if(!is.na(input$dir[2])){
    tibble(
      "QCsample" = folder.info$QCsample, 
      "Coverage" = folder.info$coverage, 
      "gender" = folder.info$gender, 
      "NGSE" = folder.info$NGSE, 
      "igv" = folder.info$igv
      )
  })
  
  output$table2 <- renderTable(if(!is.na(input$dir[2])){
    t = tibble(
      "sample" = folder.info$samples,
      "unaligned.bam" = NA,
      "bam" = NA,
      "vcf" = NA
    ) %>%
      mutate(
        #unaligned.bam = ifelse(sample %in% folder.info$unaligned.samples, HTML("<p>&#10004;</p>"), HTML("<p>&times;</p>")),
        unaligned.bam = sample %in% folder.info$unaligned.samples,
        bam = sample %in% folder.info$bam,
        vcf = sample %in% folder.info$vcf
      )
  })
  
  
  
  #########################################################
  # Other
  #########################################################
  
  #terminate this app when the browser window is closed.
  session$onSessionEnded(function() {
    stopApp()
  })
  
}

shinyApp(ui = ui, server = server)