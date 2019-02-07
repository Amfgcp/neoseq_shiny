library(shiny)
#library(shinyFiles) #this library is used, but does not need to be loaded since we only call it via shinyFiles::[function]

ui = fluidPage(
  #titlePanel("Hello Shiny!"),

  sidebarLayout(
    sidebarPanel(
      #folder input. Creates the input button with some informative texts. 
      p(style = "font-size: 96%; font-weight: bold;", "Select folder:"),
      shinyFiles::shinyDirButton("dir", "Browse...", "Select folder"), #this function is coupled to shinyDirChoose in the server section.
      textOutput("selected.folder", inline = T) #states the currently selected folder
    ),
    
    mainPanel(
      textOutput("folder")
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
     
      #list all content
      folder.info$content = list.files(folder.info$path)#, pattern = ".vcf$")
      #folder.info$n = length(folder.info$content)
    }
  )
  
  
  
  #########################################################
  # folder output
  #########################################################
  
  #states the currently selected folder
  output$selected.folder <- renderPrint(
    #input$dir only returns a second variable if a folder is selected
    if(!is.na(input$dir[2]))cat("Selected folder: ", folder.info$path) else cat("No folder selected")
  )
  
  #prints the names of all files and folders in the selected folder
  output$folder <- renderText({
    folder.info$content
  })
  
  
  
  #########################################################
  # other
  #########################################################
  
  #terminate this app when the browser window is closed.
  session$onSessionEnded(function() {
    stopApp()
  })
  
}

shinyApp(ui = ui, server = server)