library(shiny)
library(dplyr)
library(xtable) #to make html work in the tables. Doesn't work with xtable::
#library(shinyFiles) #this library is used, but does not need to be loaded since we only call it via shinyFiles::[function]
#library(stringdist) #for ClosestMatch2's amatch

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
      uiOutput("title", style = "margin-left: 30px; margin-top: 15px; margin-bottom: 20px;"),
      tabsetPanel(
        tabPanel(
          "Overview",
          tableOutput("table1")#,
          #tableOutput("table2")
        ),
        tabPanel(
          ".bam files"
        ),
        tabPanel(
          ".vcf files"
        ),
        tabPanel(
          "RNAseq"
        ),
        tabPanel(
          "HLA plots"
        )
      )
    )
  )
)

server = function(input, output, session) {
  #########################################################
  # folder input
  #########################################################
  
  #1) create a named list of places where users can start searching
  roots = c(unique(substr(list.dirs(paste0(LETTERS, ":/"), recursive = F), 1, 3)), getwd(), "O:/Transfer/Siebren/neoSeq")
  names(roots) <- c(roots[1:(length(roots)-2)], getwd(), "demo_location")
  
  #2) use the named list to allow users to select a folder. This function is coupled to shinyDirButton in the ui section.
  shinyFiles::shinyDirChoose(input, 'dir',  updateFreq = 0, roots = roots, defaultRoot = "demo_location")
  go = reactive(!is.na(input$dir[2])) #True when a folder has been selected
  
  #3) harvest data from the user's selected folder
  folder.info <- reactiveValues()
  
  observeEvent(input$dir, if(go()){
    folder.info$path <- shinyFiles::parseDirPath(roots, input$dir)
    
    #list of samples
    TNsamples = sub("_L..*", "",
      list.files(paste0(folder.info$path, "/fastq"), pattern = ".fastq.gz$")
    )
    folder.info$samples = unique(substr(TNsamples,1,nchar(TNsamples)-1))
    
    #list of bam samples
    folder.info$bam = tibble(sample = folder.info$samples) %>%
      mutate(
       normal = file.exists(paste0(folder.info$path, "/bam/", sample, "N.dedup.recal.bam")),
       tumor = file.exists(paste0(folder.info$path, "/bam/", sample, "T.dedup.recal.bam")),
       both = normal & tumor
      )
    
    #list of vcf samples
    folder.info$vcf = sub("_CombineVariants.vcf", "",
      list.files(paste0(folder.info$path, "/vcf"), pattern = ".vcf$")
    )
    
    #list of RNAseq samples
    intersect2 <- function (x, y){
      y <- as.vector(y)
      y[match(as.vector(x), y, 0L)]
    }
    common.element.in.all.sample.names = paste(Reduce(intersect2, strsplit(folder.info$samples, NULL)), collapse = '')

    ClosestMatch2 = function(string, stringVector){
      stringVector[stringdist::amatch(string, stringVector, maxDist=Inf)]
    }
    all.RNAseq.folders = list.files(paste0(folder.info$path, "/RNAseq/"))
    our.RNAseq.folder = ClosestMatch2(common.element.in.all.sample.names, all.RNAseq.folders)
    
    folder.info$RNAseq = tibble(sample = folder.info$samples) %>%
      mutate(
        normal = file.exists(paste0(folder.info$path, "/RNAseq/", our.RNAseq.folder, "/samples/", sample, "N/lib_rna/", sample, "N-rna.dedup.bam")),
        tumor = file.exists(paste0(folder.info$path, "/RNAseq/", our.RNAseq.folder, "/samples/", sample, "T/lib_rna/", sample, "T-rna.dedup.bam")),
        both = normal & tumor
      )
    
    #list of HLA plots
    folder.info$HLA = tibble(sample = folder.info$samples) %>%
      mutate(
        normal_coverage = file.exists(paste0(folder.info$path, "/HLA/plots/", sample, "N_coverage_plot.pdf")),
        tumor_coverage = file.exists(paste0(folder.info$path, "/HLA/plots/", sample, "T_coverage_plot.pdf")),
        both_coverage = normal_coverage & tumor_coverage,
        
        normal_rna = file.exists(paste0(folder.info$path, "/HLA/plots/", sample, "N-rna_coverage_plot.pdf")),
        tumor_rna = file.exists(paste0(folder.info$path, "/HLA/plots/", sample, "T-rna_coverage_plot.pdf")),
        both_rna = normal_rna & tumor_rna
      )
    
    
    

    
    #folder.info$content = list.files(folder.info$path)
    #collect the various sample names from the x-week-y folder
    #week.folder = paste0(folder.info$path, "/", folder.info$content[grep("week", folder.info$content)])
    #folder.info$vcf = sub(".vcf", "",
    #  list.files(week.folder, pattern = ".vcf$")
    #)
    #folder.info$bam = sub(".bam", "",
    #  list.files(week.folder, pattern = ".bam$")
    #)
    ##aggregate all unique sample names in this folder and the originals subfolder, regardless of file type
    #folder.info$samples = unique(
    #  gsub("-unaligned", "", 
    #       gsub("\\..*","",
    #            list.files(c(week.folder, paste0(folder.info$path, "/originals")))
    #)))
    #collect the sample names from the /originals folder
    #folder.info$samples = sub("-unaligned.bam", "",
    #  list.files(paste0(folder.info$path, "/originals"), pattern = ".bam$")
    #)
    
    #check presence of standard files. Returns TRUE or FALSE
    #checkfile = function(path, file){
    #  file.exists(paste0(path, file))
    #}
    
    #folder.info$QCsample = checkfile("/QCsample.html")
    #folder.info$coverage = checkfile("/Coverage.txt")
    #folder.info$gender   = checkfile("/gender.txt")
    #folder.info$NGSE     = checkfile("/NGSE.html")
    #folder.info$igv      = checkfile("/igv_session.xml")
  })
  
  
  
  #########################################################
  # Sidebar
  #########################################################
  
  #states the currently selected folder
  output$selected.folder <- renderPrint(
    if(go())cat("Selected folder: ", folder.info$path) else cat("No folder selected")
  )
  
  #lists the names of the .bam files found in the x-week-y folder.
  output$choose_sample <- renderUI(if(go()){
    radioButtons(
      inputId = "sampleselect",   #This does not do anything (yet).
      label = "Samples",
      choices = folder.info$samples
    )
  })
  
  
  
  #########################################################
  # Mainpanel 
  #########################################################
  
  output$title <- renderPrint(
    if(go())gsub("..*/", "", folder.info$path)
  )
  
  #custom function to change logic output to a green tick or red cross
  checklist = function(x){
    ifelse(x, "<font color=green>&#10004;</font>", "<font color=red>&times;</font>")
  }
  
  output$table1 <- renderTable(if(go()){
    tibble(
      "QCsample" = checklist(folder.info$QCsample),
      "Coverage" = checklist(folder.info$coverage),
      "gender" = checklist(folder.info$gender),
      "NGSE" = checklist(folder.info$NGSE),
      "igv" = checklist(folder.info$igv)
      )
  }, sanitize.text.function = function(x) x) #this line makes HTML work in the table
  
  output$table2 <- renderTable(if(go()){
    tibble("sample" = folder.info$samples) %>%
      mutate(
        bam = checklist(sample %in% folder.info$bam),
        vcf = checklist(sample %in% folder.info$vcf)
      )
  }, sanitize.text.function = function(x) x)
  
  
  
  #########################################################
  # Other
  #########################################################
  
  #terminate this app when the browser window is closed.
  session$onSessionEnded(function() {
    stopApp()
  })
  
}

shinyApp(ui = ui, server = server)