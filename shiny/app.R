#install.packages(c("shiny", "dplyr", "xtable", "shinyFiles", "stringdist"))
library(shiny)
library(dplyr)
library(xtable) #to make html work in the tables. Doesn't work with xtable::
#library(shinyFiles) #this library is used, but does not need to be loaded since we only call it via shinyFiles::[function]
#library(stringdist) #for ClosestMatch2's amatch
library(ggplot2)

#open a browser for the app to be displayed in. URL must be identical to the URL in run.vbs
browserpath = "C:/Program Files/Internet Explorer/iexplore.exe" #paste0(sub("/Shiny", "", getwd()), "asd")
# NOTE: disabled browser opening while developing:
# browseURL("http://127.0.0.1:7777", browserpath)

source("scripts/varCounts.R")

                    ##################
                    # User Interface #
                    ##################
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
      uiOutput("title", style = "margin-top: 15px; margin-bottom: 20px;"),
      
      tabsetPanel(id = "top_tabs",
        tabPanel(
          "Overview",
          tableOutput("overview")
        ),
        tabPanel(
          ".bam files",
          tableOutput("bams")
        ),
        tabPanel(
          ".vcf files",
          tableOutput("vcfs")
        ),
        tabPanel(
          "RNAseq",
          tableOutput("rnas")
        ),
        tabPanel(
          "HLA plots",
          tableOutput("hlas"),
          uiOutput("sel_hla"),
          uiOutput("hla_pdf")
          #tags$iframe(style="height:400px; width:100%; scrolling=yes", src="O:/Transfer/Siebren/neoSeq/HLA/plots/NIC3N-rna_coverage_plot.pdf")
        ),
        tabPanel(
          "Results",
          expressedAndTypeCountUI("var_counts")
        )
      )
    )
  )
)

                    ################
                    # Server Logic #
                    ################

server = function(input, output, session) {
  rv <- reactiveValues()
  
  #########################################################
  # folder input
  #########################################################
  
  # 1) create a named list of places where users can start searching
  # roots = c(unique(substr(list.dirs(paste0(LETTERS, ":/"), recursive = F), 1, 3)), getwd(), "O:/Transfer/Siebren/neoSeq")
  # roots = c(unique(substr(list.dirs(paste0(LETTERS, ":/"), recursive = F), 1, 3)), getwd(), "C:/Users/amfgcpaulo/OneDrive - Universidade de Lisboa/Projects/neoseq_shiny/example-samples")
  roots = c(unique(substr(list.dirs(paste0(LETTERS, ":/"), recursive = F), 1, 3)), getwd(), "C:/Cloud/OneDrive - Universidade de Lisboa/Projects/neoseq_shiny/example-samples")
  # roots = c(unique(substr(list.dirs(paste0(LETTERS, ":/"), recursive = F), 1, 3)), getwd(), "O:/ImmunoGenomics/ngsdata")
  names(roots) <- c(roots[1:(length(roots)-2)], getwd(), "demo_location")
  
  #2) use the named list to allow users to select a folder. This function is coupled to shinyDirButton in the ui section.
  shinyFiles::shinyDirChoose(input, 'dir',  updateFreq = 0, roots = roots, defaultRoot = "demo_location")
  go = reactive(!is.na(input$dir[2])) #True when a folder has been selected
  
  
  #3) harvest data from the user's selected folder
  folder <- reactiveValues()
  observeEvent(input$dir, if(go()){
    folder$path <- shinyFiles::parseDirPath(roots, input$dir)
    
    cat(file = stderr(), "#### Debug. The folder$path is: ", folder$path, "\n")
    #list of samples
    TNsamples = sub("_L..*", "",
      list.files(paste0(folder$path, "/fastq"), pattern = ".fastq.gz$")
    )
    folder$samples = unique(substr(TNsamples,1,nchar(TNsamples)-1))
    
    #table of bam samples
    folder$bam = tibble(sample = folder$samples) %>%
      mutate(
       normal = file.exists(paste0(folder$path, "/bam/", sample, "N.dedup.recal.bam")),
       tumor = file.exists(paste0(folder$path, "/bam/", sample, "T.dedup.recal.bam")),
       both = normal & tumor
      )
    
    #list of vcf samples
    folder$vcf = sub("_CombineVariants.vcf", "",
      list.files(paste0(folder$path, "/vcf"), pattern = ".vcf$")
    )
    
    #table of RNAseq samples
    intersect2 <- function (x, y){
      y <- as.vector(y)
      y[match(as.vector(x), y, 0L)]
    }
    common.element.in.all.sample.names = paste(Reduce(intersect2, strsplit(folder$samples, NULL)), collapse = '')

    ClosestMatch2 = function(string, stringVector){
      stringVector[stringdist::amatch(string, stringVector, maxDist=Inf)]
    }
    all.RNAseq.folders = list.files(paste0(folder$path, "/RNAseq/"))
    our.RNAseq.folder = ClosestMatch2(common.element.in.all.sample.names, all.RNAseq.folders)
    
    folder$RNAseq = tibble(sample = folder$samples) %>%
      mutate(
        normal = file.exists(paste0(folder$path, "/RNAseq/", our.RNAseq.folder, "/samples/", sample, "N/lib_rna/", sample, "N-rna.dedup.bam")),
        tumor = file.exists(paste0(folder$path, "/RNAseq/", our.RNAseq.folder, "/samples/", sample, "T/lib_rna/", sample, "T-rna.dedup.bam")),
        both = normal & tumor
      )
    
    #table of HLA plots
    folder$HLA = tibble(sample = folder$samples) %>%
      mutate(
        normal_coverage = file.exists(paste0(folder$path, "/HLA/plots/", sample, "N_coverage_plot.pdf")),
        tumor_coverage = file.exists(paste0(folder$path, "/HLA/plots/", sample, "T_coverage_plot.pdf")),
        both_coverage = normal_coverage & tumor_coverage,
        
        normal_rna = file.exists(paste0(folder$path, "/HLA/plots/", sample, "N-rna_coverage_plot.pdf")),
        tumor_rna = file.exists(paste0(folder$path, "/HLA/plots/", sample, "T-rna_coverage_plot.pdf")),
        both_rna = normal_rna & tumor_rna
      )
    
    #table of Results
    folder$results = tibble(sample = folder$samples) %>%
      mutate(
        present = file.exists(paste0(folder$path, "/results/", sample, ".25Lpep.txt"))
      )
      
  })
  
  
  #########################################################
  # Sidebar
  #########################################################
  
  #states the currently selected folder
  output$selected.folder <- renderPrint(
    if(go())cat("Selected folder: ", folder$path) else cat("No folder selected")
  )
  
  #lists the names of the .bam files found in the x-week-y folder.
  output$choose_sample <- renderUI(if(go()){
    radioButtons(
      inputId = "sampleselect",   #This does not do anything (yet).
      label = "Samples",
      choices = folder$samples
    )
  })
  
  
  #########################################################
  # Mainpanel 
  #########################################################
  
  output$title <- renderText(
    if(go()) paste0("<b>", gsub("..*/", "", folder$path), "</b>")
    else "<br>"
  )
  
  #change logic output to a green tick or red cross, and any other output in orange
  checklist = function(x){
    ifelse(x == T, "<font color=green>&#10004;</font>", 
           ifelse(x == F, "<font color=red>&times;</font>", 
                  paste0("<font color=orange>", x, "</font>")))
  }
  
  ####### overview tab #######
  output$overview <- renderTable(if(go()){
    tibble("sample" = folder$samples) %>%
      mutate(
        ".bam" = checklist(ifelse(folder$bam$both, T, 
                                  ifelse(folder$bam$tumor, "Tumor only",
                                         ifelse(folder$bam$normal, "Normal only",
                                                F)))),
        ".vcf" = checklist(sample %in% folder$vcf),
        "RNAseq" = checklist(ifelse(folder$RNAseq$both, T, 
                                    ifelse(folder$RNAseq$tumor, "Tumor only",
                                           ifelse(folder$RNAseq$normal, "Normal only",
                                                  F)))),
        "HLA coverage plots" = checklist(ifelse(folder$HLA$both_coverage, T, 
                                                ifelse(folder$HLA$tumor_coverage, "Tumor only", 
                                                       ifelse(folder$HLA$normal_coverage, "Normal only", 
                                                              F)))),
        "HLA rna plots" = checklist(ifelse(folder$HLA$both_rna, T, 
                                           ifelse(folder$HLA$tumor_rna, "Tumor only", 
                                                  ifelse(folder$HLA$normal_rna, "Normal only",
                                                         F))))
      )
  }, sanitize.text.function = function(x) x) #this line makes HTML work in the table
  
  ####### bam tab #######
  output$bams <- renderTable(if(go()){
    folder$bam %>%
      select(-both) %>%
      mutate(
        normal = checklist(normal),
        tumor = checklist(tumor)
      )
  }, sanitize.text.function = function(x) x)
  
  ####### vcf tab #######
  output$vcfs <- renderTable(if(go()){
    tibble("sample" = folder$samples) %>%
      mutate(
        vcf = checklist(sample %in% folder$vcf)
      )
  }, sanitize.text.function = function(x) x)
  
  ####### RNA tab #######
  output$rnas <- renderTable(if(go()){
    folder$RNAseq %>%
      select(-both) %>%
      mutate(
        normal = checklist(normal),
        tumor = checklist(tumor)
      )
  }, sanitize.text.function = function(x) x)
  
  ####### HLA tab #######
  #1) overview table
  output$hlas <- renderTable(if(go()){
    folder$HLA %>%
      select(-both_coverage, -both_rna) %>%
      mutate(
        normal_coverage = checklist(normal_coverage),
        tumor_coverage = checklist(tumor_coverage),
        normal_rna = checklist(normal_rna),
        tumor_rna = checklist(tumor_rna)
      )
  }, sanitize.text.function = function(x) x)
  
  #2) select specific pdf to view for the selected sample
  output$sel_hla <- renderUI(if(go()){
    radioButtons(
      inputId = "sel_hla2",
      label = "",
      choices = c("normal coverage", "normal rna", "tumor coverage", "tumor rna"),
      inline = T
    )
  })
  
  #3) pdf output
  output$hla_pdf <- renderUI(if(go()){
    selected_file = paste0(folder$path, "/HLA/plots/", input$sampleselect, sub("normal ", "N", sub("tumor ", "T", sub("rna", "-rna_coverage", sub("coverage", "_coverage", input$sel_hla2)))), "_plot.pdf")
    if(file.exists(selected_file)){
      tags$iframe(style="height:600px; width:100%", src=selected_file)
    } else {
      "Selected file does not exist."
    }
  })
  
  ####### Results tab #######
  var_counts <- callModule(expressedAndTypeCount, "var_counts", reactive(input$dir), reactive(folder$path))

  
  #########################################################
  # Other
  #########################################################
  
  #terminate this app when the browser window is closed.
  session$onSessionEnded(function() {
    stopApp()
  })
  
}

shinyApp(ui = ui, server = server)