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
          # Expressed count
          plotOutput("plot_results_expressed"),
          checkboxGroupInput(inputId = "checkbox_results", label = "Plot samples",
                             choices = , inline = TRUE),
          actionButton(inputId = "checkbox_results_apply_button", label = "Apply"),
          actionLink(inputId = "selectall", label = "Select All"),
          # varType count
          plotOutput("plot_results_varType")
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
  # roots = c(unique(substr(list.dirs(paste0(LETTERS, ":/"), recursive = F), 1, 3)), getwd(), "C:/Cloud/OneDrive - Universidade de Lisboa/Projects/neoseq_shiny/example-samples")
  roots = c(unique(substr(list.dirs(paste0(LETTERS, ":/"), recursive = F), 1, 3)), getwd(), "O:/ImmunoGenomics/ngsdata")
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
  # When the directory is chosen, count the expressed variants
  # and the expressed var type
  observeEvent(input$dir, if(go()){
    cat(file=stderr(), "#### Debug", "Starting expressed variants count code ", "ok", "\n")
    # get path to results folder
    pep.info <- NULL
    pep.info$path <- paste0(folder$path, "/results")
    cat(file=stderr(), "#### Debug", "Results path is: ", pep.info$path, "\n")
    
    # list all .pep.txt files
    pep.info$samples_files = list.files(pep.info$path, pattern = "^NIC.*pep.txt$")
    rv$pep.info$samples <- sub(".25Lpep.txt", "" , pep.info$samples_files)
    
    # Sort files and trim extension
    temp <- as.numeric(gsub("NIC","", rv$pep.info$samples))
    rv$pep.info$samples_sorted <- unique(rv$pep.info$samples[order(temp)])
    
    pep.info$n = length(pep.info$samples_files)
    # cat(file=stderr(), "#### Debug", "pep.info$samples_files is: ", pep.info$samples_files, "\n")
    # cat(file=stderr(), "#### Debug", "rv$pep.info$samples is: ", rv$pep.info$samples, "\n")
    cat(file=stderr(), "#### Debug", "rv$pep.info$samples_sorted is: ", rv$pep.info$samples_sorted, "\n")
    # cat(file=stderr(), "#### Debug", "pep.info$n is: ", pep.info$n, "\n")
    
    if (file.exists("data/rv$variant.gather") && file.exists("data/rv$variantTYPES_gathered")) {
      cat(file=stderr(), "#### Debug", "RDS files exist, skipping computations", "\n")
      
      rv$variant.gather <- readRDS("data/rv$variant.gather")
      rv$variant.gather_filt <- rv$variant.gather
      
      rv$variantTYPES_gathered <- readRDS(file = "data/rv$variantTYPES_gathered")
      rv$variantTYPES_gathered_filt <- rv$variantTYPES_gathered
      
    } else {
      cat(file=stderr(), "#### Debug", "RDS files do NOT exist, doing computations", "\n")
      variant.data = tibble() # expressed count
      rv$variantTYPES_gathered <- tibble() # varType count
      
      for (i in 1:pep.info$n) {
        if(i %% 3 == 0){cat(file=stderr(), "#### Debug", "iter ", i, "of pep.info$n", "\n")}
        
        # NOTE: perhaps change blank columns to NA
        # NOTE: many blank  columns on "Expressed" col
        pep_path = paste0(pep.info$path, "/", pep.info$samples_files[i])
        
        ## Expressed count
        # Read file and filter by selected & expressed variants
        pepr <- read.table(pep_path, header = TRUE, sep = "\t") %>%
          select(varID, selection, Expressed) %>%
          unique() %>%
          filter(selection=='yes')
        
        peptide_count_yes = as.numeric(nrow(pepr))
        peptide_count_expressed = as.numeric(nrow(filter(pepr, Expressed=='yes')))
        
        temp = tibble(sample=rv$pep.info$samples[i], variants_selected = peptide_count_yes, variants_expressed = peptide_count_expressed)
        # temp = tibble(sample=rv$pep.info$samples_sorted[i], variants_selected = peptide_count_yes, variants_expressed = peptide_count_expressed)
        variant.data = rbind(variant.data, temp)
        
        ## Variant Type count
        expr_pepr <- read.table(pep_path, header = TRUE, sep = "\t") %>%
          select(varID, selection, Expressed, variantTYPE) %>%
          filter(selection=='yes', Expressed=='yes') %>%
          select(varID, variantTYPE)
        
        # split variantType column values and put them into new rows by varID,
        # then remove duplicates rows (so we don't count more than once for
        # the same variant a given varType)
        expr_pepr_split <- expr_pepr %>%
          mutate(variantTYPE = strsplit(as.character(variantTYPE), ",")) %>%
          tidyr::unnest(variantTYPE) %>%
          unique()
        
        # filter unwanted varType rows
        expr_pepr_relevant <- expr_pepr_split %>%
          subset(variantTYPE != 'coding_sequence_variant' &
                   variantTYPE != 'NMD_transcript_variant'  &
                   variantTYPE != 'splice_region_variant'   &
                   variantTYPE != '5_prime_UTR_variant'     &
                   variantTYPE != '3_prime_UTR_variant'     &
                   variantTYPE != 'intron_variant'          &
                   variantTYPE != 'non_coding_transcript_variant' &
                   variantTYPE != 'non_coding_transcript_exon_variant')
        
        # count occurrences of each variant type
        expr_pepr_count <- expr_pepr_relevant %>%
          select(variantTYPE) %>%
          group_by(variantTYPE) %>%
          mutate(occurrences = n()) %>%
          unique()
        
        # count occurrences of each variant type
        expr_pepr_count <- expr_pepr_relevant %>%
          select(variantTYPE) %>%
          group_by(variantTYPE) %>%
          mutate(occurrences = n()) %>%
          unique()
        
        # create column with sample name
        expr_pepr_count$sample <- rep(sub(".25Lpep.txt", "", rv$pep.info$samples[i]), nrow(expr_pepr_count))
        
        rv$variantTYPES_gathered = rbind(rv$variantTYPES_gathered, as_tibble(expr_pepr_count))
      }
      cat(file=stderr(), "#### Debug", "for loop ended", "ok", "\n")
 
      ## Expressed count
      rv$variant.gather = variant.data %>%
        mutate(variants_selected=variants_selected - variants_expressed) %>%
        tidyr::gather(key="Expressed", value="value", -sample) %>%
        mutate(Expressed = ifelse(Expressed=="variants_expressed", "Yes", "No"))
      cat(file=stderr(), "#### Debug", "rv$variant.gather created", "ok", "\n")
      
      # switch levels order
      rv$variant.gather$Expressed <- factor(rv$variant.gather$Expressed, levels=c("Yes","No"))
      
      # expressed count reactive object to be plotted
      cat(file=stderr(), "#### Debug", "Saving file... ", "\n")
      saveRDS(rv$variant.gather, file = "data/rv$variant.gather")
      cat(file=stderr(), "#### Debug", "File saved. ", "\n")
      rv$variant.gather_filt <- rv$variant.gather
      
      ## varType count reactive object to be plotted
      saveRDS(rv$variantTYPES_gathered, file = "data/rv$variantTYPES_gathered")
      rv$variantTYPES_gathered_filt <- rv$variantTYPES_gathered
    }
    
    ## Common
    # Sort sample names to achieve a sorted xx axis
    # cat(file=stderr(), "#### Debug", "rv$variant.gather$sample is: ", rv$variant.gather$sample, "\n")
    temp <- as.numeric(gsub("NIC","", rv$variant.gather$sample))
    # cat(file=stderr(), "#### Debug", "temp is: ", temp, "\n")
    rv$sample_ordered_levels <- unique(rv$variant.gather$sample[order(temp)])
    # cat(file=stderr(), "#### Debug", "ordered levels are: ", rv$sample_ordered_levels,"\n")
    cat(file=stderr(), "#### Debug", "levels for plot created", "ok", "\n")
    
    updateCheckboxGroupInput(session, "checkbox_results", choices = rv$pep.info$samples_sorted, inline = TRUE, selected = rv$pep.info$samples_sorted )
    cat(file=stderr(), "#### Debug", "1st update to checkbox choices", "ok", "\n")

  })
  
  ## Expressed count
  output$plot_results_expressed <- renderPlot({
    ggplot(data=rv$variant.gather_filt, aes(x=factor(sample, level=rv$sample_ordered_levels), y=value, fill=Expressed)) +
      geom_bar(stat="identity") +
      theme_minimal()  +
      ylab("Number of coding change variants") +
      xlab("Sample") +
      scale_fill_grey() +
      theme_bw() +
      theme(axis.text.x=element_text(angle=90, vjust = 0.5, size=12),
            # remove the vertical grid lines
            panel.grid.major.x = element_blank(),
            panel.spacing = unit(0.02, "lines"),
            panel.border = element_blank(),
            strip.text = element_text(size=12))
  })
  
  observeEvent(input$checkbox_results_apply_button, {
    
    rv$variant.gather_filt <- filter(rv$variant.gather, sample %in% input$checkbox_results)
    # Sort sample names for sorted plot xx axis
    # temp <- as.numeric(gsub("NIC","", rv$variant.gather_filt$sample))
    # rv$sample_ordered_levels <- unique(rv$variant.gather_filt$sample[order(temp)])
  })
  
  ## varType count
  output$plot_results_varType <- renderPlot({
    ggplot(data=rv$variantTYPES_gathered_filt, aes(x=variantTYPE, y=occurrences, fill=sample)) +
      geom_bar(stat="identity") +
      theme_minimal()  +
      ylab("Number of Occurrences") +
      xlab("Variant Type") +
      scale_fill_grey() +
      theme_bw() +
      theme(axis.text.x=element_text(angle=90, vjust = 0.5, size=12),  
            # remove the vertical grid lines
            panel.grid.major.x = element_blank(), 
            panel.spacing = unit(0.02, "lines"),
            panel.border = element_blank(),
            strip.text = element_text(size=12))
  })

  # observeEvent(input$checkbox_results_apply_button,
  #              rv$variantTYPES_gathered_filt <- filter(rv$variantTYPES_gathered, sample %in% input$checkbox_results)
  # )
  
  observe({
    if (input$selectall == 0) return(NULL) 
    else if (input$selectall%%2 == 0) {
      updateCheckboxGroupInput(session, "checkbox_results", choices = rv$pep.info$samples_sorted, inline = TRUE)
    }
    else {
      updateCheckboxGroupInput(session, "checkbox_results", choices = rv$pep.info$samples_sorted, selected = rv$pep.info$samples_sorted, inline = TRUE)
    }
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