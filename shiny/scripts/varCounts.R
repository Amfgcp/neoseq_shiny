rm(list = ls(all = TRUE))

expressedAndTypeCountUI <- function(id) {

  ns <- NS(id)

  tagList(
    
    plotOutput(ns("plot_results_expressed")),
    checkboxGroupInput(ns("checkbox_results"), label = "Plot samples",
                       choices = c("yes", "no"), inline = TRUE),
    actionButton(ns("checkbox_results_apply_button"), label = "Apply"),
    actionLink(ns("selectall"), label = "Select All"),
    plotOutput(ns("plot_results_varType"))
    
  )

}

expressedAndTypeCount <- function(input, output, session, input_dir, folder_path) {
  rv <- reactiveValues()
  go = reactive(!is.na(input_dir()[2])) #True when a folder has been selected
  
  # When the directory is chosen, count the expressed variants
  # and the expressed var type
  observeEvent(input_dir(), if(go()) {
    
    cat(file=stderr(), "#### Debug", "input_dir() is: ", "\n")
    print(input_dir())
    
    cat(file=stderr(), "#### Debug_X", "Starting expressed variants count code ", "ok", "\n")
    # get path to results folder
    pep.info <- NULL
    cat(file=stderr(), "#### Debug_X", "file_path() is: ", folder_path(), "\n")
    pep.info$path <- paste0(folder_path(), "/results")
    cat(file=stderr(), "#### Debug", "Results path is: ", pep.info$path, "\n")
    
    # list all .pep.txt files
    pep.info$samples_files = list.files(pep.info$path, pattern = "^NIC.*pep.txt$")
    rv$pep.info$samples <- sub(".25Lpep.txt", "" , pep.info$samples_files)
    
    # Sort files and trim extension
    temp <- as.numeric(gsub("NIC","", rv$pep.info$samples))
    rv$pep.info$samples_sorted <- unique(rv$pep.info$samples[order(temp)])
    
    pep.info$n = length(pep.info$samples_files)
    
    cat(file=stderr(), "#### Debug", "pep.info$samples_files is: ", pep.info$samples_files, "\n")
    cat(file=stderr(), "#### Debug", "rv$pep.info$samples is: ", rv$pep.info$samples, "\n")
    cat(file=stderr(), "#### Debug", "rv$pep.info$samples_sorted is: ", rv$pep.info$samples_sorted, "\n")
    cat(file=stderr(), "#### Debug", "pep.info$n is: ", pep.info$n, "\n")
    
    
    # Check for already generated files
    if (file.exists("data/rv$variant.gather") && file.exists("data/rv$variantTYPES_gathered")) {
      
      cat(file=stderr(), "#### Debug", "RDS files exist, skipping computations", "\n")

      rv$variant.gather <- readRDS("data/rv$variant.gather")
      rv$variant.gather_filt <- rv$variant.gather

      rv$variantTYPES_gathered <- readRDS(file = "data/rv$variantTYPES_gathered")
      rv$variantTYPES_gathered_filt <- rv$variantTYPES_gathered

    } else {
    
      cat(file=stderr(), "#### Debug", "RDS files do NOT exist, doing computations", "\n")
      variant.data <- tibble() # expressed count
      rv$variantTYPES_gathered <- tibble() # varType count
      
      for (i in 1:pep.info$n) {
        if(i %% 3 == 0){cat(file=stderr(), "#### Debug", "iter ", i, "of pep.info$n", "\n")}
        
        pep_path = paste0(pep.info$path, "/", pep.info$samples_files[i])
        
        ## Expressed count
        # Read file and filter by selected & expressed variants
        pepr <- read.table(pep_path, header = TRUE, sep = "\t") %>%
          select(varID, selection, Expressed) %>%
          unique() %>%
          filter(selection == 'yes')
        
        peptide_count_yes = as.numeric(nrow(pepr))
        peptide_count_expressed = as.numeric(nrow(filter(pepr, Expressed == 'yes')))
        
        temp = tibble(sample=rv$pep.info$samples[i], variants_selected = peptide_count_yes, variants_expressed = peptide_count_expressed)
        variant.data = rbind(variant.data, temp)
        
        ## Variant Type count
        expr_pepr <- read.table(pep_path, header = TRUE, sep = "\t") %>%
          select(varID, selection, Expressed, variantTYPE) %>%
          filter(selection == 'yes', Expressed == 'yes') %>%
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
          subset( variantTYPE != 'coding_sequence_variant' &
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
      # Save file
      # cat(file=stderr(), "#### Debug", "Saving file... ", "\n")
      # saveRDS(rv$variant.gather, file = "data/rv$variant.gather")
      # cat(file=stderr(), "#### Debug", "File saved. ", "\n")
      rv$variant.gather_filt <- rv$variant.gather
      
      ## varType count reactive object to be plotted
      # Save file
      # saveRDS(rv$variantTYPES_gathered, file = "data/rv$variantTYPES_gathered")
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
  
  ## OUTPUTS ##
  
  # PLOTS #
  
  ## Expressed count
  output$plot_results_expressed <- renderPlot({
    ggplot(data=rv$variant.gather_filt, aes(x=factor(sample, level=rv$sample_ordered_levels), y=value, fill=Expressed)) +
      geom_bar(stat="identity") +
      theme_minimal()  +
      ylab("Number of coding change variants") +
      xlab("Sample") +
      scale_fill_grey() +
      theme_bw() +
      theme(axis.text.x=element_text(angle = 90, vjust = 0.5, size = 12),
            # remove the vertical grid lines
            panel.grid.major.x = element_blank(),
            panel.spacing = unit(0.02, "lines"),
            panel.border = element_blank(),
            strip.text = element_text(size = 12))
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
      theme(axis.text.x=element_text(angle = 90, vjust = 0.5, size = 12),
            # remove the vertical grid lines
            panel.grid.major.x = element_blank(),
            panel.spacing = unit(0.02, "lines"),
            panel.border = element_blank(),
            strip.text = element_text(size = 12))
  })
  
  # CHECKBOXES #
  
  observeEvent(input$checkbox_results_apply_button, {
    rv$variant.gather_filt <- filter(rv$variant.gather, sample %in% input$checkbox_results)
  })
  
  observe({
    if (input$selectall == 0) { # FIX $$
      return (NULL)
    } else if (input$selectall%%2 == 0) {
      updateCheckboxGroupInput(session, "checkbox_results", choices = rv$pep.info$samples_sorted, inline = TRUE)
    } else {
      updateCheckboxGroupInput(session, "checkbox_results", choices = rv$pep.info$samples_sorted, selected = rv$pep.info$samples_sorted, inline = TRUE)
    }
  })
  
  
}