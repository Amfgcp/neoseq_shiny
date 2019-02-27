library(dplyr)
library(ggplot2)
# Data handling
# Hardcoded pep results folder
pep.info <- list("path" = paste0('O:/ImmunoGenomics/ngsdata', '/results'))

# List all .pep.txt files
pep.info$samples = list.files(pep.info$path, pattern = "^NIC.*pep.txt$")
pep.info$n = length(pep.info$samples)

# all_varType_count <- tribble()
variantTYPES_gathered <- tibble()

for (i in 1:pep.info$n) {
  pep_path = paste0(pep.info$path, "/", pep.info$samples[i])
  
  # Read file and filter by selected & expressed variants
  pepr <- read.table(pep_path, header = TRUE, sep = "\t") %>%
          select(varID, selection, Expressed, variantTYPE) %>%
          filter(selection=='yes', Expressed=='yes') %>%
          select(varID, variantTYPE)

  # split variantType column values and put them into new rows by varID,
  # then remove duplicates rows (so we don't count more than once for 
  # the same variant a given varType)
  pepr_split <- pepr %>%
                mutate(variantTYPE = strsplit(as.character(variantTYPE), ",")) %>% 
                tidyr::unnest(variantTYPE) %>% 
                unique()
          
  # filter unwanted varType rows
  pepr_relevant <- pepr_split %>%
                   subset(variantTYPE != 'coding_sequence_variant' &
                          variantTYPE != 'NMD_transcript_variant'  &
                          variantTYPE != 'splice_region_variant'   &
                          variantTYPE != '5_prime_UTR_variant'     &
                          variantTYPE != '3_prime_UTR_variant'     &
                          variantTYPE != 'intron_variant'          &
                          variantTYPE != 'non_coding_transcript_variant' &
                          variantTYPE != 'non_coding_transcript_exon_variant')

  # count occurrences of each variant type
  pepr_count <- pepr_relevant %>%
                select(variantTYPE) %>%
                group_by(variantTYPE) %>%
                mutate(occurrences = n()) %>%
                unique()
                
  # create column with sample name
  pepr_count$sample <- rep(sub(".25Lpep.txt", "", pep.info$samples[i]), nrow(pepr_count))
  
  
  # accumulates number of variantType occurrences
  # all_varType_count <- data.table::rbindlist(list(all_varType_count, pepr_count))[, lapply(.SD, sum, na.rm = TRUE), by = variantTYPE]
  
  variantTYPES_gathered = rbind(variantTYPES_gathered, as_tibble(pepr_count))
}
    
bar_plot <- ggplot(data=variantTYPES_gathered, aes(x=variantTYPE, y=occurrences, fill=sample)) +
  geom_bar(stat="identity") +
  theme_minimal()  +
  ylab("Number of Occurrences") +
  xlab("Variant Type") # + coord_cartesian(ylim = c(0, 750))

bar_plot
