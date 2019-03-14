library(dplyr)
library(ggplot2)
## Data handling
# hardcoded pep results folder
#pep.info <- list("path" = paste0('O:/ImmunoGenomics/ngsdata', '/results'))
pep.info <- list("path" = paste0('/mnt/patharchief/ImmunoGenomics/ngsdata', '/results'))

# list all .pep.txt files
pep.info$samples = list.files(pep.info$path, pattern = "^NIC.*pep.txt$")
pep.info$n = length(pep.info$samples)

variant.data = tibble()

for (i in 1:pep.info$n) {
  # NOTE: perhaps change blank columns to NA
  # NOTE: many blank  columns on "Expressed" col
  pep_path = paste0(pep.info$path, "/", pep.info$samples[i])
  
  pepr <- read.table(pep_path, header = TRUE, sep = "\t") %>%
          dplyr::select(varID, selection, Expressed) %>%
          unique() %>%
          filter(selection=='yes')
  
  peptide_count_yes = as.numeric(nrow(pepr))
  peptide_count_expressed = as.numeric(nrow(filter(pepr, Expressed=='yes')))
  
  temp = tibble(sample=pep.info$samples[i], variants_selected = peptide_count_yes, variants_expressed = peptide_count_expressed)
  variant.data = rbind(variant.data, temp)

}

variant.gather = variant.data %>%
                 mutate(variants_selected=variants_selected - variants_expressed) %>%
                 tidyr::gather(key="Expressed", value="value", -sample) %>%
                 mutate(sample=sub(".25Lpep.txt", "", sample), Expressed = ifelse(Expressed=="variants_expressed", "Yes", "No"))

## Bar Plot

temp <- as.numeric(gsub("NIC","", variant.gather$sample))
sample_ordered_levels <- unique(variant.gather$sample[order(temp)])
variant.gather$sample <- factor(variant.gather$sample, levels=sample_ordered_levels)
variant.gather$Expressed <- factor(variant.gather$Expressed, levels=c("Yes","No"))

subset <- variant.gather[!(variant.gather$sample %in% c("NIC12", "NIC13")),]

bar_plot <- ggplot(data=subset, aes(x=sample, y=value, fill=Expressed)) +
  geom_bar(stat="identity") +
  ylab("Number of coding change variant") +
  xlab("") +
  scale_fill_grey() +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust = 1, size=6),  
        # remove the vertical grid lines
        panel.grid.major.x = element_blank(), 
        panel.spacing = unit(0.02, "lines"),
        panel.border = element_blank(),
        strip.text = element_text(size=12))

bar_plot
ggsave("~/Dropbox/Presentations/neoseq/expressedVar_Nics_MMRp.jpeg", width = 6, height = 4)


