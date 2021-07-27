library(rcdk)
library(tidyverse)
library(ape)
source("my_theme.R")
setwd("C:/Users/annel/OneDrive - Kruvelab/Kruvelab/computational/IE mudeli script ja failid/adduct_formation/data/training")

dataset = read_delim("balanced_adduct_data_training.csv",
                      delim = ",",
                      col_names = TRUE)

#clustering of the compounds based on the fingerprints
d = dist(dataset %>% select(-PubChemCID, -M_Na, -Lab), method = 'binary')
hc = hclust(d, method = 'average')
plot(hc)
y_hc = cutree(hc, 3) 
dataset_with_class <- tibble(dataset, y_hc)

dend = as.dendrogram(hc)
dend_data = ggdendro::dendro_data(dend, type = "rectangle")
names(dend_data)
dend_data$labels
head(dend_data$segments)
segments_data = dend_data$segments

dataset_with_class = dataset_with_class %>%
    mutate(label = row_number(),
           Lab = case_when(
               Lab == "Corey" ~ "Broeckling",
               Lab == "SU" ~ "Costalunga",
               Lab == "UT" ~ "Liigand",
               TRUE ~ Lab
           ))

labels_data = dend_data$labels %>%
    mutate(label = as.numeric(label)) %>%
    left_join(dataset_with_class)

p <- ggplot(data = segments_data ) + 
    geom_segment(mapping = aes(x = x, 
                               y = y, 
                               xend = xend, 
                               yend = yend),
                 color = basecolor,
                 size = 0.2) +
    geom_point(data = labels_data, 
              mapping = aes(x = x, 
                            y = y,
                            color = Lab),
              shape = 15,
              size = 5) +
    scale_color_manual(values=c("#d00000", "#ffba08", "#3f88c5", "#032b43", "#136f63")) +
    labs(x = "", y = "distance") +
    my_theme +
    theme(aspect.ratio = 1/3,
          legend.position = "bottom",
          axis.line.x = element_blank(),
          axis.text.x = element_blank())
    

print(p)

ggsave("dendrogram_balanced_data_from_fingerprints_ggplot.svg",
       width = 16,
       height = 16/773*328,
       units = "cm")

