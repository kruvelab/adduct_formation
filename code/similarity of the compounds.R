library(rcdk)
library(tidyverse)
library(ape)

dataset <- read_delim("balanced_adduct_data_training.csv",
                      delim = ",",
                      col_names = TRUE)

#clustering of the compounds based on the fingerprints
d = dist(dataset %>% select(-PubChemCID, -M_Na, -Lab), method = 'binary')
hc = hclust(d, method = 'average')
plot(hc)
y_hc = cutree(hc, 3) 
dataset_with_class <- tibble(dataset, y_hc)

colors = c("red", "blue", "green")
plot(as.phylo(hc), 
     type = "fan", 
     tip.color = colors[as.factor(dataset_with_class$Lab)],
     cex = 0.6,
     #label.offset = 1,
     #no.margin = TRUE
     )
svg("dendrogram_balanced_data_from_fingerprints.svg",
    height = 16,
    width = 16)

colors = c("red", "blue")
plot(as.phylo(hc), 
     type = "fan", 
     tip.color = colors[as.factor(dataset_with_class$M_Na)],
     cex = 0.6,
     #label.offset = 1,
     #no.margin = TRUE
)
svg("dendrogram_balanced_data_from_fingerprints_Adducts_1_0.svg",
    height = 16,
    width = 16)

library(proxy)
s = proxy::as.simil(d)
s = as_tibble(s)
s = as.matrix(s)
write_delim(s,
            "similarity.csv",
            delim = ",")
small_s = s[1:50, 1:50]
heatmap(s)

similarity = read_delim("similarity.csv",
                        delim = ",",
                        col_names = TRUE)

svg("heatmap_balanced_data_from_fingerprints.svg",
    height = 16,
    width = 16)
