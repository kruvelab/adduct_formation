library(tidyverse)
library(rcdklibs)
library(rcdk)

suppressMessages(
  list <- read_delim("data/pubchem_fingerprints_names.txt", delim = ",") %>%
    mutate(row = as.numeric(rowPC)) %>%
    select(row, description_PC)
)

pubchem_fingerprint = function(smiles) {
  tryCatch(
    {
      mol2 <- parse.smiles(smiles)[[1]]
      pubchem_fingerprints <- get.fingerprint(mol2,
                                              type = "pubchem")
      table <- as.data.frame(pubchem_fingerprints@bits) %>%
        mutate(fp = 1) %>%
        mutate(fp = as.numeric(fp))
      table = table %>%
        rename("row" = `pubchem_fingerprints@bits`)
      
      suppressMessages(
        datarow <- list %>%
          left_join(table) %>%
          mutate(row = paste("V", row, sep = "")) %>%
          select(-description_PC) %>%
          mutate(fp = case_when(
            is.na(fp) ~ 0,
            TRUE ~ fp))
      )
      
      datarow = datarow %>%
        spread(key = row, value = fp) %>%
        mutate(SMILES = smiles)
      return(datarow)
    },
    error = function(e) {
      return(tibble())
    }
  )
}