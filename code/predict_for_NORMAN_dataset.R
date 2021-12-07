library(plotly)
library(cowplot)
source("code/PubChem_fingerprint_calculation.R")
source("code/calculating_bond_properties_from_graph.R")
source("code/my_theme.R")

NORMAN = read_delim("data/NORMAN/susdat_2019-12-06-112040.csv",
                    delim = ",",
                    col_names = TRUE) %>%
  rename(PubChemCID = PubChem_CID,
         CAS = `CAS_RN Dashboard`)

NORMAN = NORMAN %>%
  group_by(Monoiso_Mass) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  filter(count > 1 & Monoiso_Mass > 100 & Monoiso_Mass < 1000) %>%
  select(Monoiso_Mass, MS_Ready_SMILES) %>%
  na.omit()

Bond_descriptors_summary <- tibble()
for (smiles in NORMAN$MS_Ready_SMILES) {
  print(smiles)
  bond_descriptors_this_SMILES = bond_descriptors(smiles)
  pubchem_descriptors_this_SMILES = pubchem_fingerprint(smiles)
  if (length(bond_descriptors_this_SMILES) != 0 & length(pubchem_descriptors_this_SMILES) != 0) {
    Bond_descriptors_summary <- Bond_descriptors_summary %>%
      bind_rows(bond_descriptors_this_SMILES %>%
                  left_join(pubchem_descriptors_this_SMILES))
  }
}

write_delim(Bond_descriptors_summary,
            "data/NORMAN/PubChem_fingerprints_bond_properties_NORMAN.csv",
            delim = ",")

Bond_descriptors_summary = read_delim("data/NORMAN/PubChem_fingerprints_bond_properties_NORMAN.csv",
                                      delim = ",")

NORMAN = NORMAN  %>%
  rename(SMILES = MS_Ready_SMILES) %>%
  left_join(Bond_descriptors_summary) %>%
  na.omit()

#making predictions----
classifier_svmPoly = readRDS("models/classifier_svmPoly.rds")

NORMAN_pred = NORMAN %>%
  mutate(M_Na_pred = predict(classifier_svmPoly, newdata = NORMAN, type = "prob")[,2])

write_delim(NORMAN_pred,
            "data/NORMAN/adduct_formation_NORMAN.csv",
            delim = ",")

