source("predict_adduct_formation.R")

SMILES = tibble("SMILES" = c("CCOc1cc(ccc1C(=O)OC)NC(=O)C", "c1cc(ccc1N)S(=O)(=O)Nc1cncc(n1)Cl"))

SMILES %>%
  group_by(SMILES) %>%
  mutate(test = predict_adduct_formation(SMILES)) %>%
  ungroup()
