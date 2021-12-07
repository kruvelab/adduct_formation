library(tidyverse)
library(caret)
library(caTools)
library(gmodels)
library(mccr)
library(yardstick)
source('code/my_theme.R')

classifier_svmPoly = readRDS("models/classifier_svmPoly.rds")

#Validation------

validation = read_delim("data/validation/ac8b04567_si_002_CID.csv",
                        delim = ",",
                        col_names = TRUE)

validation = validation %>%
  spread(key = Adduct, value = CCS)

validation = validation %>%
  mutate(M_Na = case_when(
    is.na(`[M+Na]`) ~ 0,
    TRUE ~ 1
  )) %>%
  select(Name, M_Na, PubChemCID, ExactMass)

validation_fingerprints = read_delim("data/validation/PubChem_fingerprints_bond_properties_validation.csv",
                                     delim = ",",
                                     col_names = TRUE) 

validation = validation %>%
  left_join(validation_fingerprints) %>%
  na.omit()

validation = validation %>%
  mutate(M_Na_pred_prob = predict(classifier_svmPoly, newdata = validation, type = "prob")[,2]) %>%
  mutate(M_Na_pred = case_when(
    M_Na_pred_prob > 0.5 ~ 1,
    TRUE ~ 0
  ))

write_delim(validation,
            "data/validation/Predicted_adduct_formation_validation_set.csv",
            delim = ",")

val_plot =
ggplot(data = validation %>%
         group_by(ExactMass) %>%
         mutate(count = n(),
                prob_dif = max(M_Na_pred_prob) -min(M_Na_pred_prob)) %>%
         ungroup() %>%
         filter(count > 1 & prob_dif > 0.1)) +
  geom_point(mapping = aes(x = ExactMass,
                           y = M_Na_pred_prob,
                           text = Name,
                           color = M_Na),
             size = 2) +
  geom_abline(slope = 0, intercept = 0.500) +
  my_theme

plotly::ggplotly(val_plot)

mean(validation$M_Na == validation$M_Na_pred)
#74.7%
CrossTable(validation$M_Na, validation$M_Na_pred)
#false positive 0.454
79/174
#false neg 0.166
66/398

validation2 = validation %>%
  mutate(M_Na = case_when(
    M_Na == 1 ~ 0,
    TRUE ~ 1),
    M_Na_pred = case_when(
      M_Na_pred == 1 ~ 0,
      TRUE ~ 1))


caret::precision(factor(validation2$M_Na_pred), factor(validation2$M_Na))
#0.808
caret::recall(factor(validation2$M_Na_pred), factor(validation2$M_Na))
#0.834
MLmetrics::F1_Score(y_true = validation2$M_Na, y_pred = validation2$M_Na_pred)
#82.1%
mccr::mccr(validation2$M_Na_pred, validation2$M_Na)
#0.389
bal_accuracy_vec(factor(validation2$M_Na_pred), factor(validation2$M_Na))
#0.700