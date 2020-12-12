library(tidyverse)
library(caret)
library(gmodels)
library(caTools)
source('my_theme.R')

setwd("C:/Users/annel/OneDrive - Kruvelab/Kruvelab/computational/PubChem")


#Adduct data ----
#UT data-----
UT_data <- read_delim("UT_adducts_CID.csv",
                      delim = ",",
                      col_names = TRUE) %>%
  na.omit()

UT_data <- UT_data %>%
  mutate(M_H = case_when(
    Class == 0 ~ 1,
    Class == 2 ~ 1,
    TRUE ~ 0
  ),
  M_Na = case_when(
    Class == 1 ~ 1,
    Class == 2 ~ 1,
    TRUE ~ 0
  ),
  Lab = "UT")

UT_data  %>% 
  group_by(M_Na) %>%
  summarise(n())

#there are about 2x so (224 vs 127) many compounds not forming an adduct as there are ones forming an adduct
UT_data_adduct <- sample_n(UT_data %>%
                            filter(M_Na == 1),
                          size = 100,
                          replace = FALSE)

UT_data_No_adduct <- sample_n(UT_data %>%
                             filter(M_Na == 0),
                           size = 100,
                           replace = FALSE)
UT_data_balanced <- UT_data_adduct %>%
  bind_rows(UT_data_No_adduct)

#Corey data -----

Corey_data <- read_delim("Corey_adducts_CID.csv",
                      delim = ",",
                      col_names = TRUE)

Corey_data <- Corey_data %>%
  mutate(M_Na = case_when(
                    Class == 1 ~ 1,
                    Class == 2 ~ 1,
                    TRUE ~ 0),
         Lab = "Corey")

Corey_data  %>% 
  group_by(M_Na) %>%
  summarise(n())

#there are about 6x so (72 vs 525) many compounds forming an adduct as there are ones forming an adduct
Corey_data_No_adduct <- sample_n(Corey_data %>%
                             filter(M_Na == 0),
                           size = 100,
                           replace = TRUE)

Corey_data_adduct <- sample_n(Corey_data %>%
                                   filter(M_Na == 1),
                                 size = 100,
                                 replace = FALSE)

Corey_data_balanced <- Corey_data_adduct %>%
  bind_rows(Corey_data_No_adduct)

#SU data -----

SU_data <- read_delim("SU_adducts_CID.csv",
                      delim = ",",
                      col_names = TRUE)

SU_data <- SU_data %>%
  mutate(M_Na = case_when(
                  SlopeNa == 0 ~ 0,
                  TRUE ~ 1),
        M_H = case_when(
                  SlopeH == 0 ~ 0,
                  TRUE ~ 1),
        Lab = "SU")

SU_data  %>% 
  group_by(M_Na) %>%
  summarise(n())

#this is an almost even dataset (53 vs 41)

SU_data_adduct <- sample_n(SU_data %>%
                                filter(M_Na == 1),
                              size = 100,
                              replace = TRUE)

SU_data_No_adduct <- sample_n(SU_data %>%
                                filter(M_Na == 0),
                              size = 100,
                              replace = TRUE)
SU_data_balanced <- SU_data_adduct %>%
  bind_rows(SU_data_No_adduct)

#Celma data----

Celma_data = read_delim("TableS1_CCS_RT_mz_DB_IUPA_charge.csv",
                              delim = ",",
                              col_names = TRUE)

Celma_data = Celma_data %>%
  filter(AdductSpecies == "[M+Na]" | AdductSpecies == "[M+H]") %>%
  select(PubChemCID, AdductSpecies, CCS) %>%
  group_by(PubChemCID, AdductSpecies) %>%
  summarise(CCS = mean(CCS)) %>%
  ungroup()

Celma_data = Celma_data %>%
  spread(key = AdductSpecies, value = CCS)

Celma_data = Celma_data %>%
  mutate(M_Na = case_when(
    is.na(`[M+Na]`) ~ 0,
    TRUE ~ 1
  ))

Celma_data = Celma_data %>%
  select(PubChemCID, M_Na) %>%
  mutate(Lab = "Celma")


Celma_data  %>% 
  group_by(M_Na) %>%
  summarise(n())
#this is an almost even dataset (274 vs 247)

Celma_data_adduct <- sample_n(Celma_data %>%
                             filter(M_Na == 1),
                           size = 200,
                           replace = FALSE)

Celma_data_No_adduct <- sample_n(Celma_data %>%
                                filter(M_Na == 0),
                              size = 200,
                              replace = FALSE)
Celma_data_balanced <- Celma_data_adduct %>%
  bind_rows(Celma_data_No_adduct)

#Picache data----

Picache_data = read_delim("Picache_adducts_CID.csv",
                          delim = ",",
                          col_names = TRUE)

Picache_data = Picache_data %>%
  filter(Adduct == "[M+Na]" | Adduct == "[M+H]") %>%
  select(PubChemCID, Adduct, CCS) %>%
  na.omit() %>%
  group_by(PubChemCID, Adduct) %>%
  summarise(CCS = mean(CCS)) %>%
  ungroup()

Picache_data = Picache_data %>%
  spread(key = Adduct, value = CCS)

Picache_data = Picache_data %>%
  mutate(M_Na = case_when(
    is.na(`[M+Na]`) ~ 0,
    TRUE ~ 1
  ))

Picache_data = Picache_data %>%
  select(PubChemCID, M_Na) %>%
  mutate(Lab = "Picache")


Picache_data  %>% 
  group_by(M_Na) %>%
  summarise(n())
#this is unbalance agin. 151 does not form adducts while 496 does.

Picache_data_adduct <- sample_n(Picache_data %>%
                                filter(M_Na == 1),
                              size = 150,
                              replace = FALSE)

Picache_data_No_adduct <- sample_n(Picache_data %>%
                                   filter(M_Na == 0),
                                 size = 150,
                                 replace = FALSE)

Picache_data_balanced <- Picache_data_adduct %>%
  bind_rows(Picache_data_No_adduct)

#binding all together----

data <- UT_data_balanced %>%
  bind_rows(Corey_data_balanced) %>%
  bind_rows(SU_data_balanced) %>%
  bind_rows(Celma_data_balanced) %>%
  bind_rows(Picache_data_balanced)


#Fingerprints----

fingerprints <- read_delim("PubChem_fingerprints_Bond_properties.csv",
                           delim = ",",
                           col_names = TRUE) %>%
  select(-SMILES, -X1, -SDFlineStart, -SDFlineEnd, -SDFID, -Fingerprint) %>%
  select(everything(), PubChemCID)

fingerprints <- fingerprints %>%
  select(-caret::nearZeroVar(fingerprints)) #removing fingerprints that do not change significantely between samples

fingerprints_with_NA <- fingerprints %>%
  select(everything()) %>%
  summarise_all(funs(sum(is.na(.))))
#only few columns have na-s and only for 5 rows, so I decide to drop the rows.

fingerprints <- fingerprints %>%
  na.omit()

#Putting adduct data and fingerprints together----
data <- data %>%
  select(-Name, -Class, -M_H, -SlopeH, -SlopeNa) %>%
  na.omit() %>%
  left_join(fingerprints)%>%
  na.omit()

write_delim(data,
            "balanced_adduct_data_training.csv",
            delim = ",")

compounds = data %>%
  select(PubChemCID) %>%
  unique()

set.seed(987123)
split = sample.split(compounds$PubChemCID, SplitRatio = 0.8)

compounds = compounds %>%
  mutate(split = split)

data_train = data %>%
  left_join(compounds) %>%
  filter(split == TRUE) %>%
  mutate(M_Na = as.factor(M_Na)) %>%
  select(-split)

#it turns out that for some labs there are still features with zero variance in and these need to be removed
data_train = data_train %>%
  select(-PubChemCID, -Lab, -M_Na, everything())
for (lab in levels(as.factor(data_train$Lab))) {
  data_train <- data_train %>%
    select(-caret::nearZeroVar(data_train %>%
                                 filter(Lab == lab) %>%
                                 select(-PubChemCID, -Lab, -M_Na)))
}



folds <- groupKFold(data_train$Lab, k = 5) 

fitControl <- trainControl(method = "repeatedcv",
                           number = 5,
                           repeats = 1,
                           index = folds)

classifier_svmpoly <- train(M_Na ~ .,
                    data = data_train %>%
                      select(-PubChemCID, -Lab),
                    method = "svmPoly",
                    trControl = fitControl)
saveRDS(classifier_svmpoly,
        file = "classifier_svmpoly.rds")

classifier_svmlinear <- train(M_Na ~ .,
                            data = data_train %>%
                              select(-PubChemCID, -Lab),
                            method = "svmLinear",
                            trControl = fitControl)

saveRDS(classifier_svmpoly,
        file = "classifier_svmpoly.rds")

classifier_RRF <- train(M_Na ~ .,
                            data = data_train %>%
                              select(-PubChemCID, -Lab),
                            method = "RRF",
                            trControl = fitControl)
saveRDS(classifier_RRF,
        file = "classifier_RRF.rds")

classifier_ada <- train(M_Na ~ .,
                        data = data_train %>%
                          select(-PubChemCID, -Lab),
                        method = "ada",
                        trControl = fitControl)

saveRDS(classifier_ada,
        file = "classifier_ada.rds")

classifier_DT <- train(M_Na ~ .,
                        data = data_train %>%
                          select(-PubChemCID, -Lab),
                        method = "rpart",
                        trControl = fitControl)

saveRDS(classifier_DT,
        file = "classifier_DT.rds")

classifier_knn <- train(M_Na ~ .,
                        data = data_train %>%
                          select(-PubChemCID, -Lab),
                        method = "knn",
                        trControl = fitControl)

saveRDS(classifier_knn,
        file = "classifier_knn.rds")

classifier_naive_bayes <- train(M_Na ~ .,
                        data = data_train %>%
                          select(-PubChemCID, -Lab),
                        method = "naive_bayes",
                        trControl = fitControl)

saveRDS(classifier_naive_bayes,
        file = "classifier_naive_bayes.rds")

data_pred <- data %>%
  mutate(M_Na_pred_svmLinear = predict(classifier_svmlinear, newdata = data),
         M_Na_pred_svmPoly = predict(classifier_svmpoly, newdata = data),
         M_Na_pred_RRF = predict(classifier_RRF, newdata = data),
         M_Na_pred_ada = predict(classifier_ada, newdata = data),
         M_Na_pred_DT = predict(classifier_DT, newdata = data),
         M_Na_pred_knn = predict(classifier_knn, newdata = data),
         M_Na_pred_naive_bayes = predict(classifier_naive_bayes, newdata = data)) %>%
  left_join(compounds)

data_test_pred<- data_pred %>%
  filter(split == FALSE) %>%
  select(PubChemCID, M_Na, M_Na_pred_svmLinear, M_Na_pred_svmPoly, M_Na_pred_RRF, M_Na_pred_ada, M_Na_pred_DT, M_Na_pred_knn, M_Na_pred_naive_bayes, everything())

ggplot(data = data_test_pred) +
  geom_point(mapping = aes(x = M_Na,
                           y = M_Na_pred_ada),
             size = 5,
             alpha = 1/50) +
  facet_wrap(~Lab) +
  theme_bw()

for (lab in levels(factor(data_pred$Lab))) {
  data_this <- data_test_pred %>%
    filter(Lab == lab)
  print(lab)
  print(
    mean(data_this$M_Na == data_this$M_Na_pred_naive_bayes)
  )
}

for (lab in levels(factor(data_pred$Lab))) {
  data_this <- data_test_pred %>%
    filter(Lab == lab)
  print(lab)
  print(
    CrossTable(data_this$M_Na, data_this$M_Na_pred_ada)
  )
}

classifier_ada$finalModel
pairs()
ada::varplot(classifier_ada$finalModel,  type = c("none","scores"))

#true positive and false positive rate analysis
data_pred2 = data %>%
  left_join(compounds) %>%
  filter(split == FALSE) 

data_pred2 = data_pred2 %>%
  mutate(M_Na_pred_ada = predict(classifier_ada, newdata = data_pred2, type = "prob")[,2]) %>%
  left_join(compounds) %>% 
  filter(split == FALSE)

data_pred2 = data_pred2 %>%
  select(PubChemCID, Lab, M_Na, M_Na_pred_ada)

library(ROCR)
pred <- prediction(data_pred2$M_Na_pred_ada, data_pred2$M_Na)
perf <- performance(pred,"tpr","fpr")
plot(perf,colorize=TRUE)

summary <- tibble(Cutoff = pred@cutoffs[[1]], 
                  TruePositives = pred@tp[[1]],
                  FalsePositives = pred@fp[[1]],
                  TrueNegatives = pred@tn[[1]],
                  FalseNegatives = pred@fn[[1]]) %>%
  mutate(Cutoff = case_when(
    Cutoff == "Inf" ~ 1,
    TRUE ~ Cutoff
  )) %>%
  add_row(Cutoff = 0, 
          TruePositives = sum(data_pred2$M_Na == 1), 
          FalsePositives = sum(data_pred2$M_Na != 1), 
          TrueNegatives = 0, 
          FalseNegatives = 0)

summary <- summary %>%
  mutate(TruePositiveRate = TruePositives / (TruePositives + FalseNegatives),
         FalsePositiveRate = FalsePositives / (FalsePositives + TrueNegatives ))

cutoff_plot =
ggplot(data = summary) +
  geom_line(mapping = aes(x = FalsePositiveRate,
                          y = TruePositiveRate, 
                          text = Cutoff),
            color = basecolor) +
  my_theme

plotly::ggplotly(cutoff_plot)


cutoff_plot
ggsave("TPR_FPR.svg",
       width = 8,
       height = 8,
       units = "cm")

summary <- summary %>%
  mutate(cost_pos = 1*FalsePositives + FalseNegatives,
         cost_neg = FalsePositives + 3*FalseNegatives)

ggplot(data = summary ) +
  geom_line(mapping = aes(x = FalsePositiveRate,
                          y = TruePositiveRate))+
  geom_point(data = summary %>% slice(which.min(cost_pos)),
             mapping = aes(x = FalsePositiveRate,
                           y = TruePositiveRate),
             size = 3,
             color = "red") +
  geom_point(data = summary %>% slice(which.min(cost_neg)),
             mapping = aes(x = FalsePositiveRate,
                           y = TruePositiveRate),
             size = 3,
             color = "blue")



classifier_ada = readRDS("classifier_ada.rds")







#read in the validation set and calculate for this
#after rerunnign the read-in
data_all = UT_data %>%
  bind_rows(Corey_data) %>%
  bind_rows(SU_data) %>%
  bind_rows(Celma_data) %>%
  bind_rows(Picache_data)

data_all = data_all %>%
  left_join(fingerprints)

write_delim(data_all %>%
              select(-Name, -Class, -M_H, -SlopeH, -SlopeNa) %>%
              na.omit(),
            "all_adduct_daat_training.csv",
            delim = ",")


data_all_not_train_not_test = data_all %>%
  anti_join(data)


data_all_not_train_not_test = data_all_not_train_not_test %>%
  select(-Name, -Class, -M_H, -SlopeH, -SlopeNa) %>%
  na.omit() 

data_all_not_train_not_test = data_all_not_train_not_test %>%
  mutate(M_Na_pred_prop = predict(classifier_ada, newdata = data_all_not_train_not_test, type = "prob")[,2]) %>%
  mutate(M_Na_pred = case_when(
    M_Na_pred_prop > 0.508 ~ 1,
    TRUE ~ 0
  ))

ggplot(data = data_all_not_train_not_test) +
  geom_point(mapping = aes(x = M_Na,
                           y = M_Na_pred_prop),
             size = 5,
             alpha = 1/10) +
  facet_wrap(~Lab) +
  theme_bw()

for (lab in levels(factor(data_all_not_train_not_test$Lab))) {
  data_this <- data_all_not_train_not_test %>%
    filter(Lab == lab)
  print(lab)
  print(
    CrossTable(data_this$M_Na, data_this$M_Na_pred)
  )
}

for (lab in levels(factor(data_all_not_train_not_test$Lab))) {
  data_this <- data_all_not_train_not_test %>%
    filter(Lab == lab)
  print(lab)
  print(
    mean(data_this$M_Na == data_this$M_Na_pred)
  )
}


validation = read_delim("ac8b04567_si_002_CID.csv",
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

validation_fingerprints = read_delim("PubChem_fingerprints_bond_properties_validation.csv",
                                     delim = ",",
                                     col_names = TRUE) 

validation = validation %>%
  left_join(validation_fingerprints) %>%
  na.omit()


validation = validation %>%
  mutate(M_Na_pred_prop = predict(classifier_ada, newdata = validation, type = "prob")[,2]) %>%
  mutate(M_Na_pred = case_when(
    M_Na_pred_prop > 0.508 ~ 1,
    TRUE ~ 0
  ))  

ggplot(data = validation) +
  geom_point(mapping = aes(x = M_Na,
                           y = M_Na_pred_prop),
             size = 5,
             alpha = 1/20) +
  my_theme

val_plot =
ggplot(data = validation %>%
         group_by(ExactMass) %>%
         mutate(count = n(),
                prob_dif = max(M_Na_pred_prop) -min(M_Na_pred_prop)) %>%
         ungroup() %>%
         filter(count > 1 & prob_dif > 0.1)) +
  geom_point(mapping = aes(x = ExactMass,
                           y = M_Na_pred_prop,
                           text = Name,
                           color = M_Na),
             size = 2) +
  geom_abline(slope = 0, intercept = 0.508) +
  my_theme

plotly::ggplotly(val_plot)

mean(validation$M_Na == validation$M_Na_pred)
CrossTable(validation$M_Na, validation$M_Na_pred)
#false positive 0.368
64/174
#false neg 0.183
73/398
