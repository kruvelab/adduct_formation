library(tidyverse)
library(caret)
library(gmodels)
library(caTools)
setwd("C:/Users/annel/OneDrive - Kruvelab/Kruvelab/computational/IE mudeli script ja failid/adduct_formation/code")
source('my_theme.R')

setwd("C:/Users/annel/OneDrive - Kruvelab/Kruvelab/computational/IE mudeli script ja failid/GitHub/adduct_formation/data/training")

#Adduct data ----
#UT data-----
UT_data = read_delim("UT_adducts_CID.csv",
                      delim = ",",
                      col_names = TRUE) %>%
  na.omit()

UT_data = UT_data %>%
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

set.seed(123)

#there are about 2x so (224 vs 127) many compounds not forming an adduct as there are ones forming an adduct
UT_data_adduct = sample_n(UT_data %>%
                            filter(M_Na == 1),
                          size = 100,
                          replace = FALSE)

UT_data_No_adduct = sample_n(UT_data %>%
                             filter(M_Na == 0),
                           size = 100,
                           replace = FALSE)
UT_data_balanced = UT_data_adduct %>%
  bind_rows(UT_data_No_adduct)

#Corey data -----

Corey_data = read_delim("Corey_adducts_CID.csv",
                      delim = ",",
                      col_names = TRUE)

Corey_data = Corey_data %>%
  mutate(M_Na = case_when(
                    Class == 1 ~ 1,
                    Class == 2 ~ 1,
                    TRUE ~ 0),
         Lab = "Corey")

Corey_data  %>% 
  group_by(M_Na) %>%
  summarise(n())

#there are about 6x so (72 vs 525) many compounds forming an adduct as there are ones forming an adduct
Corey_data_No_adduct = sample_n(Corey_data %>%
                             filter(M_Na == 0),
                           size = 100,
                           replace = TRUE)

Corey_data_adduct = sample_n(Corey_data %>%
                                   filter(M_Na == 1),
                                 size = 100,
                                 replace = FALSE)

Corey_data_balanced = Corey_data_adduct %>%
  bind_rows(Corey_data_No_adduct)

#SU data -----

SU_data = read_delim("SU_adducts_CID.csv",
                      delim = ",",
                      col_names = TRUE)

SU_data = SU_data %>%
  mutate(M_Na = case_when(
                  SlopeNa == 0 ~ 0,
                  TRUE ~ 1),
        M_H = case_when(
                  SlopeH == 0 ~ 0,
                  TRUE ~ 1),
        Lab = "SU")

SU_data  %>% 
  group_by(M_Na, M_H) %>%
  summarise(n())

#this is an almost even dataset (53 vs 41)

SU_data_adduct = sample_n(SU_data %>%
                                filter(M_Na == 1),
                              size = 50,
                              replace = TRUE)

SU_data_No_adduct = sample_n(SU_data %>%
                                filter(M_Na == 0),
                              size = 50,
                              replace = TRUE)
SU_data_balanced = SU_data_adduct %>%
  bind_rows(SU_data_No_adduct)

#Celma data----

Celma_data = read_delim("Celma_adduct_CID.csv",
                              delim = ",",
                              col_names = TRUE)

Celma_data  %>% 
  group_by(M_Na) %>%
  summarise(n())
#this is an almost even dataset (274 vs 247)

Celma_data_adduct = sample_n(Celma_data %>%
                             filter(M_Na == 1),
                           size = 200,
                           replace = FALSE)

Celma_data_No_adduct = sample_n(Celma_data %>%
                                filter(M_Na == 0),
                              size = 200,
                              replace = FALSE)
Celma_data_balanced = Celma_data_adduct %>%
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
#this is unbalance agin. 151 does not form adducts while 469 does.

Picache_data_adduct = sample_n(Picache_data %>%
                                filter(M_Na == 1),
                              size = 100,
                              replace = FALSE)

Picache_data_No_adduct = sample_n(Picache_data %>%
                                   filter(M_Na == 0),
                                 size = 100,
                                 replace = FALSE)

Picache_data_balanced = Picache_data_adduct %>%
  bind_rows(Picache_data_No_adduct)

#binding all together----

data_all = UT_data %>%
  bind_rows(Corey_data) %>%
  bind_rows(SU_data) %>%
  bind_rows(Celma_data) %>%
  bind_rows(Picache_data)

data = UT_data_balanced %>%
  bind_rows(Corey_data_balanced) %>%
  bind_rows(SU_data_balanced) %>%
  bind_rows(Celma_data_balanced) %>%
  bind_rows(Picache_data_balanced) %>%
  select(-Name, -Class, -M_H, -SlopeH, -SlopeNa)

data_test = data_all %>%
  anti_join(data)

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
  na.omit() %>%
  unique()

#Putting adduct data and fingerprints together----
data <- data %>%
  na.omit() %>%
  left_join(fingerprints)%>%
  na.omit()

write_delim(data,
            "balanced_adduct_data_training.csv",
            delim = ",")

data = read_delim("balanced_adduct_data_training.csv",
                  delim = ",",
                  col_names = TRUE)



compounds = data %>%
  select(PubChemCID) %>%
  unique()

set.seed(987123)

data_train = data %>%
  mutate(M_Na = case_when(
    M_Na == 1 ~ "yes",
    TRUE ~ "no"
  )) %>%
  mutate(M_Na = factor(M_Na,
                          levels = c("no", "yes")))

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
                           index = folds,
                           classProbs = TRUE)

set.seed(987123)

classifier_svmpoly <- train(M_Na ~ .,
                    data = data_train %>%
                      select(-PubChemCID, -Lab),
                    method = "svmPoly",
                    trControl = fitControl,
                    metric = "ROC")

saveRDS(classifier_svmpoly,
        file = "classifier_svmpoly.rds")

classifier_svmpoly = readRDS("classifier_svmpoly.rds")

classifier_svmpoly
MLmetrics::F1_Score(data_train$M_Na, predict(classifier_svmpoly, newdata = data_train))
mltools::mcc(predict(classifier_svmpoly, newdata = data_train), data_train$M_Na)

roc_Imp = filterVarImp(x = data_train[, -ncol(data_train)], y = data_train$M_Na)

set.seed(987123)

classifier_svmlinear <- train(M_Na ~ .,
                            data = data_train %>%
                              select(-PubChemCID, -Lab),
                            method = "svmLinear",
                            trControl = fitControl)

saveRDS(classifier_svmlinear,
        file = "classifier_svmlinear.rds")

classifier_svmlinear = readRDS("classifier_svmlinear.rds")

classifier_svmlinear
MLmetrics::F1_Score(data_train$M_Na, predict(classifier_svmlinear, newdata = data_train))
mltools::mcc(predict(classifier_svmlinear, newdata = data_train), data_train$M_Na)



#data format for RRF and others-----
data_train = data %>%
  mutate(M_Na = factor(M_Na)) 

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

set.seed(987123)

classifier_RRF <- train(M_Na ~ .,
                            data = data_train %>%
                              select(-PubChemCID, -Lab),
                            method = "RRF",
                            trControl = fitControl)
saveRDS(classifier_RRF,
        file = "classifier_RRF.rds")

classifier_RRF = readRDS("classifier_RRF.rds")

classifier_RRF
MLmetrics::F1_Score(data_train$M_Na, predict(classifier_RRF, newdata = data_train))
mltools::mcc(predict(classifier_RRF, newdata = data_train), data_train$M_Na)

set.seed(987123)

classifier_ada <- train(M_Na ~ .,
                        data = data_train %>%
                          select(-PubChemCID, -Lab),
                        method = "ada",
                        trControl = fitControl)

saveRDS(classifier_ada,
        file = "classifier_ada.rds")

classifier_ada = readRDS("classifier_ada.rds")

classifier_ada
MLmetrics::F1_Score(data_train$M_Na, predict(classifier_ada, newdata = data_train))
mltools::mcc(predict(classifier_ada, newdata = data_train), data_train$M_Na)

set.seed(987123)

classifier_DT <- train(M_Na ~ .,
                        data = data_train %>%
                          select(-PubChemCID, -Lab),
                        method = "rpart",
                        trControl = fitControl)

saveRDS(classifier_DT,
        file = "classifier_DT.rds")

classifier_DT = readRDS("classifier_DT.rds")

classifier_DT
MLmetrics::F1_Score(data_train$M_Na, predict(classifier_DT, newdata = data_train))
mltools::mcc(predict(classifier_DT, newdata = data_train), data_train$M_Na)

set.seed(987123)

classifier_knn <- train(M_Na ~ .,
                        data = data_train %>%
                          select(-PubChemCID, -Lab),
                        method = "knn",
                        trControl = fitControl)

saveRDS(classifier_knn,
        file = "classifier_knn.rds")

classifier_knn = readRDS("classifier_knn.rds")

classifier_knn
MLmetrics::F1_Score(data_train$M_Na, predict(classifier_knn, newdata = data_train))
mltools::mcc(predict(classifier_knn, newdata = data_train), data_train$M_Na)

set.seed(987123)

classifier_naive_bayes <- train(M_Na ~ .,
                        data = data_train %>%
                          select(-PubChemCID, -Lab),
                        method = "naive_bayes",
                        trControl = fitControl)

saveRDS(classifier_naive_bayes,
        file = "classifier_naive_bayes.rds")

classifier_naive_bayes = readRDS("classifier_naive_bayes.rds")

classifier_naive_bayes
MLmetrics::F1_Score(data_train$M_Na, predict(classifier_naive_bayes, newdata = data_train))
mltools::mcc(predict(classifier_naive_bayes, newdata = data_train), data_train$M_Na)

data_train_pred <- data %>%
  mutate(M_Na_pred_svmLinear = predict(classifier_svmlinear, newdata = data),
         M_Na_pred_svmPoly = predict(classifier_svmpoly, newdata = data),
         M_Na_pred_RRF = predict(classifier_RRF, newdata = data),
         M_Na_pred_ada = predict(classifier_ada, newdata = data),
         M_Na_pred_DT = predict(classifier_DT, newdata = data),
         M_Na_pred_knn = predict(classifier_knn, newdata = data),
         M_Na_pred_naive_bayes = predict(classifier_naive_bayes, newdata = data))

data_train_pred = data_train_pred %>%
  mutate(M_Na_pred_svmLinear = case_when(
    M_Na_pred_svmLinear == "yes" ~ 1,
    TRUE ~ 0
  ),
  M_Na_pred_svmPoly = case_when(
    M_Na_pred_svmPoly == "yes" ~ 1,
    TRUE ~ 0
  ))

write_delim(data_train_pred,
            "Predicted_adduct_formation_training_set.csv",
            delim = ",")

data_test_pred <- data_test %>%
  select(PubChemCID, Lab, M_Na) %>%
  unique() %>%
  left_join(fingerprints) %>%
  na.omit() %>%
  unique()


data_test_pred <- data_test_pred %>%
  mutate(M_Na_pred_svmLinear = predict(classifier_svmlinear, newdata = data_test_pred),
         M_Na_pred_svmPoly = predict(classifier_svmpoly, newdata = data_test_pred),
         M_Na_pred_RRF = predict(classifier_RRF, newdata = data_test_pred),
         M_Na_pred_ada = predict(classifier_ada, newdata = data_test_pred),
         M_Na_pred_DT = predict(classifier_DT, newdata = data_test_pred),
         M_Na_pred_knn = predict(classifier_knn, newdata = data_test_pred),
         M_Na_pred_naive_bayes = predict(classifier_naive_bayes, newdata = data_test_pred))

data_test_pred = data_test_pred %>%
  mutate(M_Na_pred_svmLinear = case_when(
    M_Na_pred_svmLinear == "yes" ~ 1,
    TRUE ~ 0
  ),
  M_Na_pred_svmPoly = case_when(
    M_Na_pred_svmPoly == "yes" ~ 1,
    TRUE ~ 0
  ))


data_test_pred<- data_test_pred %>%
  select(PubChemCID, M_Na, M_Na_pred_svmLinear, M_Na_pred_svmPoly, M_Na_pred_RRF, M_Na_pred_ada, M_Na_pred_DT, M_Na_pred_knn, M_Na_pred_naive_bayes, everything())

write_delim(data_test_pred,
            "Predicted_adduct_formation_test_set.csv",
            delim = ",")


ggplot(data = data_test_pred) +
  geom_point(mapping = aes(x = M_Na,
                           y = M_Na_pred_ada),
             size = 5,
             alpha = 1/50) +
  facet_wrap(~Lab) +
  theme_bw()

for (lab in levels(factor(data_test_pred$Lab))) {
  this_lab = data_test_pred %>%
    filter(Lab == lab)
  print(lab)
  CrossTable(this_lab$M_Na, this_lab$M_Na_pred_svmPoly)
}



accuracy = data_test_pred %>%
  group_by(Lab) %>%
  summarise(count_pos = sum(M_Na),
            count_neg = n() - sum(M_Na),
            accuracy_svmLinear = mean(M_Na == M_Na_pred_svmLinear),
            accuracy_svmPoly = mean(M_Na == M_Na_pred_svmPoly),
            accuracy_RRF = mean(M_Na == M_Na_pred_RRF),
            accuracy_ada = mean(M_Na == M_Na_pred_ada),
            accuracy_DT = mean(M_Na == M_Na_pred_DT),
            accuracy_knn = mean(M_Na == M_Na_pred_knn),
            accuracy_naive_bayes = mean(M_Na == M_Na_pred_naive_bayes)) %>%
  ungroup() %>%
  bind_rows(data_test_pred %>%
              group_by() %>%
              summarise(count_pos = sum(M_Na),
                        count_neg = n() - sum(M_Na),
                        accuracy_svmLinear = mean(M_Na == M_Na_pred_svmLinear),
                        accuracy_svmPoly = mean(M_Na == M_Na_pred_svmPoly),
                        accuracy_RRF = mean(M_Na == M_Na_pred_RRF),
                        accuracy_ada = mean(M_Na == M_Na_pred_ada),
                        accuracy_DT = mean(M_Na == M_Na_pred_DT),
                        accuracy_knn = mean(M_Na == M_Na_pred_knn),
                        accuracy_naive_bayes = mean(M_Na == M_Na_pred_naive_bayes)) %>%
              ungroup())

balanced_accuracy = data_test_pred %>%
  mutate(M_Na = as.factor(M_Na),
         M_Na_pred_svmLinear = as.factor(M_Na_pred_svmLinear),
         M_Na_pred_svmPoly = as.factor(M_Na_pred_svmPoly)) %>%
  group_by(Lab) %>%
  summarise(accuracy_svmLinear = bal_accuracy_vec(M_Na, M_Na_pred_svmLinear),
            accuracy_svmPoly = bal_accuracy_vec(M_Na, M_Na_pred_svmPoly),
            accuracy_RRF = bal_accuracy_vec(M_Na, M_Na_pred_RRF),
            accuracy_ada = bal_accuracy_vec(M_Na, M_Na_pred_ada),
            accuracy_DT = bal_accuracy_vec(M_Na, M_Na_pred_DT),
            accuracy_knn = bal_accuracy_vec(M_Na, M_Na_pred_knn),
            accuracy_naive_bayes = bal_accuracy_vec(M_Na, M_Na_pred_naive_bayes)) %>%
  ungroup() %>%
  bind_rows(data_test_pred %>%
              mutate(M_Na = as.factor(M_Na),
                     M_Na_pred_svmLinear = as.factor(M_Na_pred_svmLinear),
                     M_Na_pred_svmPoly = as.factor(M_Na_pred_svmPoly)) %>%
              group_by() %>%
              summarise(accuracy_svmLinear = bal_accuracy_vec(M_Na, M_Na_pred_svmLinear),
                        accuracy_svmPoly = bal_accuracy_vec(M_Na, M_Na_pred_svmPoly),
                        accuracy_RRF = bal_accuracy_vec(M_Na, M_Na_pred_RRF),
                        accuracy_ada = bal_accuracy_vec(M_Na, M_Na_pred_ada),
                        accuracy_DT = bal_accuracy_vec(M_Na, M_Na_pred_DT),
                        accuracy_knn = bal_accuracy_vec(M_Na, M_Na_pred_knn),
                        accuracy_naive_bayes = bal_accuracy_vec(M_Na, M_Na_pred_naive_bayes)) %>%
              ungroup())

write_delim(balanced_accuracy,
            "balanced_accuracy.csv",
            delim = ",")

data_test_pred2 = data_test_pred %>%
  mutate(M_Na = case_when(
    M_Na == 1 ~ 0,
    TRUE ~ 1),
    M_Na_pred_svmLinear = case_when(
      M_Na_pred_svmLinear == 1 ~ 0,
      TRUE ~ 1),
    M_Na_pred_svmPoly = case_when(
      M_Na_pred_svmPoly == 1 ~ 0,
      TRUE ~ 1),
    M_Na_pred_RRF = case_when(
      M_Na_pred_RRF == 1 ~ 0,
      TRUE ~ 1),
    M_Na_pred_ada = case_when(
      M_Na_pred_ada == 1 ~ 0,
      TRUE ~ 1),
    M_Na_pred_DT = case_when(
      M_Na_pred_DT == 1 ~ 0,
      TRUE ~ 1),
    M_Na_pred_knn = case_when(
      M_Na_pred_knn == 1 ~ 0,
      TRUE ~ 1),
    M_Na_pred_naive_bayes = case_when(
      M_Na_pred_naive_bayes == 1 ~ 0,
      TRUE ~ 1),)


F1_scores = data_test_pred2 %>%
  group_by(Lab) %>%
  summarise(count_pos = sum(M_Na),
            count_neg = n() - sum(M_Na),
            F1_svmLinear = MLmetrics::F1_Score(y_true = M_Na, y_pred = M_Na_pred_svmLinear), 
            F1_svmPoly = MLmetrics::F1_Score(y_true = M_Na, y_pred = M_Na_pred_svmPoly),
            F1_RRF = MLmetrics::F1_Score(y_true = M_Na, y_pred = M_Na_pred_RRF),
            F1_ada = MLmetrics::F1_Score(y_true = M_Na, y_pred = M_Na_pred_ada),
            F1_DT = MLmetrics::F1_Score(y_true = M_Na, y_pred = M_Na_pred_DT),
            F1_knn = MLmetrics::F1_Score(y_true = M_Na, y_pred = M_Na_pred_knn),
            F1_naive_bayes = MLmetrics::F1_Score(y_true = M_Na, y_pred = M_Na_pred_naive_bayes)) %>%
  ungroup() %>%
  bind_rows(data_test_pred2 %>%
              group_by() %>%
              summarise(count_pos = sum(M_Na),
                        count_neg = n() - sum(M_Na),
                        F1_svmLinear = MLmetrics::F1_Score(M_Na, M_Na_pred_svmLinear),
                        F1_svmPoly = MLmetrics::F1_Score(M_Na, M_Na_pred_svmPoly),
                        F1_RRF = MLmetrics::F1_Score(M_Na, M_Na_pred_RRF),
                        F1_ada = MLmetrics::F1_Score(M_Na, M_Na_pred_ada),
                        F1_DT = MLmetrics::F1_Score(M_Na, M_Na_pred_DT),
                        F1_knn = MLmetrics::F1_Score(M_Na, M_Na_pred_knn),
                        F1_naive_bayes = MLmetrics::F1_Score(M_Na, M_Na_pred_naive_bayes)) %>%
              ungroup())

Matthews = data_test_pred2 %>%
  group_by(Lab) %>%
  summarise(count_pos = sum(M_Na),
            count_neg = n() - sum(M_Na),
            Matthews_svmLinear = mccr::mccr(M_Na, M_Na_pred_svmLinear),
            Matthews_svmPoly = mccr::mccr(M_Na, M_Na_pred_svmPoly),
            Matthews_RRF = mccr::mccr(M_Na, M_Na_pred_RRF),
            Matthews_ada = mccr::mccr(M_Na, M_Na_pred_ada),
            Matthews_DT = mccr::mccr(M_Na, M_Na_pred_DT),
            Matthews_knn = mccr::mccr(M_Na, M_Na_pred_knn),
            Matthews_naive_bayes = mccr::mccr(M_Na, M_Na_pred_naive_bayes)) %>%
  ungroup() %>%
  bind_rows(data_test_pred2 %>%
              group_by() %>%
              summarise(count_pos = sum(M_Na),
                        count_neg = n() - sum(M_Na),
                        Matthews_svmLinear = mccr::mccr(M_Na, M_Na_pred_svmLinear),
                        Matthews_svmPoly = mccr::mccr(M_Na, M_Na_pred_svmPoly),
                        Matthews_RRF = mccr::mccr(M_Na, M_Na_pred_RRF),
                        Matthews_ada = mccr::mccr(M_Na, M_Na_pred_ada),
                        Matthews_DT = mccr::mccr(M_Na, M_Na_pred_DT),
                        Matthews_knn = mccr::mccr(M_Na, M_Na_pred_knn),
                        Matthews_naive_bayes = mccr::mccr(M_Na, M_Na_pred_naive_bayes)) %>%
              ungroup())

Recall = data_test_pred2 %>%
  group_by(Lab) %>%
  summarise(count_pos = sum(M_Na),
            count_neg = n() - sum(M_Na),
            recall_svmLinear = caret::recall(factor(M_Na_pred_svmLinear), factor(M_Na)),
            recall_svmPoly = caret::recall(factor(M_Na_pred_svmPoly), factor(M_Na)),
            recall_RRF = caret::recall(factor(M_Na_pred_RRF), factor(M_Na)),
            recall_ada = caret::recall(factor(M_Na_pred_ada), factor(M_Na)),
            recall_DT = caret::recall(factor(M_Na_pred_DT), factor(M_Na)),
            recall_knn = caret::recall(factor(M_Na_pred_knn), factor(M_Na)),
            recall_naive_bayes = caret::recall(factor(M_Na_pred_naive_bayes), factor(M_Na))) %>%
  ungroup() %>%
  bind_rows(data_test_pred2 %>%
              group_by() %>%
              summarise(count_pos = sum(M_Na),
                        count_neg = n() - sum(M_Na),
                        recall_svmLinear = caret::recall(factor(M_Na_pred_svmLinear), factor(M_Na)),
                        recall_svmPoly = caret::recall(factor(M_Na_pred_svmPoly), factor(M_Na)),
                        recall_RRF = caret::recall(factor(M_Na_pred_RRF), factor(M_Na)),
                        recall_ada = caret::recall(factor(M_Na_pred_ada), factor(M_Na)),
                        recall_DT = caret::recall(factor(M_Na_pred_DT), factor(M_Na)),
                        recall_knn = caret::recall(factor(M_Na_pred_knn), factor(M_Na)),
                        recall_naive_bayes = caret::recall(factor(M_Na_pred_naive_bayes), factor(M_Na))) %>%
              ungroup() )

Precision = data_test_pred2 %>%
  group_by(Lab) %>%
  summarise(count_pos = sum(M_Na),
            count_neg = n() - sum(M_Na),
            precision_svmLinear = caret::precision(factor(M_Na_pred_svmLinear), factor(M_Na)),
            precision_svmPoly = caret::precision(factor(M_Na_pred_svmPoly), factor(M_Na)),
            precision_RRF = caret::precision(factor(M_Na_pred_RRF), factor(M_Na)),
            precision_ada = caret::precision(factor(M_Na_pred_ada), factor(M_Na)),
            precision_DT = caret::precision(factor(M_Na_pred_DT), factor(M_Na)),
            precision_knn = caret::precision(factor(M_Na_pred_knn), factor(M_Na)),
            precision_naive_bayes = caret::precision(factor(M_Na_pred_naive_bayes), factor(M_Na))) %>%
  ungroup() %>%
  bind_rows(data_test_pred2 %>%
              group_by() %>%
              summarise(count_pos = sum(M_Na),
                        count_neg = n() - sum(M_Na),
                        precision_svmLinear = caret::precision(factor(M_Na_pred_svmLinear), factor(M_Na)),
                        precision_svmPoly = caret::precision(factor(M_Na_pred_svmPoly), factor(M_Na)),
                        precision_RRF = caret::precision(factor(M_Na_pred_RRF), factor(M_Na)),
                        precision_ada = caret::precision(factor(M_Na_pred_ada), factor(M_Na)),
                        precision_DT = caret::precision(factor(M_Na_pred_DT), factor(M_Na)),
                        precision_knn = caret::precision(factor(M_Na_pred_knn), factor(M_Na)),
                        precision_naive_bayes = caret::precision(factor(M_Na_pred_naive_bayes), factor(M_Na))) %>%
              ungroup())

write_delim(accuracy,
            "accuracy_test_set.csv",
            delim = ",")

write_delim(F1_scores,
            "F1_score_test_set.csv",
            delim = ",")

write_delim(Matthews,
            "Matthews_correlation_test_set.csv",
            delim = ",")

write_delim(Recall,
            "Recall_test_set.csv",
            delim = ",")

write_delim(Precision,
            "Precision_test_set.csv",
            delim = ",")

classifier_svmpoly$finalModel
caret::varImp(classifier_RRF, scale = FALSE)


#true positive and false positive rate analysis
data_pred2 = data_test %>%
  select(PubChemCID, M_Na, Lab) %>%
  left_join(fingerprints) %>%
  na.omit()

data_pred2 = data_pred2 %>%
  mutate(M_Na_pred_svmPoly = predict(classifier_svmpoly, newdata = data_pred2, type = "prob")[,2]) 

data_pred2 = data_pred2 %>%
  select(PubChemCID, Lab, M_Na, M_Na_pred_svmPoly)

library(ROCR)
pred <- prediction(data_pred2$M_Na_pred_svmPoly, data_pred2$M_Na)
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



classifier_svmpoly = readRDS("classifier_svmpoy.rds")


#Validation------

setwd("C:/Users/annel/OneDrive - Kruvelab/Kruvelab/computational/IE mudeli script ja failid/adduct_formation/data/validation")
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
  mutate(M_Na_pred_prob = predict(classifier_svmpoly, newdata = validation, type = "prob")[,2]) %>%
  mutate(M_Na_pred = case_when(
    M_Na_pred_prob > 0.5 ~ 1,
    TRUE ~ 0
  ))

ggplot(data = validation) +
  geom_point(mapping = aes(x = M_Na,
                           y = M_Na_pred_prob),
             size = 5,
             alpha = 1/20) +
  my_theme

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
library(mccr)
mccr::mccr(validation2$M_Na_pred, validation2$M_Na)
#0.389