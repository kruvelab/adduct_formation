library(tidyverse)
library(caret)
library(caTools)
library(gmodels)
library(ROCR)
library(yardstick)
setwd("C:/Users/annel/OneDrive - Kruvelab/Kruvelab/computational/IE mudeli script ja failid/GitHub/adduct_formation/code")
source('my_theme.R')

setwd("C:/Users/annel/OneDrive - Kruvelab/Kruvelab/computational/IE mudeli script ja failid/GitHub/adduct_formation/data/training")

#reading in pretrained models----
classifier_svmpoly = readRDS("classifier_svmpoly.rds")
classifier_svmlinear = readRDS("classifier_svmlinear.rds")
classifier_RRF = readRDS("classifier_RRF.rds")
classifier_ada = readRDS("classifier_ada.rds")
classifier_DT = readRDS("classifier_DT.rds")
classifier_knn = readRDS("classifier_knn.rds")
classifier_naive_bayes = readRDS("classifier_naive_bayes.rds")

#reading in test set and making predictions -------

data_test = read_delim("adduct_data_test.csv",
                       delim = ",")

data_test_pred <- data_test %>%
  mutate(M_Na_pred_svmLinear = predict(classifier_svmlinear, newdata = data_test),
         M_Na_pred_svmPoly = predict(classifier_svmpoly, newdata = data_test),
         M_Na_pred_RRF = predict(classifier_RRF, newdata = data_test),
         M_Na_pred_ada = predict(classifier_ada, newdata = data_test),
         M_Na_pred_DT = predict(classifier_DT, newdata = data_test),
         M_Na_pred_knn = predict(classifier_knn, newdata = data_test),
         M_Na_pred_naive_bayes = predict(classifier_naive_bayes, newdata = data_test)) %>%
  mutate(M_Na_pred_svmLinear = case_when(
    M_Na_pred_svmLinear == "yes" ~ 1,
    TRUE ~ 0
  ),
  M_Na_pred_svmPoly = case_when(
    M_Na_pred_svmPoly == "yes" ~ 1,
    TRUE ~ 0
  )) %>%
  select(PubChemCID, M_Na, M_Na_pred_svmLinear, M_Na_pred_svmPoly, M_Na_pred_RRF, M_Na_pred_ada, M_Na_pred_DT, M_Na_pred_knn, M_Na_pred_naive_bayes, everything())

write_delim(data_test_pred,
            "Predicted_adduct_formation_test_set.csv",
            delim = ",")

#model evaluation with different metrics-------

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

write_delim(balanced_accuracy,
            "balanced_accuracy.csv",
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


#true positive and false positive rate analysis-------

for (lab in levels(factor(data_test$Lab))) {
  data_test_pred2 = data_test %>%
    filter(Lab == lab) 
  
  data_test_pred2 = data_test_pred2 %>%
    mutate(M_Na_pred_svmPoly = predict(classifier_svmpoly, newdata = data_test_pred2, type = "prob")[,2]) %>%
    select(PubChemCID, Lab, M_Na, M_Na_pred_svmPoly)
  
  pred <- prediction(data_test_pred2$M_Na_pred_svmPoly, data_test_pred2$M_Na)
  perf <- performance(pred,"tpr","fpr")
  plot(perf, colorize=TRUE)
  
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
            TruePositives = sum(data_test_pred2$M_Na == 1), 
            FalsePositives = sum(data_test_pred2$M_Na != 1), 
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
  
  cutoff_plot
  ggsave(paste(lab, "_TPR_FPR.svg", sep = ""),
         width = 8,
         height = 8,
         units = "cm")
}
