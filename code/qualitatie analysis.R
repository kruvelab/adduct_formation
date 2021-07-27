

desc_names = read_delim("PubChem_descs_names.csv",
                        delim = ",",
                        col_names = TRUE) %>%
  rename(desc = Position)


SU_data2 = SU_data %>%
  left_join(fingerprints)

test_SU_data = SU_data2 %>%
  select(5, 8:305) %>%
  na.omit()

significant_SU = tibble()
for (variable in 2:length(test_SU_data)) {
  ct = CrossTable(test_SU_data$M_Na, unname(c(test_SU_data[variable]))[[1]])$t
  ct = matrix(ct,
              nrow = 2,
              dimnames = list(Guess = c("Milk", "Tea"),
                              Truth = c("Milk", "Tea")))
  
  ft = fisher.test(x = ct)
  p = ft$p.value
  if (p < 0.05) {
    significant_SU = significant_SU %>%
      bind_rows(tibble("desc" = colnames(test_SU_data[variable]),
                       "p_Costalunga" = p))
  }
}

UT_data2 = UT_data %>%
  left_join(fingerprints)


test_UT_data = UT_data2 %>%
  select(5, 7:310) %>%
  na.omit()


significant_UT = tibble()
for (variable in 2:length(test_UT_data)) {
  ct = CrossTable(test_UT_data$M_Na, unname(c(test_UT_data[variable]))[[1]])$t
  ct = matrix(ct,
              nrow = 2,
              dimnames = list(Guess = c("Milk", "Tea"),
                              Truth = c("Milk", "Tea")))
  
  ft = fisher.test(x = ct)
  p = ft$p.value
  if (p < 0.05) {
    significant_UT = significant_UT %>%
      bind_rows(tibble("desc" = colnames(test_UT_data[variable]),
                       "p_Liigand" = p))
  }
}


Corey_data2 = Corey_data %>%
  left_join(fingerprints)


test_Corey_data = Corey_data2 %>%
  select(4, 6:309) %>%
  na.omit()


significant_Corey = tibble()
for (variable in 2:length(test_Corey_data)) {
  ct = CrossTable(test_Corey_data$M_Na, unname(c(test_Corey_data[variable]))[[1]])$t
  ct = matrix(ct,
              nrow = 2,
              dimnames = list(Guess = c("Milk", "Tea"),
                              Truth = c("Milk", "Tea")))
  
  ft = fisher.test(x = ct)
  p = ft$p.value
  if (p < 0.05) {
    significant_Corey = significant_Corey %>%
      bind_rows(tibble("desc" = colnames(test_Corey_data[variable]),
                       "p_Broeckling" = p))
  }
}


Celma_data2 = Celma_data %>%
  left_join(fingerprints)


test_Celma_data = Celma_data2 %>%
  select(2, 4:307) %>%
  na.omit()


significant_Celma = tibble()
for (variable in 2:length(test_Celma_data)) {
  ct = CrossTable(test_Celma_data$M_Na, unname(c(test_Celma_data[variable]))[[1]])$t
  ct = matrix(ct,
              nrow = 2,
              dimnames = list(Guess = c("Milk", "Tea"),
                              Truth = c("Milk", "Tea")))
  
  ft = fisher.test(x = ct)
  p = ft$p.value
  if (p < 0.05) {
    significant_Celma = significant_Celma %>%
      bind_rows(tibble("desc" = colnames(test_Celma_data[variable]),
                       "p_Celma" = p))
  }
}


Picache_data2 = Picache_data %>%
  left_join(fingerprints)


test_Picache_data = Picache_data2 %>%
  select(2, 4:307) %>%
  na.omit()


significant_Picache = tibble()
for (variable in 2:length(test_Picache_data)) {
  ct = CrossTable(test_Picache_data$M_Na, unname(c(test_Picache_data[variable]))[[1]])$t
  ct = matrix(ct,
              nrow = 2,
              dimnames = list(Guess = c("Milk", "Tea"),
                              Truth = c("Milk", "Tea")))
  
  ft = fisher.test(x = ct)
  p = ft$p.value
  if (p < 0.05) {
    significant_Picache = significant_Picache %>%
      bind_rows(tibble("desc" = colnames(test_Picache_data[variable]),
                       "p_Picache" = p))
  }
}

significant_SU %>%
  inner_join(significant_UT) %>%
  inner_join(significant_Corey) %>% 
  inner_join(significant_Celma) %>%
  inner_join(significant_Picache) %>%
  left_join(desc_names)


summary_significance = significant_SU %>%
  full_join(significant_UT) %>%
  full_join(significant_Corey) %>% 
  full_join(significant_Celma) %>%
  full_join(significant_Picache) %>%
  left_join(desc_names)

write_delim(summary_significance,
            "summary_significant_descs.csv",
            delim = ",")



Collected_adduct = UT_data %>%
  select(PubChemCID, M_Na) %>%
  rename(M_Na_UT = M_Na) %>%
  left_join(SU_data %>%
              select(Name, PubChemCID, M_Na) %>%
              rename(M_Na_SU = M_Na)) %>%
  left_join(Corey_data %>%
              select(PubChemCID, M_Na) %>%
              rename(M_Na_Broeckling = M_Na))

Collected_adduct = Collected_adduct %>%
  rowwise() %>%
  mutate(count =3 - sum(is.na(c(M_Na_UT, M_Na_SU, M_Na_Broeckling)))) %>%
  filter(count > 1)

Collected_adduct = Collected_adduct %>%
  rowwise() %>%
  mutate(mean = mean(c(M_Na_UT, M_Na_SU, M_Na_Broeckling), na.rm = TRUE),
         sd = sd(c(M_Na_UT, M_Na_SU, M_Na_Broeckling), na.rm = TRUE))


test_SU_data2 = test_SU_data %>%
  gather(key = "desc", value = "value", 2:4)

ggplot(data = test_SU_data %>%
         filter(desc != "max_dist_2O" & 
                  desc != "max_dist_carbonyl" &
                  desc != "max_path_non_rotatabale_bonds_2O" &
                  desc != "max_path_non_rotatabale_bonds_carbonyl" & 
                  desc != "min_dist_2O" & 
                  desc != "min_dist_carbonyl" &
                  desc != "min_path_non_rotatabale_bonds_2O" &
                  desc != "min_path_non_rotatabale_bonds_carbonyl" )) +
  geom_histogram(mapping = aes(x = value, fill = factor(M_Na)),
                 binwidth = 0.5) +
  facet_wrap(~desc)

ggsave("PubChem_SU_descs2.svg",
       width = 16,
       units = "cm")



UT_data = UT_data %>%
  mutate(M_Na_UT = M_Na) %>%
  select(PubChemCID, M_Na_UT) %>%
  na.omit()

Corey_data = Corey_data %>%
  mutate(M_Na_Corey = M_Na) %>%
  select(PubChemCID, M_Na_Corey) %>%
  na.omit()

SU_data = SU_data %>%
  mutate(M_Na_SU = M_Na) %>%
  select(PubChemCID, M_Na_SU) %>%
  na.omit()

Celma_data = Celma_data %>%
  mutate(M_Na_Celma = M_Na) %>%
  select(PubChemCID, M_Na_Celma) %>%
  na.omit()

Picache_data = Picache_data %>%
  mutate(M_Na_Picache = M_Na) %>%
  select(PubChemCID, M_Na_Picache) %>%
  na.omit()

data = UT_data %>%
  full_join(Corey_data) %>%
  full_join(SU_data) %>%
  full_join(Celma_data) %>%
  full_join(Picache_data)

data = data %>%
  rowwise() %>%
  mutate(count = 5 - sum(is.na(c(M_Na_UT, M_Na_SU, M_Na_Corey, M_Na_Celma, M_Na_Picache))),
         mean = mean(na.omit(c(M_Na_UT, M_Na_SU, M_Na_Corey, M_Na_Celma, M_Na_Picache))),
         stdev = sd(na.omit(c(M_Na_UT, M_Na_SU, M_Na_Corey, M_Na_Celma, M_Na_Picache)))) %>%
  filter(count > 1)

#68% näitavad ideaalset kooskõla



