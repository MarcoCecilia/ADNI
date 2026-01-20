
setwd("C:/Users/HP/OneDrive/Documenti/APPLIED STATISTICS/PROGETTO")
# ---- MANOVA ----
library(dplyr)
library(caret)

library(readxl)
df <- read_excel("Dataset_Alz_feat_eng.xlsx")   # di default legge il primo foglio

names(df)

# Estrai solo l'ultima visita per ogni soggetto
df_last <- df %>%
  group_by(Subject_ID) %>%
  filter(Visit == max(Visit)) %>%
  ungroup()

# Verifica bilanciamento delle classi originali
cat("\n Frequenze DIAGNOSIS:\n")
print(table(df_last$DIAGNOSIS))

# MANOVA con volumi cerebrali
manova_model <- manova(cbind(Hippocampus_Total, Thalamus_Total, SuperiorTemporal_Total,
                             Insula_Total, Precuneus_Total, Mammilare_Total,
                             CaudalMiddleFrontal_Total, InfLatVentricle_Total) 
                       ~ DIAGNOSIS, data = df_last)

# Test globale MANOVA
cat("\n MANOVA - Test globale (Wilks):\n")
print(summary(manova_model, test = "Wilks"))

# ANOVA univariata per ciascuna variabile
cat("\n ANOVA per singola variabile:\n")
print(summary.aov(manova_model))

#---- logit---

# LIBRERIE
library(dplyr)
library(caret)

# Seleziona ultima visita
df_last <- df %>%
  group_by(Subject_ID) %>%
  filter(Visit == max(Visit)) %>%
  ungroup()

# Crea variabile binaria: AD vs CN/MCI
df_last$DIAGNOSIS_BIN <- factor(
  ifelse(df_last$DIAGNOSIS == 3, "AD", "NoAD"),
  levels = c("NoAD", "AD")
)

# Variabili selezionate
selected_vars <- c("DIAGNOSIS_BIN",
                   "Hippocampus_Total", "Precuneus_Total", "InfLatVentricle_Total",
                   "genotype", 
                   "rs744373_C", "rs11767557_C", "rs11771145_A", "rs11136000_T",
                   "rs3851179_A", "rs17125944_C", "rs3764650_G")

df_model <- df_last[, selected_vars] %>%
  na.omit()

# Mantieni le variabili categoriche come fattori
categorical_vars <- c("genotype", "rs744373_C", "rs11767557_C", "rs11771145_A", 
                      "rs11136000_T", "rs3851179_A", "rs17125944_C", "rs3764650_G")

df_model[categorical_vars] <- lapply(df_model[categorical_vars], factor)

# REGRESSIONE LOGISTICA
logit_model <- glm(DIAGNOSIS_BIN ~ ., data = df_model, family = binomial)

# Risultati
summary(logit_model)

# Odds ratio
cat("\n Odds Ratio:\n")
print(round(exp(coef(logit_model)), 3))




#---- logit genetico: AD vs NonAD-----
# LIBRERIE
library(dplyr)
library(car)

# DATI E PREPROCESSING
df_last <- df %>%
  group_by(Subject_ID) %>%
  filter(Visit == max(Visit)) %>%
  ungroup()

selected_vars <- c("DIAGNOSIS", "genotype", 
                   "rs744373_C", "rs11767557_C", "rs11771145_A", "rs11136000_T",
                   "rs3851179_A", "rs17125944_C", "rs3764650_G")

df_model <- df_last[, selected_vars] %>% na.omit()

# Codifica binaria SNPs (0/1/2 → 0/1)
snp_vars <- selected_vars[-c(1, 2)]  # tutti tranne DIAGNOSIS e genotype
df_model[snp_vars] <- lapply(df_model[snp_vars], function(x) {
  x <- as.numeric(as.character(x))  # converte factor → numeric se serve
  ifelse(x == 0, 0, 1)              # 0 = wild-type, 1 = almeno una variante
})

# Raggruppa genotype in categorie APOE
df_model$genotype_group <- case_when(
  df_model$genotype %in% c("2/2", "3/2", "3/3") ~ "non-carrier",
  df_model$genotype %in% c("4/2", "4/3") ~ "heterozygous",
  df_model$genotype == "4/4" ~ "homozygous",
  TRUE ~ NA_character_
)
df_model$genotype_group <- factor(df_model$genotype_group, levels = c("non-carrier", "heterozygous", "homozygous"))

# Codifica binaria della diagnosi: AD vs NoAD (valori 1,2,3 → NoAD/AD)
df_model$DIAGNOSIS_BIN <- factor(ifelse(df_model$DIAGNOSIS == 3, "AD", "NoAD"),
                                 levels = c("NoAD", "AD"))

#Rimuovi variabili inutili
df_model$DIAGNOSIS <- NULL
df_model$genotype <- NULL

# REGRESSIONE LOGISTICA GENETICA
logit_genetico <- glm(DIAGNOSIS_BIN ~ ., 
                      data = df_model,
                      family = binomial)

# RISULTATI
cat("\n Coefficienti della regressione logistica:\n")
print(summary(logit_genetico))

# Odds ratio con intervallo di confidenza
cat("\n Odds Ratio e intervalli di confidenza:\n")
print(exp(cbind(OR = coef(logit_genetico), confint(logit_genetico))))

# Controllo multicollinearità
cat("\n Variance Inflation Factors (VIF):\n")
print(vif(logit_genetico))

# Frequenze delle variabili categoriali
cat("\n Frequenze delle categorie (dopo la trasformazione):\n")

cat("\n Genotype group:\n")
print(table(df_model$genotype_group))

for (snp in snp_vars) {
  cat(sprintf("\n %s (0=wild-type, 1=variante):\n", snp))
  print(table(df_model[[snp]]))
}


#----logit genetica rimuovendo le non significative----
# LIBRERIE
library(dplyr)
library(car)

# DATI E PREPROCESSING
df_last <- df %>%
  group_by(Subject_ID) %>%
  filter(Visit == max(Visit)) %>%
  ungroup()

selected_vars <- c("DIAGNOSIS", "genotype", 
                   "rs744373_C",  "rs11771145_A", 
                   "rs3851179_A")

df_model <- df_last[, selected_vars] %>% na.omit()

# Codifica binaria SNPs (0/1/2 → 0/1)
snp_vars <- selected_vars[-c(1, 2)]  # tutti tranne DIAGNOSIS e genotype
df_model[snp_vars] <- lapply(df_model[snp_vars], function(x) {
  x <- as.numeric(as.character(x))  # converte factor → numeric se serve
  ifelse(x == 0, 0, 1)              # 0 = wild-type, 1 = almeno una variante
})

# Raggruppa genotype in categorie APOE
df_model$genotype_group <- case_when(
  df_model$genotype %in% c("2/2", "3/2", "3/3") ~ "non-carrier",
  df_model$genotype %in% c("4/2", "4/3") ~ "heterozygous",
  df_model$genotype == "4/4" ~ "homozygous",
  TRUE ~ NA_character_
)
df_model$genotype_group <- factor(df_model$genotype_group, levels = c("non-carrier", "heterozygous", "homozygous"))

# Codifica binaria della diagnosi: AD vs NoAD (valori 1,2,3 → NoAD/AD)
df_model$DIAGNOSIS_BIN <- factor(ifelse(df_model$DIAGNOSIS == 3, "AD", "NoAD"),
                                 levels = c("NoAD", "AD"))

# Rimuovi variabili inutili
df_model$DIAGNOSIS <- NULL
df_model$genotype <- NULL

# REGRESSIONE LOGISTICA GENETICA
logit_genetico <- glm(DIAGNOSIS_BIN ~ ., 
                      data = df_model,
                      family = binomial)
# RISULTATI
cat("\n Coefficienti della regressione logistica:\n")
print(summary(logit_genetico))

# Odds ratio con intervallo di confidenza
cat("\n Odds Ratio e intervalli di confidenza:\n")
print(exp(cbind(OR = coef(logit_genetico), confint(logit_genetico))))

#  Controllo multicollinearità
cat("\n Variance Inflation Factors (VIF):\n")
print(vif(logit_genetico))

#  Frequenze delle variabili categoriali
cat("\n Frequenze delle categorie (dopo la trasformazione):\n")

cat("\n🔹 Genotype group:\n")
print(table(df_model$genotype_group))

for (snp in snp_vars) {
  cat(sprintf("\n🔹 %s (0=wild-type, 1=variante):\n", snp))
  print(table(df_model[[snp]]))
}

#libreria
library(ggplot2)

#  Estrai coefficienti, CI e OR
coef_table <- summary(logit_genetico)$coefficients
conf_int <- confint(logit_genetico)  # CI al 95%
OR <- exp(coef(logit_genetico))
OR_CI <- exp(conf_int)

# Crea data frame per il grafico
df_or <- data.frame(
  Variable = rownames(coef_table),
  OR = OR,
  CI_lower = OR_CI[, 1],
  CI_upper = OR_CI[, 2],
  p_value = coef_table[, 4]
)

# Rimuovi intercetta
df_or <- df_or[df_or$Variable != "(Intercept)", ]

# Forest plot con ggplot2
ggplot(df_or, aes(x = reorder(Variable, OR), y = OR)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip() +
  scale_y_log10() +
  labs(
    title = "Odds Ratio con Intervalli di Confidenza (95%)",
    x = "Variabile",
    y = "Odds Ratio (scala log)"
  ) +
  theme_minimal()


## come le variabili genetiche migliorano la classificazione con modello multinomiale
#-----Codice per testare una variabile genetica alla volta----
library(nnet)
library(dplyr)
library(caret)

# Variabili cerebrali significative dalla MANOVA
brain_vars <- c("Hippocampus_Total", "Precuneus_Total", "InfLatVentricle_Total")

# Prepara il dataset iniziale
df_model <- df_last %>%
  select(DIAGNOSIS, all_of(brain_vars),
         genotype, rs744373_C, rs11767557_C, rs11771145_A, 
         rs11136000_T, rs3851179_A, rs17125944_C, rs3764650_G)

# Codifica binaria SNPs (0/1/2 -> 0/1)
snp_vars <- c("rs744373_C", "rs11767557_C", "rs11771145_A", 
              "rs11136000_T", "rs3851179_A", "rs17125944_C", "rs3764650_G")
df_model[snp_vars] <- lapply(df_model[snp_vars], function(x) ifelse(x == 0, 0, 1))

# Raggruppa genotype in categorie APOE
df_model$genotype_group <- case_when(
  df_model$genotype %in% c("2/2", "3/2", "3/3") ~ "non-carrier",
  df_model$genotype %in% c("4/2", "4/3") ~ "heterozygous",
  df_model$genotype == "4/4" ~ "homozygous",
  TRUE ~ NA_character_
)
df_model$genotype_group <- factor(df_model$genotype_group, levels = c("non-carrier", "heterozygous", "homozygous"))

# Rimuovi la colonna originale del genotype
df_model$genotype <- NULL

# Modello base
base_data <- df_model %>% select(DIAGNOSIS, all_of(brain_vars)) %>% na.omit()
model_base <- multinom(DIAGNOSIS ~ ., data = base_data, trace = FALSE)
aic_base <- AIC(model_base)
acc_base <- mean(predict(model_base) == base_data$DIAGNOSIS)

cat("MODELLO BASE:\n")
cat("AIC:", aic_base, "\n")
cat("Accuracy:", round(acc_base, 3), "\n\n")

# Lista delle variabili genetiche da testare una alla volta
genetic_vars <- c(snp_vars, "genotype_group")

# Ciclo: una variabile genetica alla volta
for (var in genetic_vars) {
  formula_vars <- c("DIAGNOSIS", brain_vars, var)
  data_sub <- df_model[, formula_vars] %>% na.omit()
  
  model <- multinom(DIAGNOSIS ~ ., data = data_sub, trace = FALSE)
  preds <- predict(model)
  acc <- mean(preds == data_sub$DIAGNOSIS)
  aic <- AIC(model)
  
  cat("MODELLO CON:", var, "\n")
  cat("AIC:", aic, "\n")
  cat("Accuracy:", round(acc, 3), "\n")
  print(confusionMatrix(factor(preds), factor(data_sub$DIAGNOSIS)))
  cat("-----------------------------------------------------\n\n")
}


#-----modelli congiunti (genotype+1 snippet)-----

library(nnet)
library(dplyr)
library(caret)

# Variabili cerebrali selezionate
brain_vars <- c("Hippocampus_Total", "Precuneus_Total", "InfLatVentricle_Total")

# SNP da testare
snp_vars <- c("rs744373_C", "rs11767557_C", "rs11771145_A", 
              "rs11136000_T", "rs3851179_A", "rs17125944_C", "rs3764650_G")

# Prepara base comune
df_model_base <- df_last %>%
  select(DIAGNOSIS, all_of(brain_vars), genotype, all_of(snp_vars)) %>%
  na.omit()

# Codifica binaria SNPs
df_model_base[snp_vars] <- lapply(df_model_base[snp_vars], function(x) ifelse(x == 0, 0, 1))

# Raggruppa genotype
df_model_base$genotype_group <- case_when(
  df_model_base$genotype %in% c("2/2", "3/2", "3/3") ~ "non-carrier",
  df_model_base$genotype %in% c("4/3", "4/2") ~ "heterozygous",
  df_model_base$genotype == "4/4" ~ "homozygous",
  TRUE ~ NA_character_
)
df_model_base$genotype_group <- factor(df_model_base$genotype_group,
                                       levels = c("non-carrier", "heterozygous", "homozygous"))

# Rimuovi genotype originale
df_model_base$genotype <- NULL

# Rimuovi righe con NA
df_model_base <- na.omit(df_model_base)

# Testa combinazioni: 1 SNP + genotype_group + brain_vars
for (snp in snp_vars) {
  formula_comb <- as.formula(paste("DIAGNOSIS ~", 
                                   paste(c(brain_vars, snp, "genotype_group"), collapse = " + ")))
  
  model <- multinom(formula_comb, data = df_model_base, trace = FALSE)
  preds <- predict(model)
  acc <- mean(preds == df_model_base$DIAGNOSIS)
  aic <- AIC(model)
  
  cat("\n🧬 MODELLO CONGIUNTO CON:", snp, "+ genotype_group\n")
  cat("AIC:", aic, "\n")
  cat("Accuracy:", round(acc, 3), "\n")
  print(confusionMatrix(factor(preds), factor(df_model_base$DIAGNOSIS)))
  cat("--------------------------------------------------\n")
}

#-----combinazione di 3 snippets+genotype----
library(nnet)
library(dplyr)
library(caret)
library(gtools)  # per combinations()

# Variabili cerebrali
brain_vars <- c("Hippocampus_Total", "Precuneus_Total", "InfLatVentricle_Total")

# SNP disponibili
snp_vars <- c("rs744373_C", "rs11767557_C", "rs11771145_A", 
              "rs11136000_T", "rs3851179_A", "rs17125944_C", "rs3764650_G")

# Base dataset
df_model_base <- df_last %>%
  select(DIAGNOSIS, all_of(brain_vars), genotype, all_of(snp_vars)) %>%
  na.omit()

# Codifica SNP binaria
df_model_base[snp_vars] <- lapply(df_model_base[snp_vars], function(x) ifelse(x == 0, 0, 1))

# Genotype grouping
df_model_base$genotype_group <- case_when(
  df_model_base$genotype %in% c("2/2", "3/2", "3/3") ~ "non-carrier",
  df_model_base$genotype %in% c("4/3", "4/2") ~ "heterozygous",
  df_model_base$genotype == "4/4" ~ "homozygous",
  TRUE ~ NA_character_
)
df_model_base$genotype_group <- factor(df_model_base$genotype_group,
                                       levels = c("non-carrier", "heterozygous", "homozygous"))

df_model_base$genotype <- NULL
df_model_base <- na.omit(df_model_base)

# Scegli combinazioni da testare (es. 2 SNP alla volta)
combi_snp <- combinations(n = length(snp_vars), r = 3, v = snp_vars)

# Loop su combinazioni
for (i in 1:nrow(combi_snp)) {
  snps <- combi_snp[i, ]
  vars <- c(brain_vars, snps, "genotype_group")
  formula_comb <- as.formula(paste("DIAGNOSIS ~", paste(vars, collapse = " + ")))
  
  model <- multinom(formula_comb, data = df_model_base, trace = FALSE)
  preds <- predict(model)
  acc <- mean(preds == df_model_base$DIAGNOSIS)
  aic <- AIC(model)
  
  cat("\n MODELLO CONGIUNTO CON:", paste(snps, collapse = " + "), "+ genotype_group\n")
  cat("AIC:", round(aic, 3), "\n")
  cat("Accuracy:", round(acc, 3), "\n")
  cat("--------------------------------------------------\n")
}


#-----visualizzazione-----
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# Tabella delle combinazioni testate e relative performance
results <- tribble(
  ~combination, ~accuracy, ~aic,
  "rs3851179_A + genotype_group", 0.674, 232.9,
  "rs11771145_A + rs17125944_C + rs3851179_A", 0.667, 229.852,
  "rs11767557_C + rs11771145_A + rs3851179_A", 0.667, 234.05,
  "rs11767557_C + rs3851179_A + rs744373_C", 0.667, 237.86,
  "rs17125944_C + rs3764650_G + rs3851179_A", 0.667, 235.36,
  "rs11771145_A + rs3851179_A + rs744373_C", 0.667, 232.251,
  "rs11136000_T + rs11767557_C + rs3851179_A", 0.651, 239.117,
  "rs11771145_A + rs3764650_G + rs3851179_A", 0.651, 232.5,
  "rs11136000_T + rs17125944_C + rs3851179_A", 0.659, 235.942,
  "rs11767557_C + rs17125944_C + rs3851179_A", 0.659, 237.012,
  "rs17125944_C + rs3851179_A + rs744373_C", 0.659, 234.29,
  "rs3764650_G + rs3851179_A + rs744373_C", 0.659, 236.522
)

# Ordina per accuracy
results <- results %>%
  arrange(desc(accuracy)) %>%
  mutate(combination = factor(combination, levels = combination))

#  Barplot delle accuracy
ggplot(results, aes(x = combination, y = accuracy)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Accuracy dei modelli genetici", x = "Combinazione SNP + APOE", y = "Accuracy") +
  theme_minimal(base_size = 13)

# Heatmap: frequenza SNP tra i top performer (Accuracy ≥ 0.659)
snp_list <- c("rs744373_C", "rs11767557_C", "rs11771145_A", 
              "rs11136000_T", "rs3851179_A", "rs17125944_C", "rs3764650_G")

# Espandi combinazioni in SNP singoli
heatmap_data <- results %>%
  filter(accuracy >= 0.659) %>%
  rowwise() %>%
  mutate(snp_components = list(str_extract_all(combination, paste(snp_list, collapse = "|"))[[1]])) %>%
  unnest(snp_components) %>%
  count(snp_components, name = "freq")

# Heatmap
ggplot(heatmap_data, aes(x = "", y = snp_components, fill = freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "Frequenza SNP nei migliori modelli", x = NULL, y = "SNP") +
  theme_minimal(base_size = 13)


#----volcano plot---
# Required packages
library(nnet)
library(dplyr)
library(broom)
library(ggplot2)

# 1. Select and prepare your dataset
brain_vars <- c("Hippocampus_Total", "Precuneus_Total", "InfLatVentricle_Total")
snp_vars <- c("rs744373_C", "rs11767557_C", "rs11771145_A", 
              "rs11136000_T", "rs3851179_A", "rs17125944_C", "rs3764650_G")

# Create a binary outcome for volcano plot (e.g., AD vs CN)
df_binary <- df_last %>%
  filter(DIAGNOSIS %in% c("1", "3")) %>%  # CN = 1, AD = 3
  mutate(DIAGNOSIS = factor(DIAGNOSIS, levels = c("1", "3"))) %>%
  select(DIAGNOSIS, all_of(brain_vars), all_of(snp_vars)) %>%
  na.omit()

# Binarize SNPs
df_binary[snp_vars] <- lapply(df_binary[snp_vars], function(x) ifelse(x == 0, 0, 1))

# 2. Fit logistic regression model (AD vs CN)
model_bin <- glm(DIAGNOSIS ~ ., data = df_binary, family = binomial)

# 3. Extract coefficients and p-values
coefs <- tidy(model_bin)

# Filter only SNPs
coefs_snps <- coefs %>%
  filter(term %in% snp_vars) %>%
  mutate(logp = -log10(p.value))

# 4. Volcano plot
ggplot(coefs_snps, aes(x = estimate, y = logp, label = term)) +
  geom_point(size = 3) +
  geom_text(vjust = -1, size = 3.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  xlab("Log-Odds (Effect Size)") +
  ylab("-log10(p-value)") +
  ggtitle("Volcano Plot: SNP Effect Sizes vs Significance (AD vs CN)") +
  theme_minimal()

#CI 95% bootstrap

boot_acc_multinom_oob <- function(data, formula, B = 2000, seed = 123){
  set.seed(seed)
  n <- nrow(data)
  acc <- numeric(B)
  
  for(b in 1:B){
    idx <- sample.int(n, n, replace = TRUE)
    oob <- setdiff(seq_len(n), unique(idx))
    
    # se OOB troppo piccolo, salta (succede con N piccoli)
    if(length(oob) < 5){
      acc[b] <- NA
      next
    }
    
    d_train <- data[idx, , drop = FALSE]
    d_test  <- data[oob, , drop = FALSE]
    
    m <- multinom(formula, data = d_train, trace = FALSE)
    pred <- predict(m, newdata = d_test)
    acc[b] <- mean(pred == d_test$DIAGNOSIS)
  }
  
  acc <- acc[!is.na(acc)]
  c(
    acc_hat = mean(acc),
    ci_low  = quantile(acc, 0.025, names = FALSE),
    ci_high = quantile(acc, 0.975, names = FALSE),
    B_used  = length(acc)
  )
}

res_base_oob <- boot_acc_multinom_oob(base_data, DIAGNOSIS ~ ., B = 2000)
res_base_oob

vars <- c(brain_vars, "rs3851179_A", "genotype_group")
d_sub <- df_model_base %>% select(DIAGNOSIS, all_of(vars)) %>% na.omit()

res_joint_oob <- boot_acc_multinom_oob(d_sub, DIAGNOSIS ~ ., B = 2000)
res_joint_oob

#CI repeated cv

library(caret)
library(nnet)
library(dplyr)
# Forza DIAGNOSIS come fattore (classificazione)
df_last <- df_last %>%
  mutate(DIAGNOSIS = factor(DIAGNOSIS, levels = c(1, 2, 3),
                            labels = c("CN", "MCI", "AD")))
base_data <- df_last %>%
  select(DIAGNOSIS, all_of(brain_vars)) %>%
  na.omit()
res_base <- cv_acc_ci_multinom(base_data, DIAGNOSIS ~ ., k = 10, repeats = 50, seed = 123)

cv_acc_ci_multinom <- function(data, formula, k = 10, repeats = 50, seed = 123){
  
  data <- na.omit(data)
  
  ctrl <- trainControl(
    method = "repeatedcv",
    number = k,
    repeats = repeats,
    savePredictions = "final",
    classProbs = FALSE
  )
  
  set.seed(seed)
  fit <- train(
    formula,
    data = data,
    method = "multinom",
    trControl = ctrl,
    trace = FALSE
  )
  
  acc_vals <- fit$resample$Accuracy
  m <- length(acc_vals)
  acc_mean <- mean(acc_vals)
  acc_sd <- sd(acc_vals)
  
  ci <- acc_mean + qt(c(0.025, 0.975), df = m - 1) * acc_sd / sqrt(m)
  
  list(
    acc_mean = acc_mean,
    ci_low   = ci[1],
    ci_high  = ci[2],
    acc_sd   = acc_sd,
    m        = m,
    fit      = fit
  )
}

brain_vars <- c("Hippocampus_Total", "Precuneus_Total", "InfLatVentricle_Total")

base_data <- df_last %>%
  select(DIAGNOSIS, all_of(brain_vars)) %>%
  na.omit()

res_base <- cv_acc_ci_multinom(base_data, DIAGNOSIS ~ ., k = 10, repeats = 50, seed = 123)

cat("MODELLO BASE (brain vars):\n")
cat("Accuracy media:", round(res_base$acc_mean, 3),
    " | 95% CI [", round(res_base$ci_low, 3), ",", round(res_base$ci_high, 3), "]\n\n")

# Assicurati che l'outcome sia factor (classificazione)
df_model_base <- df_model_base %>%
  mutate(DIAGNOSIS = factor(DIAGNOSIS, levels = c(1, 2, 3), labels = c("CN", "MCI", "AD")))

results_joint_ci <- data.frame(
  model = character(),
  acc_mean = numeric(),
  ci_low = numeric(),
  ci_high = numeric(),
  stringsAsFactors = FALSE
)

for (snp in snp_vars) {
  
  d_sub <- df_model_base %>%
    select(DIAGNOSIS, all_of(brain_vars), all_of(snp), genotype_group) %>%
    na.omit()
  
  res <- cv_acc_ci_multinom(d_sub, DIAGNOSIS ~ ., k = 10, repeats = 50, seed = 123)
  
  results_joint_ci <- rbind(results_joint_ci, data.frame(
    model = paste("brain_vars +", snp, "+ genotype_group"),
    acc_mean = res$acc_mean,
    ci_low = res$ci_low,
    ci_high = res$ci_high
  ))
  
  cat("MODELLO CONGIUNTO:", snp, "+ genotype_group\n")
  cat("Accuracy media:", round(res$acc_mean, 3),
      " | 95% CI [", round(res$ci_low, 3), ",", round(res$ci_high, 3), "]\n")
  cat("--------------------------------------------------\n\n")
}

results_joint_ci <- results_joint_ci %>% arrange(desc(acc_mean))
print(results_joint_ci)

