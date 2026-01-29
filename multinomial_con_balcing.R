

library(nnet)
library(cluster)
library(ggplot2)
library(tidyr)
library(dplyr)
library(smotefamily)
library(tibble)

#In che modo le diverse misure volumetriche cerebrali influenzano la probabilità di 
#appartenenza a ciascuna classe diagnostica (CN, MCI, AD)?

#----fase1-----

# Carica il dataset
alz <- read_excel("Dataset_Alz_ICV_final.xlsx")

# Trasforma la variabile DIAGNOSIS in fattore con etichette esplicite
alz$DIAGNOSIS <- factor(alz$DIAGNOSIS,
                        levels = c(1, 2, 3),
                        labels = c("CN", "MCI", "AD"))


# Conta i valori mancanti per variabile
colSums(is.na(alz))

# Rimuove le righe con NA (opzione consigliata per 2 casi mancanti)
alz <- na.omit(alz)

# Trasforma anche gender in fattore (se non lo è già)
alz$gender <- as.factor(alz$gender)

# Statistiche descrittive per le variabili numeriche
summary(select(alz, where(is.numeric)))

ggplot(alz, aes(x = DIAGNOSIS, fill = DIAGNOSIS)) +
  geom_bar() +
  scale_fill_manual(values = c("CN" = "#90ee90",     # verde chiaro
                               "MCI" = "#006400",    # verde scuro
                               "AD" = "#800080")) +  # viola
  labs(title = "Distribuzione delle diagnosi",
       x = "Diagnosi", y = "Frequenza") +
  theme_minimal()


#grafici esplorativi (boxplot)

# Elenca le variabili numeriche cerebrali
brain_vars <- c("Hippocampus_Total", "Thalamus_Total", "SuperiorTemporal_Total", 
                "Insula_Total", "Precuneus_Total", "Mammilare_Total", 
                "CaudalMiddleFrontal_Total", "InfLatVentricle_Total")

# Filtra eventuali NA
alz_clean <- na.omit(alz)

# Crea un long dataframe per ggplot
alz_long <- alz_clean %>%
  select(DIAGNOSIS, all_of(brain_vars)) %>%
  pivot_longer(cols = -DIAGNOSIS, names_to = "Region", values_to = "Value")

# Boxplot per ciascuna regione cerebrale
ggplot(alz_long, aes(x = DIAGNOSIS, y = Value, fill = DIAGNOSIS)) +
  geom_boxplot() +
  facet_wrap(~Region, scales = "free", ncol = 2) +
  labs(title = "Distribuzione delle misure cerebrali per diagnosi",
       x = "Diagnosi", y = "Valore (standardizzato)") +
  theme_minimal() +
  theme(legend.position = "none")

# Fase 2: Modello multinomiale con class weights

# Calcola la frequenza per classe
class_freq <- table(alz_clean$DIAGNOSIS)

# Crea un vettore di pesi della stessa lunghezza del dataset
class_weights <- alz_clean$DIAGNOSIS %>%
  as.character() %>%           # trasforma in carattere
  sapply(function(x) 1 / class_freq[x]) %>%   # assegna peso 1/frequenza
  as.numeric()                # assicurati che siano numeri


# Fit del modello multinomiale pesato
model_weighted <- multinom(
  DIAGNOSIS ~ Hippocampus_Total + Thalamus_Total + SuperiorTemporal_Total +
    Insula_Total + Precuneus_Total + Mammilare_Total +
    CaudalMiddleFrontal_Total + InfLatVentricle_Total,
  data = alz,
  weights = class_weights
)

# Sommario del modello
summary(model_weighted)

# Calcolo dei p-value
z <- summary(model_weighted)$coefficients / summary(model_weighted)$standard.errors
p <- 2 * (1 - pnorm(abs(z)))
print(p)

#nessuna variabile è significativa :(, cambiamo metodo per balancing!

#fase2: SMOTE + Multinom

# Rimuovi variabili non numeriche (Subject_ID, gender)
alz_smote <- alz_clean %>% 
  select(-Subject_ID, -gender)

# Codifica numerica della diagnosi
alz_smote$DIAGNOSIS_NUM <- as.numeric(alz_clean$DIAGNOSIS)

# Seleziona solo variabili numeriche per SMOTE
names(alz_smote)
X <- alz_smote %>% select(-genotype,-rs744373_C,-rs11767557_C,-rs11771145_A, -rs11136000_T,-rs3851179_A,-rs17125944_C, -rs3764650_G,-Visit, -DIAGNOSIS, -DIAGNOSIS_NUM)
target <- alz_smote$DIAGNOSIS_NUM

# Applica SMOTE multiclass
set.seed(123)
smote_result <- SMOTE(X = X, target = target, K = 5)

# Ricostruiamo il nuovo dataframe bilanciato
alz_balanced <- smote_result$data
alz_balanced$DIAGNOSIS <- factor(alz_balanced$class,
                                 levels = c(1, 2, 3),
                                 labels = c("CN", "MCI", "AD"))

# Fit modello multinomiale sul dataset bilanciato
model_smote <- multinom(
  DIAGNOSIS ~ Hippocampus_Total + Thalamus_Total + SuperiorTemporal_Total +
    Insula_Total + Precuneus_Total + Mammilare_Total +
    CaudalMiddleFrontal_Total + InfLatVentricle_Total,
  data = alz_balanced
)

# Riassunto
summary(model_smote)

# p-value
z <- summary(model_smote)$coefficients / summary(model_smote)$standard.errors
p <- 2 * (1 - pnorm(abs(z)))
print(p)
#alcune variabili sono significative! soprattutto per distinguere AD da CN (per MCI no) ma a noi interessa questo!

# Codice completo per plottare tutte le variabili cerebrali
# Funzione per plottare l'effetto di una covariata
plot_multinom_effect <- function(model, data, covariate) {
  x_seq <- seq(min(data[[covariate]], na.rm = TRUE),
               max(data[[covariate]], na.rm = TRUE),
               length.out = 100)
  
  newdata <- data[rep(1, length(x_seq)), ]
  newdata[] <- lapply(newdata, function(x) if(is.numeric(x)) mean(x, na.rm = TRUE) else x[1])
  newdata[[covariate]] <- x_seq
  
  probs <- predict(model, newdata = newdata, type = "probs")
  probs_df <- as.data.frame(probs)
  probs_df[[covariate]] <- x_seq
  
  probs_long <- pivot_longer(probs_df, cols = -all_of(covariate),
                             names_to = "DIAGNOSIS", values_to = "Probability")
  
  ggplot(probs_long, aes_string(x = covariate, y = "Probability", color = "DIAGNOSIS")) +
    geom_line(size = 1.2) +
    labs(x = covariate, y = "Probability") +
    theme_minimal()
}

# Lista delle variabili cerebrali nel modello
covariate_list <- c("Hippocampus_Total", "Thalamus_Total", "SuperiorTemporal_Total", 
                    "Insula_Total", "Precuneus_Total", "Mammilare_Total", 
                    "CaudalMiddleFrontal_Total", "InfLatVentricle_Total")

# Loop su ogni variabile e stampa il grafico
for (cov in covariate_list) {
  print(plot_multinom_effect(model_smote, alz, cov))
}


#----- Modello 2 smote: Ridotto (5 variabili significative)----
#Hippocampus_Total, InfLatVentricle_Total, Precuneus_Total,Insula_Total, Thalamus_Total
# Modello 1: 5 variabili
model_5vars <- multinom(
  DIAGNOSIS ~ Hippocampus_Total + InfLatVentricle_Total + 
    Precuneus_Total + Insula_Total + Thalamus_Total,
  data = alz_balanced
)

# Sommario + p-value
summary(model_5vars)
z5 <- summary(model_5vars)$coefficients / summary(model_5vars)$standard.errors
p5 <- 2 * (1 - pnorm(abs(z5)))
print(p5)
#tutte significative per AD

#Modello 3 smote: Parsimonioso (3 variabili forti)
# Modello 2: 3 variabili chiave
model_3vars <- multinom(
  DIAGNOSIS ~ Hippocampus_Total + InfLatVentricle_Total + Precuneus_Total,
  data = alz_balanced
)

# Sommario + p-value
summary(model_3vars)
z3 <- summary(model_3vars)$coefficients / summary(model_3vars)$standard.errors
p3 <- 2 * (1 - pnorm(abs(z3)))
print(p3)

# Effect plot per le 5 variabili del modello finale

# Funzione per plottare gli effetti marginali
plot_multinom_effect <- function(model, data, covariate) {
  x_seq <- seq(min(data[[covariate]], na.rm = TRUE),
               max(data[[covariate]], na.rm = TRUE),
               length.out = 100)
  
  newdata <- data[rep(1, length(x_seq)), ]
  newdata[] <- lapply(newdata, function(x) if (is.numeric(x)) mean(x, na.rm = TRUE) else x[1])
  newdata[[covariate]] <- x_seq
  
  probs <- predict(model, newdata = newdata, type = "probs")
  probs_df <- as.data.frame(probs)
  probs_df[[covariate]] <- x_seq
  
  probs_long <- pivot_longer(probs_df, cols = -all_of(covariate),
                             names_to = "DIAGNOSIS", values_to = "Probability")
  
  ggplot(probs_long, aes_string(x = covariate, y = "Probability", color = "DIAGNOSIS")) +
    geom_line(size = 1.2) +
    labs(title = paste("Effetto di", covariate),
         x = covariate, y = "Probabilità Predetta") +
    theme_minimal()
}

# Lista delle variabili nel modello ridotto
vars_5 <- c("Hippocampus_Total", "InfLatVentricle_Total", "Precuneus_Total",
            "Insula_Total", "Thalamus_Total")

# Ciclo per generare i grafici
for (v in vars_5) {
  print(plot_multinom_effect(model_5vars, alz_balanced, v))
}


# smote distribuzione nuova

ggplot(alz_balanced, aes(x = DIAGNOSIS, fill = DIAGNOSIS)) +
  geom_bar() +
  labs(title = "Distribuzione delle diagnosi dopo SMOTE",
       x = "Diagnosi", y = "Frequenza") +
  theme_minimal()



#-----undersampling----
# Sottocampioniamo
alz_under <- alz_clean %>%
  group_by(DIAGNOSIS) %>%
  sample_n(size = 37) %>%
  ungroup()

# Verifica
table(alz_under$DIAGNOSIS)

# Modello completo su dataset sottocampionato (undersampling)
model_under_full <- multinom(
  DIAGNOSIS ~ Hippocampus_Total + Thalamus_Total + SuperiorTemporal_Total +
    Insula_Total + Precuneus_Total + Mammilare_Total +
    CaudalMiddleFrontal_Total + InfLatVentricle_Total,
  data = alz_under
)

# Sommario + z e p-value
summary(model_under_full)

z_under <- summary(model_under_full)$coefficients / summary(model_under_full)$standard.errors
p_under <- 2 * (1 - pnorm(abs(z_under)))
print(p_under)

#Modello ridotto (4 variabili significative per AD)

model_under_5vars <- multinom(
  DIAGNOSIS ~ Hippocampus_Total + InfLatVentricle_Total +
    Precuneus_Total + Thalamus_Total,
  data = alz_under
)

# Sommario + p-value
summary(model_under_5vars)
z5_under <- summary(model_under_5vars)$coefficients / summary(model_under_5vars)$standard.errors
p5_under <- 2 * (1 - pnorm(abs(z5_under)))
print(p5_under)


# Modello minimo a 3 variabili
model_under_3vars <- multinom(
  DIAGNOSIS ~ Hippocampus_Total + InfLatVentricle_Total + Precuneus_Total,
  data = alz_under
)

# Sommario + p-value
summary(model_under_3vars)
z3_under <- summary(model_under_3vars)$coefficients / summary(model_under_3vars)$standard.errors
p3_under <- 2 * (1 - pnorm(abs(z3_under)))
print(p3_under)

# Funzione per plottare l'effetto marginale di una covariata
plot_multinom_effect <- function(model, data, covariate) {
  x_seq <- seq(min(data[[covariate]], na.rm = TRUE),
               max(data[[covariate]], na.rm = TRUE),
               length.out = 100)
  
  newdata <- data[rep(1, length(x_seq)), ]
  newdata[] <- lapply(newdata, function(x) if (is.numeric(x)) mean(x, na.rm = TRUE) else x[1])
  newdata[[covariate]] <- x_seq
  
  probs <- predict(model, newdata = newdata, type = "probs")
  probs_df <- as.data.frame(probs)
  probs_df[[covariate]] <- x_seq
  
  probs_long <- pivot_longer(probs_df, cols = -all_of(covariate),
                             names_to = "DIAGNOSIS", values_to = "Probability")
  
  ggplot(probs_long, aes_string(x = covariate, y = "Probability", color = "DIAGNOSIS")) +
    geom_line(size = 1.2) +
    labs(title = paste("Effetto di", covariate),
         x = covariate, y = "Probabilità predetta") +
    theme_minimal()
}

# Variabili finali selezionate
final_vars_under <- c("Hippocampus_Total", "InfLatVentricle_Total",
                      "Precuneus_Total",  "Thalamus_Total")

# Ciclo per generare e mostrare i grafici
for (v in final_vars_under) {
  print(plot_multinom_effect(model_under_5vars, alz_under, v))
}


# K-Medoids undersampling + Multinomial Model
# Librerie
library(dplyr)
library(cluster)
library(tibble)
library(nnet)

# Funzione corretta per selezionare i medoids in ciascun gruppo
get_medoids_sample <- function(df, group_key, k) {
  df_numeric <- df %>% select(where(is.numeric))
  
  # Se ci sono meno di 2 righe, ritorna il gruppo originale
  if (nrow(df_numeric) < 2) return(df)
  
  # Imposta k al massimo come n - 1 per evitare errori in pam()
  k_effettivo <- min(k, nrow(df_numeric) - 1)
  
  pam_fit <- pam(df_numeric, k = k_effettivo)
  df[pam_fit$id.med, ]
}

# Caricamento e preparazione dati
alz_kmed <- read_excel("Dataset_Alz_ICV_final.xlsx")
alz_kmed$DIAGNOSIS <- factor(alz_kmed$DIAGNOSIS, levels = c(1, 2, 3), labels = c("CN", "MCI", "AD"))
alz_kmed <- na.omit(alz_kmed)

# Selezione delle variabili
alz_kmed <- alz_kmed %>%
  select(DIAGNOSIS, Hippocampus_Total, InfLatVentricle_Total,
         Precuneus_Total, Thalamus_Total, Insula_Total, SuperiorTemporal_Total, Mammilare_Total, CaudalMiddleFrontal_Total)

# Trova numero minimo soggetti per gruppo
min_n <- alz_kmed %>% count(DIAGNOSIS) %>% pull(n) %>% min()

# Applica undersampling con k-medoids
alz_kmed_bal <- alz_kmed %>%
  group_by(DIAGNOSIS) %>%
  group_modify(~ get_medoids_sample(.x, .y, min_n)) %>%
  ungroup()

# Verifica bilanciamento
print(table(alz_kmed_bal$DIAGNOSIS))

#modello completo
model_kmed <- multinom(
  DIAGNOSIS ~ Hippocampus_Total + InfLatVentricle_Total + Precuneus_Total +
    Thalamus_Total + Insula_Total+ SuperiorTemporal_Total + Mammilare_Total + CaudalMiddleFrontal_Total,
  data = alz_kmed_bal
)

# Calcolo dei p-value
summary(model_kmed)
z_kmed <- summary(model_kmed)$coefficients / summary(model_kmed)$standard.errors
p_kmed <- 2 * (1 - pnorm(abs(z_kmed)))
print(p_kmed)

# Modello multinomiale con 4 variabili
model_kmed_5 <- multinom(
  DIAGNOSIS ~ Hippocampus_Total + InfLatVentricle_Total + Precuneus_Total +
    Thalamus_Total ,
  data = alz_kmed_bal
)

# Calcolo dei p-value
summary(model_kmed_5)
z_kmed_5 <- summary(model_kmed_5)$coefficients / summary(model_kmed_5)$standard.errors
p_kmed_5 <- 2 * (1 - pnorm(abs(z_kmed_5)))
print(p_kmed_5)

# Funzione per plottare effetto marginale
plot_multinom_effect <- function(model, data, covariate) {
  x_seq <- seq(min(data[[covariate]], na.rm = TRUE),
               max(data[[covariate]], na.rm = TRUE),
               length.out = 100)
  
  newdata <- data[rep(1, length(x_seq)), ]
  newdata[] <- lapply(newdata, function(x) if (is.numeric(x)) mean(x, na.rm = TRUE) else x[1])
  newdata[[covariate]] <- x_seq
  
  probs <- predict(model, newdata = newdata, type = "probs")
  probs_df <- as.data.frame(probs)
  probs_df[[covariate]] <- x_seq
  
  probs_long <- pivot_longer(probs_df, cols = -all_of(covariate),
                             names_to = "DIAGNOSIS", values_to = "Probability")


ggplot(probs_long, aes_string(x = covariate, y = "Probability", color = "DIAGNOSIS")) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c(
    "CN" = "#90ee90",     # light green
    "MCI" = "#006400",    # dark green
    "AD" = "#800080"      # purple
  )) +
  labs(title = paste("Effetto di", covariate),
       x = covariate, y = "Probabilità predetta") +
  theme_minimal()
}

# Variabili finali del modello ridotto
final_vars_kmed <- c("Hippocampus_Total", "InfLatVentricle_Total", "Precuneus_Total", "Thalamus_Total")

# Plot per ciascuna variabile
for (v in final_vars_kmed) {
  print(plot_multinom_effect(model_kmed_5, alz_kmed_bal, v))
}

#----#
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(dplyr)

# PCA sui dati originali
alz_pca <- alz_kmed %>% select(-DIAGNOSIS)
pca_result <- PCA(alz_pca, scale.unit = TRUE, graph = FALSE)

# PCA plot: aggiungi etichette per medoids
alz_kmed$subset <- "Original"
alz_kmed_bal$subset <- "Medoid"

alz_pca_plot <- bind_rows(alz_kmed, alz_kmed_bal)

fviz_pca_ind(pca_result,
             habillage = alz_pca_plot$DIAGNOSIS,
             addEllipses = TRUE,
             col.ind = ifelse(alz_pca_plot$subset == "Medoid", "red", "grey"),
             label = "none",
             pointsize = 2.5) +
  ggtitle("PCA - Medoid Samples vs Full Dataset") +
  theme_minimal()


#### Mann–Whitney U test

library(dplyr)
install.packages("effsize")
library(effsize)

# 1. Sottoinsieme: solo CN e AD sul dataset reale (no SMOTE)
real_data <- alz_clean %>%
  filter(DIAGNOSIS %in% c("CN", "AD")) %>%
  mutate(
    # assicuro che CN sia il riferimento (così d > 0 significa AD > CN)
    DIAGNOSIS = factor(DIAGNOSIS),
    DIAGNOSIS = relevel(DIAGNOSIS, ref = "CN")
  )

# 2. Statistiche descrittive
thalamus_summary <- real_data %>%
  group_by(DIAGNOSIS) %>%
  summarise(
    n        = n(),
    mean_T   = mean(Thalamus_Total, na.rm = TRUE),
    sd_T     = sd(Thalamus_Total,   na.rm = TRUE)
  )

thalamus_summary

# 3. Mann–Whitney U test (non parametrico)
mw_thal <- wilcox.test(
  Thalamus_Total ~ DIAGNOSIS,
  data   = real_data,
  exact  = FALSE,      # per evitare warning su campioni grandi
  conf.int = TRUE
)

mw_thal

library(dplyr)

#calcolare il coefficiente di correlazione di Pearson tra il volume del talamo e il volume del ventricolo laterale inferiore, 
#separatamente per i gruppi AD e CN

# Sottoinsieme reale: solo CN e AD
ventricular_data <- alz %>%
  filter(DIAGNOSIS %in% c("CN", "AD")) %>%
  mutate(DIAGNOSIS = factor(DIAGNOSIS))

# Correlazioni Pearson e p-value per ciascun gruppo
ventricular_corr <- ventricular_data %>%
  group_by(DIAGNOSIS) %>%
  summarise(
    n        = n(),
    r        = cor(Thalamus_Total, InfLatVentricle_Total,
                   use = "complete.obs", method = "pearson"),
    p_value  = cor.test(Thalamus_Total, InfLatVentricle_Total)$p.value,
    .groups  = "drop"
  )

ventricular_corr