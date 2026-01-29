# ================================
# LIBRERIE
# ================================
library(mvtnorm)
library(rgl)
library(car)
library(MVN)
library(biotools)
library(ggplot2)

# ================================
# LETTURA DATI
# ================================
alz <- read.csv("alz_finale_unico.csv")

alz$DIAGNOSIS <- factor(alz$DIAGNOSIS)
alz$gender <- factor(alz$gender)

# ================================
# VARIABILI MRI
# ================================
brain_vars <- c(
  "Hippocampus_Total",
  "Thalamus_Total",
  "SuperiorTemporal_Total",
  "Insula_Total",
  "Precuneus_Total",
  "Mammilare_Total",
  "CaudalMiddleFrontal_Total",
  "InfLatVentricle_Total"
)

# ================================
# DATASET COERENTE (MRI + DIAGNOSIS)
# ================================
alz_sub <- alz[, c(brain_vars, "DIAGNOSIS")]

# rimozione NA coerente
alz_sub <- na.omit(alz_sub)

# ================================
# STANDARDIZZAZIONE (IMPORTANTE)
# ================================
alz_sub_scaled <- scale(alz_sub[, brain_vars])

# ================================
# CLUSTERING GERARCHICO
# ================================
dist.matrix <- dist(alz_sub_scaled, method = "euclidean")
cluster.dl <- hclust(dist.matrix, method = "ward.D2")

alz_sub$cluster <- factor(cutree(cluster.dl, k = 2))

# ================================
# DENDROGRAMMA
# ================================
plot(cluster.dl, main = "Dendrogram", hang = -0.1, labels = FALSE)
rect.hclust(cluster.dl, k = 2)

# ================================
# COPHENETIC CORRELATION
# ================================
coph.coeff <- cor(dist.matrix, cophenetic(cluster.dl))
cat("Cophenetic correlation coefficient:", coph.coeff, "\n")

# ================================
# ASSOCIAZIONE CLUSTER vs DIAGNOSIS
# ================================
tab_cluster_diag <- table(alz_sub$cluster, alz_sub$DIAGNOSIS)
print(tab_cluster_diag)

# scegli automaticamente il test corretto
if (any(chisq.test(tab_cluster_diag)$expected < 5)) {
  assoc_test <- fisher.test(tab_cluster_diag)
  cat("Fisher exact test\n")
} else {
  assoc_test <- chisq.test(tab_cluster_diag)
  cat("Chi-square test\n")
}

print(assoc_test)

# ================================
# MANOVA
# ================================
fit.manova <- manova(as.matrix(alz_sub_scaled) ~ cluster, data = alz_sub)

summary(fit.manova, test = "Wilks")
summary(fit.manova, test = "Pillai")

# ================================
# TEST ASSUNZIONI
# ================================
for (k in levels(alz_sub$cluster)) {
  cat("Cluster", k, "\n")
  print(
    mvn(
      alz_sub_scaled[alz_sub$cluster == k, ],
      mvnTest = "hz"
    )
  )
}

boxM(alz_sub_scaled, alz_sub$cluster)

# ================================
# BOXPLOT
# ================================
alz_sub$DIAGNOSIS_GROUPED <- ifelse(
  alz_sub$DIAGNOSIS %in% c(1, 2),
  "Non-Alzheimer",
  "Alzheimer"
)
alz_sub$DIAGNOSIS_GROUPED <- factor(alz_sub$DIAGNOSIS_GROUPED)

for (var in brain_vars) {
  
  p1 <- ggplot(alz_sub, aes(x = cluster, y = .data[[var]], fill = cluster)) +
    geom_boxplot(alpha = 0.7) +
    labs(title = paste(var, "by Cluster")) +
    theme_minimal()
  
  p2 <- ggplot(alz_sub, aes(x = DIAGNOSIS_GROUPED, y = .data[[var]], fill = DIAGNOSIS_GROUPED)) +
    geom_boxplot(alpha = 0.7) +
    labs(title = paste(var, "by Diagnosis")) +
    theme_minimal()
  
  print(p1)
  print(p2)
}

