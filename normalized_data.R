setwd("C:/Users/HP/OneDrive/Documenti/APPLIED STATISTICS/PROGETTO")

library(dplyr)
library(caret)
library(readxl)
library(purrr)

# percorso file
file_path <- "Alzhemeir_Dataset.xlsx"  

# prendo i nomi di tutti i fogli
sheet_names <- excel_sheets(file_path)

# leggo e concateno, aggiungendo visit = 1,2,3,...
df_final <- map2_dfr(
  .x = sheet_names,
  .y = seq_along(sheet_names),
  ~ read_excel(path = file_path, sheet = .x) %>%
    mutate(Visit = .y)
)

library(writexl)
write_xlsx(df_final, "output_concatenato.xlsx")


df <-df_final

df <- df %>%
  mutate(
    # Strutture limbiche e talamo
    Hippocampus_Total        = `Left-Hippocampus`        + `Right-Hippocampus`,
    Thalamus_Total           = `Left-Thalamus`           + `Right-Thalamus`,
    
    # Cortecce temporali / insula / precuneus
    SuperiorTemporal_Total   = `ctx-lh-superiortemporal` + `ctx-rh-superiortemporal`,
    Insula_Total             = `ctx-lh-insula`           + `ctx-rh-insula`,
    Precuneus_Total          = `ctx-lh-precuneus`        + `ctx-rh-precuneus`,
    
    # Corpi mammillari
    Mammilare_Total          = `L-C.mammilare`           + `R-C.mammilare`,
    
    # Caudal middle frontal
    CaudalMiddleFrontal_Total = `ctx-lh-caudalmiddlefrontal` + `ctx-rh-caudalmiddlefrontal`,
    
    # Ventricoli inferiori laterali
    InfLatVentricle_Total    = `Left-Inf-Lat-Vent`       + `Right-Inf-Lat-Vent`
  )

df <- df %>%
  mutate(
    Third_Vent_combined = coalesce(`Third-Ventricle`, `3rd-Ventricle`)
  )

df <- df %>%
  mutate(
    ICV_proxy = `Left-Cerebral-White-Matter` +
      `Right-Cerebral-White-Matter` +
      CSF +
      `Left-Lateral-Ventricle` +
      `Right-Lateral-Ventricle` +
      `Left-Inf-Lat-Vent` +
      `Right-Inf-Lat-Vent` +
      Third_Vent_combined
  )

df <- df %>%
  mutate(
    Hippocampus_ICV          = Hippocampus_Total        / ICV_proxy,
    Thalamus_ICV             = Thalamus_Total           / ICV_proxy,
    SuperiorTemporal_ICV     = SuperiorTemporal_Total   / ICV_proxy,
    Insula_ICV               = Insula_Total             / ICV_proxy,
    Precuneus_ICV            = Precuneus_Total          / ICV_proxy,
    Mammilare_ICV            = Mammilare_Total          / ICV_proxy,
    CaudalMiddleFrontal_ICV  = CaudalMiddleFrontal_Total / ICV_proxy,
    InfLatVentricle_ICV      = InfLatVentricle_Total    / ICV_proxy
  )

df <- df %>%
  mutate(
    Hippocampus_ICVz         = as.numeric(scale(Hippocampus_ICV)),
    Thalamus_ICVz            = as.numeric(scale(Thalamus_ICV)),
    SuperiorTemporal_ICVz    = as.numeric(scale(SuperiorTemporal_ICV)),
    Insula_ICVz              = as.numeric(scale(Insula_ICV)),
    Precuneus_ICVz           = as.numeric(scale(Precuneus_ICV)),
    Mammilare_ICVz           = as.numeric(scale(Mammilare_ICV)),
    CaudalMiddleFrontal_ICVz = as.numeric(scale(CaudalMiddleFrontal_ICV)),
    InfLatVentricle_ICVz     = as.numeric(scale(InfLatVentricle_ICV))
  )
final_df <- df %>%
  select(
    Subject_ID,
    DIAGNOSIS,
    gender,
    genotype,
    rs744373_C,
    rs11767557_C,
    rs11771145_A,
    rs11136000_T,
    rs3851179_A,
    rs17125944_C,
    rs3764650_G,
    Visit,
    Hippocampus_Total        = Hippocampus_ICVz,
    Thalamus_Total           = Thalamus_ICVz,
    SuperiorTemporal_Total   = SuperiorTemporal_ICVz,
    Insula_Total             = Insula_ICVz,
    Precuneus_Total          = Precuneus_ICVz,
    Mammilare_Total          = Mammilare_ICVz,
    CaudalMiddleFrontal_Total = CaudalMiddleFrontal_ICVz,
    InfLatVentricle_Total    = InfLatVentricle_ICVz
  )

#Salva su file
write_xlsx(final_df, "Dataset_Alz_ICV_final.xlsx")

