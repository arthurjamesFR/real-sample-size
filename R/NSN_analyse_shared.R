# ==============================================================================
# HETEROGENEITY OF BASELINE RISK AND SAMPLE SIZE OF MORTALITY TRIALS
# MULTICENTER VERSION — Traumabase 

# ==============================================================================
# STEP 0 — PACKAGES, PARAMETERS, LABELS, GRAPHIC THEME
# ==============================================================================

library(readxl); library(tidyverse); library(mice); library(patchwork)

MODE <- "FINAL"

# ---- Paths -------------------------------------------------------------------
DATA_PATH  <- "path/to/traumabase_extract.xlsx"
DATA_SHEET <- "Golden"
CSV_DELIM  <- ";"
OUTPUT_DIR <- "outputs"
LOCAL_DIR  <- "outputs_local"

# ---- Simulation parameters ---------------------------------------------------
if (MODE == "FINAL") {
  M         <- 30
  MAXIT     <- 10
  NSIM_ADJ  <- 2000
  NSIM_VAL  <- 2000
  NSIM_DIAG <- 2000
  B_CENTRE  <- 200
  RUN_ADJ_CENTRE <- TRUE
} else {
  M         <- 5
  MAXIT     <- 3
  NSIM_ADJ  <- 200
  NSIM_VAL  <- 500
  NSIM_DIAG <- 200
  B_CENTRE  <- 20
  RUN_ADJ_CENTRE <- FALSE
}

POWER  <- 0.80
ALPHA  <- 0.05
OR_REF <- 0.70
SEED   <- 42

AUC_CIBLE  <- 0.90
AUC_SENS   <- c(0, 0.92, 0.90, 0.85)
Q_QUAD     <- 5000
K_QUAD     <- 201

# ---- Multicenter and confidentiality parameters ------------------------------
AGE_MIN        <- 18
SEUIL_CELLULE  <- 10
MIN_N_CENTRE_INCL <- 100
MIN_N_CENTRE   <- 200
MIN_DC_CENTRE  <- 20
SEUIL_NR       <- 0.50
MAX_COMP_MANQ  <- 2
N_CORES        <- max(1, parallel::detectCores() - 1)

# ---- Category labels (CHECK against the output of STEP 1) --------------------
OUI <- "Oui"; NON <- "Non"
LIB_MASCULIN    <- "Masculin"
MECA_PENETRANTS <- c("Arme à feu", "Arme blanche", "Autre traumatisme pénétrant")
CODES_INCONNU   <- c("Inconnu", "inconnu", "Non renseigné", "Non réalisé",
                     "6 (Inconnu)", "6, (Inconnu)")
HARMONISER      <- c("Mydriase Bilatérale" = "Mydriase bilatérale")

# ---- Meaning of "NA" (not applicable) by variable ----------------------------
NA_VAUT <- list(
  CGR_6H             = "0",
  PIC                = NON,
  REA_CATECHO        = NON,
  ACIDE_TRANEXAMIQUE = NON,
  IOT_PREHOSP        = NON,
  BLOC_RI_24H        = NON)

set.seed(SEED)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOCAL_DIR,  showWarnings = FALSE, recursive = TRUE)
sink(file.path(OUTPUT_DIR, "run_log.txt"), split = TRUE)

theme_journal <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_blank(), plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      axis.line   = element_line(color = "black", linewidth = 0.3),
      axis.ticks  = element_line(color = "black", linewidth = 0.3),
      axis.text   = element_text(color = "black", size = rel(0.9)),
      axis.title  = element_text(color = "black", size = rel(1.0)),
      panel.grid  = element_blank(),
      legend.position = "bottom", legend.title = element_blank(),
      legend.text = element_text(size = rel(0.85)), legend.key = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = rel(0.95), hjust = 0))
}

logit <- function(p) log(p / (1 - p))
expit <- function(x) 1 / (1 + exp(-x))
clamp <- function(p, eps = 1e-4) pmin(pmax(p, eps), 1 - eps)
num   <- function(x) suppressWarnings(as.numeric(gsub(",", ".", x)))
oui   <- function(x) ifelse(is.na(x), NA, x == OUI)
echant <- function(x, n = 2e5) if (length(x) > n) sample(x, n) else x

masquer_effectifs <- function(tb) {
  out <- as.character(tb)
  out[tb > 0 & tb < SEUIL_CELLULE] <- paste0("<", SEUIL_CELLULE)
  setNames(out, names(tb))
}

# ==============================================================================
# STEP 1 — LOADING, VARIABLE SELECTION, CATEGORY CHECKS
# ==============================================================================

lire_extraction <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls"))
    return(read_excel(path, sheet = DATA_SHEET, col_types = "text"))
  if (ext %in% c("csv", "txt"))
    return(read_delim(path, delim = CSV_DELIM,
                      col_types = cols(.default = col_character()),
                      na = "", locale = locale(encoding = "UTF-8"),
                      progress = FALSE))
  stop("Unsupported format: ", ext)
}

raw <- lire_extraction(DATA_PATH)
cat(sprintf("\n[1] Extract loaded: %d rows, %d columns\n",
            nrow(raw), ncol(raw)))

VARS_ID <- c("SITE_ID", "NUM_PATIENT")
VARS_ANALYSE <- c(
  "ANNEE", "BASE_SOURCE", "CRITERE_INCLUSION",
  "AGE", "AGE_PED", "SEXE", "MECANISME_CAUSE",
  "GLASGOW_INITIAL", "GLASGOW_HOSPIT",
  "PAS_ARRIVEE_MEDIC", "PAS_MIN_TRANSPORT", "PAS_HOSPIT", "FC_ARRIVEE_MEDIC",
  "REA_CATECHO", "ISS_SCORE", "AIS_HEAD_NECK", "TRISS", "DECES",
  "CGR_6H", "HEMOCUE_INITIAL", "HEMOCUE_ARRIVEE_HOP", "TOTAL_IGS",
  "LACTATES_HOSPIT", "PLAQUETTES", "FIBRINOGENE",
  "ANOMAL_PUPIL_PREHOSP", "CHOC_HEMORRAG_V1", "ASAP_PS", "ACR_PREHOSP",
  "PIC", "TRAUMA_CRANIEN", "ACIDE_TRANEXAMIQUE", "IOT_PREHOSP", "BLOC_RI_24H")

VARS_OPTIONNELLES <- c("CRITERE_INCLUSION", "AGE_PED", "AIS_HEAD_NECK")

manquantes <- setdiff(c(VARS_ID, VARS_ANALYSE), names(raw))
manq_opt   <- intersect(manquantes, VARS_OPTIONNELLES)
manq_oblig <- setdiff(manquantes, VARS_OPTIONNELLES)

if (length(manq_oblig))
  stop("REQUIRED variables missing: ", paste(manq_oblig, collapse = ", "))
for (v in manq_opt) {
  cat(sprintf("[1] Optional variable missing, created empty: %s\n", v))
  raw[[v]] <- NA_character_
}

raw <- raw[, c(VARS_ID, VARS_ANALYSE)]

VARS_NUM <- c("AGE", "AGE_PED", "ANNEE", "GLASGOW_INITIAL", "GLASGOW_HOSPIT",
              "PAS_ARRIVEE_MEDIC", "PAS_MIN_TRANSPORT", "PAS_HOSPIT",
              "FC_ARRIVEE_MEDIC", "ISS_SCORE", "AIS_HEAD_NECK", "TRISS",
              "CGR_6H", "HEMOCUE_INITIAL", "HEMOCUE_ARRIVEE_HOP", "TOTAL_IGS",
              "LACTATES_HOSPIT", "PLAQUETTES", "FIBRINOGENE")
VARS_CAT <- setdiff(VARS_ANALYSE, VARS_NUM)

# ---- 1.1 Status of each cell (before cleaning) -------------------------------
NIV_STATUT <- c("Observed", "NA", "NR", "ND", "Empty")
statut <- function(x) {
  x <- trimws(x)
  case_when(is.na(x) | x == ""                    ~ "Empty",
            x == "NR"                             ~ "NR",
            x == "NA"                             ~ "NA",
            x == "ND" | x %in% CODES_INCONNU      ~ "ND",
            TRUE                                  ~ "Observed")
}
st <- raw %>% mutate(across(all_of(VARS_ANALYSE), statut))

# ---- 1.2 Cleaning ------------------------------------------------------------
nettoyer <- function(x, v) {
  x <- trimws(as.character(x))
  a_harm <- x %in% names(HARMONISER)
  x[a_harm] <- HARMONISER[x[a_harm]]
  if (!is.null(NA_VAUT[[v]])) x[x %in% "NA"] <- NA_VAUT[[v]]
  x[x %in% c("NA", "NR", "ND", "", CODES_INCONNU)] <- NA
  x
}
df <- raw
for (v in VARS_ANALYSE) df[[v]] <- nettoyer(df[[v]], v)
df <- df %>% mutate(across(all_of(VARS_NUM), num))

# ---- 1.3 Observed categories (to compare with the labels of STEP 0) ----------
cat(sprintf("\n[1] Categories of categorical variables (counts < %d masked)\n",
            SEUIL_CELLULE))
for (v in VARS_CAT) {
  cat(sprintf("\n-- %s --\n", v))
  print(masquer_effectifs(table(df[[v]], useNA = "ifany")))
}

verif <- function(v, attendues) if (!any(df[[v]] %in% attendues))
  stop(sprintf("%s: none of the expected categories (%s). Adapt STEP 0.",
               v, paste(attendues, collapse = " / ")))
verif("DECES", c(OUI, NON)); verif("SEXE", LIB_MASCULIN)
verif("MECANISME_CAUSE", MECA_PENETRANTS); verif("TRAUMA_CRANIEN", OUI)
verif("REA_CATECHO", OUI)

# ---- 1.4 Scale of the registry TRISS (probability of survival) ---------------
if (max(df$TRISS, na.rm = TRUE) > 1) {
  cat("\n[1] Registry TRISS expressed in % -> converted to a proportion\n")
  df$TRISS <- df$TRISS / 100
}

# ==============================================================================
# STEP 2 — TRISS RECOMPUTATION AND FLOW CHART
# ==============================================================================

df <- df %>% mutate(
  triss = TRISS,
  meca_pene = case_when(MECANISME_CAUSE %in% MECA_PENETRANTS ~ 1,
                        is.na(MECANISME_CAUSE) ~ NA_real_,
                        TRUE ~ 0),
  pas_triss = case_when(REA_CATECHO %in% OUI ~ 55,
                        TRUE ~ pmin(PAS_ARRIVEE_MEDIC, PAS_MIN_TRANSPORT,
                                    PAS_HOSPIT, na.rm = TRUE)),
  pas_triss = if_else(is.infinite(pas_triss), NA_real_, pas_triss),
  pas_rts = case_when(pas_triss >= 90 ~ 4, pas_triss >= 76 ~ 3,
                      pas_triss >= 50 ~ 2, pas_triss >= 1 ~ 1,
                      pas_triss == 0 ~ 0, TRUE ~ NA_real_),
  gcs_rts = case_when(GLASGOW_INITIAL >= 13 ~ 4, GLASGOW_INITIAL >= 9 ~ 3,
                      GLASGOW_INITIAL >= 6 ~ 2, GLASGOW_INITIAL %in% c(4, 5) ~ 1,
                      GLASGOW_INITIAL == 3 ~ 0, TRUE ~ NA_real_),
  fr_rts        = 4,
  rts_admission = 0.9368 * gcs_rts + 0.7326 * pas_rts + 0.2908 * fr_rts,
  age_index     = ifelse(AGE >= 55, 1, 0),
  b_triss = case_when(
    meca_pene == 1 ~ -0.6029 + 1.1430 * rts_admission - 0.1516 * ISS_SCORE
    - 2.6676 * age_index,
    meca_pene == 0 ~ -1.2470 + 0.9544 * rts_admission - 0.0768 * ISS_SCORE
    - 1.9052 * age_index,
    TRUE ~ NA_real_),
  triss_calc = 1 / (1 + exp(-b_triss)),
  triss2     = coalesce(triss, triss_calc),
  DECES_bin  = case_when(DECES == OUI ~ 1L, DECES == NON ~ 0L,
                         TRUE ~ NA_integer_))

r_dir <- cor(df$triss, df$triss_calc, use = "complete.obs")
cat(sprintf("\n[2] Correlation registry TRISS / recomputed Ps: %.3f\n", r_dir))
if (is.na(r_dir) || r_dir < 0.5)
  stop("Registry TRISS inconsistent with the recomputed value: check that it is ",
       "a probability of SURVIVAL.")

# ---- 2.1 Flow chart ----------------------------------------------------------
n0    <- nrow(df)
keep1 <- !duplicated(df[, VARS_ID])
keep2 <- keep1 & !(df$CRITERE_INCLUSION %in% NON)
pediatrique <- (!is.na(df$AGE) & df$AGE < AGE_MIN) |
  (is.na(df$AGE) & !is.na(df$AGE_PED))
keep3 <- keep2 & !pediatrique
n_comp_manq <- rowSums(is.na(df[, c("GLASGOW_INITIAL", "pas_triss",
                                    "ISS_SCORE", "AGE", "meca_pene")]))
insuffisant <- is.na(df$triss) & n_comp_manq > MAX_COMP_MANQ &
  is.na(df$DECES_bin)
keep4 <- keep3 & !insuffisant

n_site   <- table(df$SITE_ID[keep4])
petits   <- names(n_site)[n_site < MIN_N_CENTRE_INCL]
keep5    <- keep4 & !(df$SITE_ID %in% petits)

flow_table <- data.frame(
  etape = c("Records extracted", "After removal of duplicates",
            "After Traumabase inclusion criteria",
            sprintf("After exclusion of patients < %d years", AGE_MIN),
            "After exclusion (risk not computable AND unknown outcome)",
            sprintf("After exclusion of centers < %d patients (%d centers excluded)",
                    MIN_N_CENTRE_INCL, length(petits))),
  n = c(n0, sum(keep1), sum(keep2), sum(keep3), sum(keep4), sum(keep5)))
flow_table$exclus <- c(NA, -diff(flow_table$n))

cat("\n=== TABLE S0a: flow chart ===\n")
print(flow_table, row.names = FALSE, right = FALSE)

df <- df[keep5, ]
st <- st[keep5, ]

# ---- 2.2 Center anonymization ------------------------------------------------
sites  <- sort(unique(df$SITE_ID))
codes  <- sprintf("C%02d", sample(seq_along(sites)))
write.csv2(data.frame(SITE_ID = sites, centre = codes),
           file.path(LOCAL_DIR, "center_key_CONFIDENTIAL.csv"),
           row.names = FALSE)
df$centre <- factor(codes[match(df$SITE_ID, sites)], levels = sort(codes))
df$SITE_ID <- NULL; df$NUM_PATIENT <- NULL
st$SITE_ID <- NULL; st$NUM_PATIENT <- NULL

eff_c <- as.numeric(table(df$centre))
desc_centres_global <- data.frame(
  n_centres = nlevels(df$centre), n_patients = nrow(df),
  patients_par_centre_mediane = median(eff_c),
  q25 = unname(quantile(eff_c, .25)), q75 = unname(quantile(eff_c, .75)),
  min = min(eff_c), max = max(eff_c))
cat("\n=== Contributing centers ===\n")
print(desc_centres_global, row.names = FALSE)

activite <- df %>% filter(!is.na(ANNEE)) %>% group_by(centre) %>%
  summarise(annees = max(ANNEE) - min(ANNEE) + 1, n = n(), .groups = "drop") %>%
  mutate(patients_par_an = n / annees)
q3 <- function(x) round(unname(quantile(x, c(.5, .25, .75))), 1)
desc_activite <- data.frame(
  indicateur = c("Years of activity per center", "Patients per center per year"),
  mediane = c(q3(activite$annees)[1], q3(activite$patients_par_an)[1]),
  q25     = c(q3(activite$annees)[2], q3(activite$patients_par_an)[2]),
  q75     = c(q3(activite$annees)[3], q3(activite$patients_par_an)[3]))
cat("\n=== Center activity (median, IQR) ===\n")
print(desc_activite, row.names = FALSE, right = FALSE)

# ---- 2.3 Data status (NR / NA / ND) ------------------------------------------
table_statut <- do.call(rbind, lapply(VARS_ANALYSE, function(v) {
  pct <- round(100 * prop.table(table(factor(st[[v]], levels = NIV_STATUT))), 1)
  data.frame(variable = v,
             as.list(setNames(as.numeric(pct), paste0("pct_", NIV_STATUT))),
             check.names = FALSE)
}))
cat("\n=== TABLE S0b: data status by variable (%) ===\n")
print(table_statut, row.names = FALSE, right = FALSE)

nr_rate <- setNames(table_statut$pct_NR / 100, table_statut$variable)

manq_centres <- df %>% group_by(centre) %>%
  summarise(n = n(),
            pct_triss_registre = round(100 * mean(!is.na(triss)), 1),
            across(c(GLASGOW_INITIAL, pas_triss, ISS_SCORE, AGE, meca_pene,
                     DECES_bin),
                   ~ round(100 * mean(is.na(.x)), 1), .names = "manq_{.col}"),
            .groups = "drop") %>%
  filter(n >= SEUIL_CELLULE)
cat("\n=== TABLE S0c: missing data by center (%) ===\n")
print(as.data.frame(manq_centres), row.names = FALSE, right = FALSE)

# ==============================================================================
# STEP 3 — MULTIPLE IMPUTATION OF TRISS COMPONENTS (MICE)
# ==============================================================================

num_aux <- c("GLASGOW_HOSPIT", "CGR_6H", "HEMOCUE_INITIAL", "TOTAL_IGS",
             "LACTATES_HOSPIT", "PLAQUETTES")
cat_aux <- c("SEXE", "ANOMAL_PUPIL_PREHOSP", "CHOC_HEMORRAG_V1", "REA_CATECHO",
             "ASAP_PS", "ACR_PREHOSP", "PIC", "TRAUMA_CRANIEN")

aux_exclus <- c(num_aux, cat_aux)[nr_rate[c(num_aux, cat_aux)] > SEUIL_NR]
if (length(aux_exclus))
  cat("\n[3] Auxiliary variables excluded (NR > threshold):",
      paste(aux_exclus, collapse = ", "), "\n")

comp_nr <- nr_rate[c("GLASGOW_INITIAL", "PAS_HOSPIT", "ISS_SCORE", "AGE",
                     "MECANISME_CAUSE")]
if (any(comp_nr > 0.05))
  warning("TRISS components with > 5% NR: ",
          paste(names(comp_nr)[comp_nr > 0.05], collapse = ", "))

num_vars <- c("GLASGOW_INITIAL", "pas_triss", "ISS_SCORE", "AGE", "ANNEE",
              setdiff(num_aux, aux_exclus))
cat_vars <- c("MECANISME_CAUSE", "DECES_bin", "centre", "BASE_SOURCE",
              setdiff(cat_aux, aux_exclus))

imp_data <- df %>%
  select(all_of(c(num_vars, cat_vars))) %>%
  mutate(across(all_of(cat_vars), ~ factor(as.character(.x)))) %>%
  as.data.frame()

cat("\n[3] Missing data before imputation (%):\n")
print(round(colMeans(is.na(imp_data)) * 100, 1))

stopifnot(nrow(imp_data) == nrow(df))

if (N_CORES > 1 && "futuremice" %in% getNamespaceExports("mice")) {
  imp <- futuremice(imp_data, m = M, maxit = MAXIT, parallelseed = SEED,
                    n.core = min(N_CORES, M))
} else {
  imp <- mice(imp_data, m = M, maxit = MAXIT, seed = SEED, printFlag = FALSE)
}

print(imp$loggedEvents)
print(plot(imp))
print(densityplot(imp, ~ GLASGOW_INITIAL + ISS_SCORE + pas_triss))

# ---- 3.1 TRISS recomputation in each imputed dataset -------------------------

long <- complete(imp, action = "long", include = FALSE) %>%
  mutate(
    id_row    = as.integer(as.character(.id)),
    meca_pene = if_else(MECANISME_CAUSE %in% MECA_PENETRANTS, 1, 0),
    pas_rts = case_when(pas_triss >= 90 ~ 4, pas_triss >= 76 ~ 3,
                        pas_triss >= 50 ~ 2, pas_triss >= 1 ~ 1,
                        pas_triss == 0 ~ 0, TRUE ~ NA_real_),
    gcs_rts = case_when(GLASGOW_INITIAL >= 13 ~ 4, GLASGOW_INITIAL >= 9 ~ 3,
                        GLASGOW_INITIAL >= 6 ~ 2, GLASGOW_INITIAL %in% c(4, 5) ~ 1,
                        GLASGOW_INITIAL == 3 ~ 0, TRUE ~ NA_real_),
    rts     = 0.9368 * gcs_rts + 0.7326 * pas_rts + 0.2908 * 4,
    age_idx = ifelse(AGE >= 55, 1, 0),
    b = ifelse(meca_pene == 1,
               -0.6029 + 1.1430 * rts - 0.1516 * ISS_SCORE - 2.6676 * age_idx,
               -1.2470 + 0.9544 * rts - 0.0768 * ISS_SCORE - 1.9052 * age_idx),
    triss_calc_imp = 1 / (1 + exp(-b))) %>%
  arrange(.imp, id_row)

stopifnot(all(long$id_row == rep(seq_len(nrow(df)), M)))

long <- long %>% mutate(
  triss_registre = df$triss[id_row],
  triss_final    = coalesce(triss_registre, triss_calc_imp),
  risk_death     = 1 - triss_final,
  centre         = df$centre[id_row])

risk_mat <- matrix(long$risk_death, nrow = nrow(df), ncol = M)

# ==============================================================================
# STEP 4 — VALIDATION OF THE RECOMPUTED RISK
# ==============================================================================

auc_fast <- function(y, p) {
  r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

val <- df %>% filter(!is.na(triss), !is.na(triss_calc))
chk <- df %>% mutate(risk = 1 - triss2) %>%
  filter(!is.na(risk), DECES_bin %in% c(0, 1))

AUC_EMPIRIQUE <- auc_fast(chk$DECES_bin, chk$risk)
BRIER         <- mean((chk$risk - chk$DECES_bin)^2)
MORT_OBS      <- mean(chk$DECES_bin)
OE_RATIO      <- MORT_OBS / mean(chk$risk)

table_validation <- data.frame(
  indicateur = c("n (recomputed vs registry)", "Pearson correlation",
                 "Mean bias", "n (vs observed mortality)", "AUC",
                 "Brier score", "Observed-to-expected ratio",
                 "Observed mortality (%)", "Predicted mortality (%)"),
  valeur = c(nrow(val),
             round(cor(val$triss, val$triss_calc), 3),
             round(mean(val$triss_calc - val$triss), 3),
             nrow(chk), round(AUC_EMPIRIQUE, 3), round(BRIER, 3),
             round(OE_RATIO, 2), round(100 * MORT_OBS, 1),
             round(100 * mean(chk$risk), 1)))

cat("\n=== TABLE S1: validation of the recomputed risk ===\n")
print(table_validation, row.names = FALSE, right = FALSE)

fig_S1 <- ggplot(val, aes(triss, triss_calc)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  stat_bin_2d(aes(fill = after_stat(ifelse(count < SEUIL_CELLULE, NA, count))),
              bins = 50) +
  scale_fill_gradient(low = "grey85", high = "black", trans = "log10",
                      na.value = "transparent") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Registry TRISS probability of survival", y = "Recomputed probability of survival") +
  theme_journal() + theme(legend.position = "right")
print(fig_S1)

# ==============================================================================
# STEP 5 — RISK POOLS (SUBGROUPS)
# ==============================================================================

long$cgr6h_obs <- df$CGR_6H[long$id_row]
long$meca_obs  <- df$meca_pene[long$id_row]
long$tc_obs    <- df$TRAUMA_CRANIEN[long$id_row] %in% OUI

membre <- function(mask) long$risk_death[mask & !is.na(long$risk_death)]

pools <- list(
  "Overall" = long$risk_death[!is.na(long$risk_death)],
  "Hemorrhagic shock"   = membre(long$cgr6h_obs >= 4 & !is.na(long$cgr6h_obs)),
  "Penetrating" = membre(long$meca_obs == 1 & !is.na(long$meca_obs)),
  "Blunt"     = membre(long$meca_obs == 0 & !is.na(long$meca_obs)),
  "Traumatic brain injury" = membre(long$tc_obs))

pools <- lapply(pools, function(x) x[is.finite(x)])

table_pools <- data.frame(
  sous_groupe  = names(pools),
  n_patients   = sapply(pools, function(x) round(length(x) / M)),
  risque_moyen = round(sapply(pools, mean), 3),
  ecart_type   = round(sapply(pools, sd), 3),
  q25          = round(sapply(pools, quantile, .25), 3),
  mediane      = round(sapply(pools, median), 3),
  q75          = round(sapply(pools, quantile, .75), 3),
  row.names    = NULL)

cat("\n=== TABLE S2: populations and risk distributions ===\n")
print(table_pools, row.names = FALSE, right = FALSE)

# ==============================================================================
# STEP 6 — TABLE 1: PATIENT CHARACTERISTICS
# ==============================================================================

tab <- df %>% mutate(
  SEX_male     = ifelse(is.na(SEXE), NA, SEXE == LIB_MASCULIN),
  PAS_preH     = PAS_ARRIVEE_MEDIC,
  FC_preH      = FC_ARRIVEE_MEDIC,
  GCS          = GLASGOW_INITIAL,
  TXA_preH     = oui(ACIDE_TRANEXAMIQUE),
  IOT_preH     = oui(IOT_PREHOSP),
  Hb           = HEMOCUE_ARRIVEE_HOP,
  Lactates     = LACTATES_HOSPIT,
  Plaquettes   = PLAQUETTES / 1000,
  Fibrinogene  = FIBRINOGENE,
  AIS_head_gt2 = AIS_HEAD_NECK > 2,
  ISS          = ISS_SCORE,
  ISS_gt15     = ISS_SCORE > 15,
  Fib_lt1.5    = FIBRINOGENE < 1.5,
  RTS          = rts_admission,
  Bloc_H24     = replace_na(BLOC_RI_24H, NON) == OUI,
  Hemorragie   = oui(CHOC_HEMORRAG_V1),
  CGR_ge4      = CGR_6H >= 4)

grps <- list(
  "Overall" = rep(TRUE, nrow(tab)),
  "Hemorrhagic shock"   = tab$CGR_ge4 %in% TRUE,
  "Penetrating" = tab$meca_pene %in% 1,
  "Blunt"     = tab$meca_pene %in% 0,
  "Traumatic brain injury" = tab$TRAUMA_CRANIEN %in% OUI)

cont_vars   <- c("AGE", "PAS_preH", "FC_preH", "GCS", "Hb", "Lactates",
                 "Plaquettes", "Fibrinogene", "ISS", "RTS")
cont_labels <- c("Age (years)", "Prehospital SBP (mmHg)", "Prehospital heart rate (bpm)",
                 "Initial Glasgow Coma Scale", "Hemoglobin (g/dL)", "Lactates (mmol/L)",
                 "Platelets (G/L)", "Fibrinogen (g/L)", "ISS", "RTS")
bin_vars    <- c("SEX_male", "meca_pene", "TXA_preH", "IOT_preH",
                 "AIS_head_gt2", "ISS_gt15", "Fib_lt1.5", "Bloc_H24",
                 "Hemorragie", "CGR_ge4")
bin_labels  <- c("Male sex", "Penetrating", "Prehospital tranexamic acid",
                 "Prehospital intubation", "Head AIS > 2", "ISS > 15",
                 "Fibrinogen < 1.5 g/L", "Surgery within 24 h", "Hemorrhagic shock (registry)",
                 "≥ 4 RBC units within 6 h")

vide <- function(v) all(is.na(tab[[v]]))
garder_c <- !vapply(cont_vars, vide, logical(1))
garder_b <- !vapply(bin_vars,  vide, logical(1))
cont_vars <- cont_vars[garder_c]; cont_labels <- cont_labels[garder_c]
bin_vars  <- bin_vars[garder_b];  bin_labels  <- bin_labels[garder_b]

fmt_cont <- function(x) { x <- x[!is.na(x)]
if (!length(x)) return("—")
sprintf("%.1f [%.1f–%.1f]", median(x), quantile(x, .25), quantile(x, .75)) }
fmt_bin  <- function(x) { x <- x[!is.na(x)]; k <- sum(x)
if (!length(x)) return("—")
if (k > 0 && k < SEUIL_CELLULE) return(sprintf("<%d", SEUIL_CELLULE))
sprintf("%d (%.1f%%)", k, 100 * mean(x)) }
fmt_miss <- function(x) { n <- sum(is.na(x))
if (n == 0) "—" else sprintf("%d (%.1f%%)", n, 100 * n / length(x)) }

rows <- list(
  n = c(Variable = "n",
        sapply(grps, function(m) as.character(sum(m))), Missing = ""),
  centres = c(Variable = "Contributing centers",
              sapply(grps, function(m) as.character(n_distinct(tab$centre[m]))),
              Missing = ""))
for (i in seq_along(cont_vars))
  rows[[cont_vars[i]]] <- c(Variable = cont_labels[i],
                            sapply(grps, function(m) fmt_cont(tab[[cont_vars[i]]][m])),
                            Missing = fmt_miss(tab[[cont_vars[i]]]))
for (i in seq_along(bin_vars))
  rows[[bin_vars[i]]] <- c(Variable = bin_labels[i],
                           sapply(grps, function(m) fmt_bin(tab[[bin_vars[i]]][m])),
                           Missing = fmt_miss(tab[[bin_vars[i]]]))

table1 <- as.data.frame(do.call(rbind, rows))
names(table1) <- c("Variable", names(grps), "Missing (whole population)")

cat("\n=== TABLE 1: patient characteristics ===\n")
print(table1, right = FALSE, row.names = FALSE)

# ==============================================================================
# STEP 7 — CALCULATION ENGINE
# ==============================================================================

# ---- 7.1 Effect models -------------------------------------------------------
apply_OR  <- function(p, OR)  expit(logit(clamp(p)) + log(OR))
apply_RR  <- function(p, RR)  pmin(RR * p, 1)
apply_ARR <- function(p, ARR) pmax(p - ARR, 0)

# ---- 7.2 Quadrature ----------------------------------------------------------
Z_QUAD <- qnorm((seq_len(K_QUAD) - 0.5) / K_QUAD)

pool_logit <- function(pool, Q = Q_QUAD) {
  n <- min(Q, length(pool))
  logit(clamp(unname(quantile(pool, probs = (seq_len(n) - 0.5) / n, type = 1))))
}

marg_risk <- function(lp, delta, sigma, effect = NULL) {
  if (sigma <= 0) {
    p <- expit(lp + delta)
    return(mean(if (is.null(effect)) p else effect(p)))
  }
  acc <- 0
  for (z in Z_QUAD) {
    p   <- expit(lp + delta + sigma * z)
    acc <- acc + mean(if (is.null(effect)) p else effect(p))
  }
  acc / K_QUAD
}

# ---- 7.3 Recentering and calibration -----------------------------------------
solve_delta <- function(lp, sigma, cible) {
  if (sigma <= 0) return(0)
  uniroot(function(d) marg_risk(lp, d, sigma) - cible, c(-6, 6))$root
}

auc_oracle <- function(pool, n = 20000L, seed = SEED) {
  set.seed(seed)
  p <- sample(pool, n, TRUE)
  y <- as.integer(runif(n) < p)
  if (sum(y) == 0 || sum(y) == n) return(NA_real_)
  auc_fast(y, p)
}

calibrate_sigma <- function(pool, target_auc = AUC_CIBLE, n = 20000L,
                            seed = SEED) {
  if (is.na(target_auc) || target_auc <= 0) return(0)
  a_max <- auc_oracle(pool, n, seed)
  if (is.na(a_max) || target_auc >= a_max) {
    warning(sprintf("Target AUC %.3f >= DGP ceiling %.3f -> sigma = 0",
                    target_auc, a_max))
    return(0)
  }
  set.seed(seed)
  p    <- sample(pool, n, TRUE)
  lp   <- logit(clamp(p))
  z    <- rnorm(n)
  u    <- runif(n)
  pbar <- mean(p)
  auc_at <- function(sigma) {
    d <- solve_delta(lp, sigma, pbar)
    y <- as.integer(u < expit(lp + d + sigma * z))
    if (sum(y) == 0 || sum(y) == n) return(NA_real_)
    auc_fast(y, p)
  }
  a_lo <- auc_at(0.05)
  if (is.na(a_lo) || a_lo <= target_auc) {
    warning(sprintf("Target AUC %.3f indistinguishable from the ceiling (AUC at sigma = 0.05: %.3f) -> sigma = 0",
                    target_auc, a_lo))
    return(0)
  }
  uniroot(function(s) auc_at(s) - target_auc, c(0.05, 6))$root
}

# ---- 7.4 Formula and unadjusted sample sizes ---------------------------------
nsn_formule <- function(pbar, OR, power = POWER, alpha = ALPHA)
  ceiling(power.prop.test(p1 = pbar, p2 = apply_OR(pbar, OR),
                          power = power, sig.level = alpha)$n)

nsn_unadj <- function(lp, delta, sigma, effect, power = POWER, alpha = ALPHA) {
  p1 <- marg_risk(lp, delta, sigma)
  p2 <- marg_risk(lp, delta, sigma, effect)
  ceiling(power.prop.test(p1 = p1, p2 = p2, power = power,
                          sig.level = alpha)$n)
}

# ---- 7.5 Simulated power -----------------------------------------------------
sim_trial <- function(pool, effect, N, sigma, delta) {
  p   <- sample(pool, 2 * N, TRUE)
  lp  <- logit(clamp(p))
  trt <- rep(0:1, each = N)
  ps  <- if (sigma > 0) expit(lp + delta + sigma * rnorm(2 * N))
  else expit(lp + delta)
  list(y   = rbinom(2 * N, 1, ifelse(trt == 1, effect(ps), ps)),
       trt = trt, lp = lp)
}

power_unadj_sim <- function(pool, effect, N, nsim, sigma, delta,
                            alpha = ALPHA) {
  hits <- 0L
  for (s in seq_len(nsim)) {
    d  <- sim_trial(pool, effect, N, sigma, delta)
    pv <- suppressWarnings(prop.test(
      c(sum(d$y[d$trt == 1]), sum(d$y[d$trt == 0])), c(N, N),
      correct = FALSE)$p.value)
    if (!is.na(pv) && pv < alpha) hits <- hits + 1L
  }
  hits / nsim
}

fit_trt <- function(y, trt, lp, spec = c("libre", "offset"), robust = FALSE) {
  spec <- match.arg(spec)
  X   <- if (spec == "libre") cbind(1, trt, lp) else cbind(1, trt)
  off <- if (spec == "libre") NULL else lp
  fit <- suppressWarnings(glm.fit(X, y, family = binomial(), offset = off,
                                  control = glm.control(epsilon = 1e-10,
                                                        maxit = 100)))
  if (!isTRUE(fit$converged)) return(c(b = NA, se = NA, p = NA))
  k  <- ncol(X)
  V  <- chol2inv(fit$qr$qr[seq_len(k), seq_len(k), drop = FALSE])
  se <- sqrt(V[2, 2])
  if (robust) {
    meat <- crossprod(X * (y - fit$fitted.values))
    se   <- sqrt((V %*% meat %*% V)[2, 2])
  }
  b <- unname(fit$coefficients[2])
  c(b = b, se = se, p = 2 * pnorm(abs(b / se), lower.tail = FALSE))
}

power_adj <- function(pool, effect, N, nsim, sigma, delta,
                      spec = "libre", robust = FALSE, alpha = ALPHA) {
  hits <- 0L
  for (s in seq_len(nsim)) {
    d  <- sim_trial(pool, effect, N, sigma, delta)
    pv <- fit_trt(d$y, d$trt, d$lp, spec, robust)["p"]
    if (!is.na(pv) && pv < alpha) hits <- hits + 1L
  }
  hits / nsim
}

# ---- 7.6 Bisection on N ------------------------------------------------------
find_nsn <- function(power_of_N, target = POWER, lo = 10, hi = 2e5,
                     seed = SEED) {
  pw <- function(N) { set.seed(seed); power_of_N(N) }
  if (pw(hi) < target) { warning("upper bound too low -> NA"); return(NA_integer_) }
  while (hi - lo > 1) {
    mid <- (lo + hi) %/% 2L
    if (pw(mid) >= target) hi <- mid else lo <- mid
  }
  as.integer(hi)
}

nsn_ci <- function(power_of_N, target = POWER, R = 20, hi = 2e5,
                   base_seed = SEED) {
  v <- vapply(seq_len(R), function(r)
    find_nsn(power_of_N, target = target, hi = hi, seed = base_seed + r),
    numeric(1))
  v <- v[is.finite(v)]
  c(mediane = median(v), lo = unname(quantile(v, .025)),
    hi = unname(quantile(v, .975)), moyenne = mean(v), ecart_type = sd(v))
}

# ==============================================================================
# STEP 8 — CALIBRATION OF THE DATA-GENERATING PROCESS
# ==============================================================================

# ---- 8.1 Observed discrimination by population -------------------------------
masques_obs <- list(
  "Overall" = rep(TRUE, nrow(df)),
  "Hemorrhagic shock"   = (df$CGR_6H >= 4) %in% TRUE,
  "Penetrating" = df$meca_pene %in% 1,
  "Blunt"     = df$meca_pene %in% 0,
  "Traumatic brain injury" = df$TRAUMA_CRANIEN %in% OUI)

table_auc_sg <- do.call(rbind, lapply(names(masques_obs), function(nm) {
  d <- df %>% mutate(risk = 1 - triss2) %>%
    filter(masques_obs[[nm]], !is.na(risk), DECES_bin %in% c(0, 1))
  data.frame(sous_groupe = nm, n = nrow(d), deces = sum(d$DECES_bin),
             AUC_observee = round(auc_fast(d$DECES_bin, d$risk), 3),
             AUC_plafond_DGP = round(auc_oracle(pools[[nm]]), 3),
             OE = round(mean(d$DECES_bin) / mean(d$risk), 2),
             R2_brier = round(1 - mean((d$risk - d$DECES_bin)^2) /
                                (mean(d$DECES_bin) * (1 - mean(d$DECES_bin))), 3),
             stringsAsFactors = FALSE)
}))

cat("\n=== TABLE S3bis: observed discrimination by population ===\n")
print(table_auc_sg, row.names = FALSE, right = FALSE)

# ---- 8.2 Calibration of sigma and delta --------------------------------------
SIGMA_EPS <- calibrate_sigma(pools[["Overall"]], AUC_CIBLE)

lp_pool <- lapply(pools, pool_logit)
pbar    <- sapply(pools, mean)
delta   <- mapply(function(lp, m) solve_delta(lp, SIGMA_EPS, m), lp_pool, pbar)

table_dgp <- data.frame(
  sous_groupe    = names(pools),
  risque_moyen   = round(pbar, 4),
  delta          = round(delta, 4),
  moyenne_p_star = round(mapply(function(lp, d)
    marg_risk(lp, d, SIGMA_EPS), lp_pool, delta), 4),
  row.names = NULL)

cat("\n=== TABLE S3: calibration of the DGP ===\n")
cat(sprintf("In-sample AUC = %.3f | target AUC = %.3f | SIGMA_EPS = %.3f\n",
            AUC_EMPIRIQUE, AUC_CIBLE, SIGMA_EPS))
print(table_dgp, row.names = FALSE, right = FALSE)
stopifnot(max(abs(table_dgp$moyenne_p_star - table_dgp$risque_moyen)) < 1e-3)

# ==============================================================================
# STEP 9 — VALIDATION: ANALYTICAL vs SIMULATED UNADJUSTED SAMPLE SIZE
# ==============================================================================

set.seed(SEED)
table_validation_analytique <- do.call(rbind, lapply(names(pools), function(nm) {
  n_ana <- nsn_unadj(lp_pool[[nm]], delta[nm], SIGMA_EPS,
                     function(p) apply_OR(p, OR_REF))
  n_sim <- find_nsn(function(N)
    power_unadj_sim(pools[[nm]], function(p) apply_OR(p, OR_REF), N,
                    NSIM_VAL, SIGMA_EPS, delta[nm]), hi = 3e5)
  data.frame(sous_groupe = nm, NSN_analytique = n_ana, NSN_simule = n_sim,
             ecart_relatif = round(n_sim / n_ana - 1, 3))
}))

cat("\n=== TABLE S4: validation of the analytical result (OR = 0.70) ===\n")
print(table_validation_analytique, row.names = FALSE, right = FALSE)

# ==============================================================================
# STEP 10 — MAIN SAMPLE-SIZE CALCULATION (single source of truth)
# ==============================================================================

OR_grid <- sort(unique(c(OR_REF, seq(0.50, 0.90, by = 0.05))), decreasing = TRUE)

set.seed(SEED)
nsn_master <- do.call(rbind, lapply(names(pools), function(nm) {
  cat(sprintf("  [master] %-22s\n", nm))
  lp <- lp_pool[[nm]]; d <- delta[nm]; p1 <- pbar[nm]
  do.call(rbind, lapply(OR_grid, function(or) {
    eff <- function(p) apply_OR(p, or)
    data.frame(
      sous_groupe    = nm,
      OR             = or,
      n_patients     = round(length(pools[[nm]]) / M),
      risque_moyen   = round(p1, 4),
      ARR_formule    = round(p1 - apply_OR(p1, or), 4),
      ARR_marginal   = round(marg_risk(lp, d, SIGMA_EPS)
                             - marg_risk(lp, d, SIGMA_EPS, eff), 4),
      NSN_formule    = nsn_formule(p1, or),
      NSN_non_ajuste = nsn_unadj(lp, d, SIGMA_EPS, eff),
      NSN_ajuste     = if (isTRUE(all.equal(or, OR_REF)))
        find_nsn(function(N) power_adj(pools[[nm]], eff, N, NSIM_ADJ,
                                       SIGMA_EPS, d, "libre"), hi = 3e4)
      else NA_integer_,
      stringsAsFactors = FALSE)
  }))
}))

if (anyNA(nsn_master$NSN_non_ajuste) ||
    any(is.na(nsn_master$NSN_ajuste) & abs(nsn_master$OR - OR_REF) < 1e-9))
  warning("nsn_master: sample size not found -> increase the upper bound.")

# ==============================================================================
# STEP 11 — ANALYSIS 1: FORMULA vs UNADJUSTED vs ADJUSTED (OR = 0.70)
# ==============================================================================

res1 <- nsn_master %>%
  filter(abs(OR - OR_REF) < 1e-9) %>%
  mutate(perte_ARR        = round(ARR_marginal / ARR_formule, 2),
         ratio_non_ajuste = round(NSN_non_ajuste / NSN_formule, 2),
         ratio_ajuste     = round(NSN_ajuste / NSN_formule, 2)) %>%
  select(sous_groupe, n_patients, risque_moyen, ARR_formule, ARR_marginal,
         perte_ARR, NSN_formule, NSN_non_ajuste, NSN_ajuste,
         ratio_non_ajuste, ratio_ajuste)

cat("\n=== TABLE 2: sample size by population (OR = 0.70) ===\n")
print(res1, row.names = FALSE, right = FALSE)

cat("\n--- Internal consistency: observed ratio vs (ARR_formula/ARR_marginal)^2 ---\n")
print(data.frame(
  sous_groupe = res1$sous_groupe,
  ratio_observe = res1$ratio_non_ajuste,
  ratio_attendu = round((res1$ARR_formule / res1$ARR_marginal)^2, 2)),
  row.names = FALSE, right = FALSE)

# ---- 11.1 Predictive index: R² of baseline risk ------------------------------
indice_R2 <- function(pool, lp, d, sigma) {
  pb <- mean(pool)
  m2 <- marg_risk(lp, d, sigma, function(p) p^2)
  c(R2_score  = var(pool) / (pb * (1 - pb)),
    R2_latent = (m2 - pb^2) / (pb * (1 - pb)))
}

table_R2 <- do.call(rbind, lapply(names(pools), function(nm) {
  lp <- lp_pool[[nm]]; d <- delta[nm]
  ix <- indice_R2(pools[[nm]], lp, d, SIGMA_EPS)
  p1 <- marg_risk(lp, d, SIGMA_EPS)
  p2 <- marg_risk(lp, d, SIGMA_EPS, function(p) apply_OR(p, OR_REF))
  r  <- res1[res1$sous_groupe == nm, ]
  data.frame(
    sous_groupe          = nm,
    R2_score             = round(unname(ix["R2_score"]), 3),
    R2_latent            = round(unname(ix["R2_latent"]), 3),
    OR_conditionnel      = OR_REF,
    OR_marginal          = round((p2 / (1 - p2)) / (p1 / (1 - p1)), 3),
    ratio_observe        = r$ratio_non_ajuste,
    ratio_predit_score   = round(1 / (1 - unname(ix["R2_score"]))^2, 2),
    ratio_predit_latent  = round(1 / (1 - unname(ix["R2_latent"]))^2, 2),
    ratio_ajuste_observe = r$ratio_ajuste,
    stringsAsFactors = FALSE)
}))

cat("\n=== TABLE 2bis: R² index and predicted ratio (OR = 0.70) ===\n")
print(table_R2, row.names = FALSE, right = FALSE)

z_a <- qnorm(1 - ALPHA / 2); z_b <- qnorm(POWER)
puiss_approx <- function(ratio) round(pnorm(sqrt(1 / ratio) * (z_a + z_b) - z_a), 3)
table_puissance <- res1 %>%
  transmute(sous_groupe, NSN_formule,
            puissance_non_ajustee = puiss_approx(ratio_non_ajuste),
            puissance_ajustee     = puiss_approx(ratio_ajuste))
cat("\n=== Approximate power of formula-sized trials ===\n")
print(table_puissance, row.names = FALSE, right = FALSE)

if (MODE == "FINAL") {
  cat("\n--- Monte Carlo range of the adjusted sample size (whole population) ---\n")
  print(round(nsn_ci(function(N)
    power_adj(pools[["Overall"]],
              function(p) apply_OR(p, OR_REF), N, NSIM_ADJ,
              SIGMA_EPS, delta["Overall"], "libre"), hi = 3e4), 1))
}

# ==============================================================================
# STEP 12 — TYPE I ERROR AND EFFICIENCY OF ESTIMATORS
# ==============================================================================

N_typeI  <- res1$NSN_non_ajuste[res1$sous_groupe == "Overall"]
pool_gen <- pools[["Overall"]]
lp_gen   <- lp_pool[["Overall"]]
d_gen    <- delta["Overall"]
eff_nul  <- function(p) p

# ---- 12.1 Type I error under OR = 1 ------------------------------------------
methodes <- list(
  "Unadjusted"              = function() power_unadj_sim(
    pool_gen, eff_nul, N_typeI, NSIM_VAL, SIGMA_EPS, d_gen),
  "Adjusted, estimated coefficient"    = function() power_adj(
    pool_gen, eff_nul, N_typeI, NSIM_VAL, SIGMA_EPS, d_gen, "libre"),
  "Adjusted, estimated coefficient, robust variance" = function() power_adj(
    pool_gen, eff_nul, N_typeI, NSIM_VAL, SIGMA_EPS, d_gen, "libre",
    robust = TRUE),
  "Offset, model variance"  = function() power_adj(
    pool_gen, eff_nul, N_typeI, NSIM_VAL, SIGMA_EPS, d_gen, "offset"),
  "Offset, robust variance" = function() power_adj(
    pool_gen, eff_nul, N_typeI, NSIM_VAL, SIGMA_EPS, d_gen, "offset",
    robust = TRUE))

table_typeI <- data.frame(
  methode = names(methodes),
  erreur_type_I = sapply(methodes, function(f) { set.seed(SEED); round(f(), 3) }),
  cible = ALPHA, row.names = NULL)
mc_se <- sqrt(table_typeI$erreur_type_I * (1 - table_typeI$erreur_type_I) / NSIM_VAL)
table_typeI$IC95_bas  <- round(table_typeI$erreur_type_I - 1.96 * mc_se, 3)
table_typeI$IC95_haut <- round(table_typeI$erreur_type_I + 1.96 * mc_se, 3)
table_typeI <- table_typeI[, c("methode", "erreur_type_I", "IC95_bas",
                               "IC95_haut", "cible")]

cat(sprintf("\n=== TABLE S5: type I error (OR = 1, N = %d) ===\n", N_typeI))
print(table_typeI, row.names = FALSE, right = FALSE)

# ---- 12.2 Efficiency ---------------------------------------------------------
diag_estimateurs <- function(pool, N, sigma, delta, nsim, seed = SEED) {
  set.seed(seed)
  out <- vector("list", nsim)
  for (s in seq_len(nsim)) {
    d <- sim_trial(pool, eff_nul, N, sigma, delta)
    f_libre <- fit_trt(d$y, d$trt, d$lp, "libre")
    f_lrob  <- fit_trt(d$y, d$trt, d$lp, "libre", robust = TRUE)
    f_off   <- fit_trt(d$y, d$trt, d$lp, "offset")
    f_rob   <- fit_trt(d$y, d$trt, d$lp, "offset", robust = TRUE)
    out[[s]] <- data.frame(b_libre = f_libre["b"], se_libre = f_libre["se"],
                           se_libre_rob = f_lrob["se"],
                           b_off = f_off["b"], se_off = f_off["se"],
                           se_off_rob = f_rob["se"])
  }
  do.call(rbind, out)
}

dg <- diag_estimateurs(pool_gen, N_typeI, SIGMA_EPS, d_gen, NSIM_DIAG)

table_efficience <- data.frame(
  estimateur   = c("Estimated coefficient", "Estimated coefficient, robust variance",
                   "Offset, model variance", "Offset, robust variance"),
  biais        = round(c(mean(dg$b_libre), mean(dg$b_libre),
                         mean(dg$b_off), mean(dg$b_off)), 4),
  sd_empirique = round(c(sd(dg$b_libre), sd(dg$b_libre),
                         sd(dg$b_off), sd(dg$b_off)), 4),
  se_declaree  = round(c(mean(dg$se_libre), mean(dg$se_libre_rob),
                         mean(dg$se_off), mean(dg$se_off_rob)), 4),
  row.names = NULL)
table_efficience$ratio_sd_se <- round(
  table_efficience$sd_empirique / table_efficience$se_declaree, 2)
table_efficience$NSN_relatif <- round(
  (table_efficience$sd_empirique / table_efficience$sd_empirique[1])^2, 2)

cat("\n=== TABLE S6: efficiency and calibration of standard errors ===\n")
print(table_efficience, row.names = FALSE, right = FALSE)

# ==============================================================================
# STEP 13 — ANALYSIS 2: STABILITY OF THE RATIO ACROSS THE OR GRID
# ==============================================================================

res2 <- nsn_master %>%
  transmute(sous_groupe, OR, NSN_formule, NSN_non_ajuste,
            ratio = round(NSN_non_ajuste / NSN_formule, 2))

cat("\n=== TABLE S7: sample size and ratio according to the target odds ratio ===\n")
print(as.data.frame(res2), row.names = FALSE, right = FALSE)
cat(sprintf("Range of the ratio, all populations and effects: %.2f to %.2f\n",
            min(res2$ratio), max(res2$ratio)))

# ==============================================================================
# STEP 14 — ANALYSIS 3: EFFECT SCALE (OR vs RR vs ARR)
# ==============================================================================

res3 <- do.call(rbind, lapply(names(pools), function(nm) {
  lp <- lp_pool[[nm]]; d <- delta[nm]; p1 <- pbar[nm]
  g  <- apply_OR(p1, OR_REF); RR <- g / p1; ARR <- p1 - g
  m  <- nsn_master[nsn_master$sous_groupe == nm &
                     abs(nsn_master$OR - OR_REF) < 1e-9, ]
  data.frame(sous_groupe = nm, NSN_formule = m$NSN_formule,
             OR_constant  = m$NSN_non_ajuste,
             RR_constant  = nsn_unadj(lp, d, SIGMA_EPS, function(p) apply_RR(p, RR)),
             ARR_constant = nsn_unadj(lp, d, SIGMA_EPS, function(p) apply_ARR(p, ARR)),
             stringsAsFactors = FALSE)
})) %>%
  mutate(ratio_OR  = round(OR_constant / NSN_formule, 2),
         ratio_RR  = round(RR_constant / NSN_formule, 2),
         ratio_ARR = round(ARR_constant / NSN_formule, 2))

cat("\n=== TABLE 3: sample size according to the assumed effect scale ===\n")
print(res3, row.names = FALSE, right = FALSE)

cat("\n--- Ranking of populations by scale (1 = highest ratio) ---\n")
print(data.frame(
  sous_groupe = res3$sous_groupe,
  rang_OR  = rank(-res3$ratio_OR),
  rang_ARR = rank(-res3$ratio_ARR)), row.names = FALSE, right = FALSE)

# ==============================================================================
# STEP 15 — SENSITIVITY ANALYSES
# ==============================================================================

# ---- 15.1 Sensitivity to the target AUC (whole population) -------------------
set.seed(SEED)
sens_auc <- do.call(rbind, lapply(AUC_SENS, function(a) {
  s  <- calibrate_sigma(pool_gen, a)
  d  <- solve_delta(lp_gen, s, pbar["Overall"])
  nf <- nsn_formule(pbar["Overall"], OR_REF)
  nu <- nsn_unadj(lp_gen, d, s, function(p) apply_OR(p, OR_REF))
  na <- find_nsn(function(N) power_adj(pool_gen, function(p) apply_OR(p, OR_REF),
                                       N, NSIM_ADJ, s, d, "libre"), hi = 3e4)
  data.frame(scenario = ifelse(a == 0, "oracle (perfect score)",
                               sprintf("AUC %.2f", a)),
             sigma = round(s, 3), NSN_formule = nf,
             NSN_non_ajuste = nu, NSN_ajuste = na,
             ratio_non_ajuste = round(nu / nf, 2),
             ratio_ajuste     = round(na / nf, 2), stringsAsFactors = FALSE)
}))

cat("\n=== TABLE S8: sensitivity to the assumed discrimination of the score ===\n")
print(sens_auc, row.names = FALSE, right = FALSE)

# ---- 15.1bis Sigma calibrated in each population -----------------------------
auc_obs   <- setNames(table_auc_sg$AUC_observee, table_auc_sg$sous_groupe)
sigma_pop <- sapply(names(pools),
                    function(nm) calibrate_sigma(pools[[nm]], auc_obs[nm]))
delta_pop <- mapply(function(lp, s, m) solve_delta(lp, s, m),
                    lp_pool, sigma_pop, pbar)

set.seed(SEED)
sens_sigma_pop <- do.call(rbind, lapply(names(pools), function(nm) {
  lp <- lp_pool[[nm]]; s <- sigma_pop[nm]; d <- delta_pop[nm]
  nf <- nsn_formule(pbar[nm], OR_REF)
  nu <- nsn_unadj(lp, d, s, function(p) apply_OR(p, OR_REF))
  na <- find_nsn(function(N) power_adj(pools[[nm]],
                                       function(p) apply_OR(p, OR_REF), N, NSIM_ADJ, s, d, "libre"),
                 hi = 3e4)
  data.frame(sous_groupe = nm, AUC_observee = unname(auc_obs[nm]),
             sigma = round(s, 3), NSN_formule = nf,
             NSN_non_ajuste = nu, NSN_ajuste = na,
             ratio_non_ajuste = round(nu / nf, 2),
             ratio_ajuste     = round(na / nf, 2), stringsAsFactors = FALSE)
}))

cat("\n=== TABLE S10: sigma calibrated on the discrimination of each population ===\n")
print(sens_sigma_pop, row.names = FALSE, right = FALSE)

# ---- 15.2 Sensitivity to recalibration on observed mortality -----------------
sens_calib <- do.call(rbind, lapply(names(pools), function(nm) {
  lp <- lp_pool[[nm]]
  cible_cal <- pbar[nm] * OE_RATIO
  d_cal <- solve_delta(lp, SIGMA_EPS, cible_cal)
  nf <- nsn_formule(cible_cal, OR_REF)
  nu <- nsn_unadj(lp, d_cal, SIGMA_EPS, function(p) apply_OR(p, OR_REF))
  data.frame(sous_groupe = nm,
             risque_brut = round(pbar[nm], 4),
             risque_recale = round(cible_cal, 4),
             NSN_formule = nf, NSN_non_ajuste = nu,
             ratio = round(nu / nf, 2), stringsAsFactors = FALSE)
}))

cat(sprintf("\n=== TABLE S9: sensitivity to recalibration on observed mortality (O:E = %.2f) ===\n",
            OE_RATIO))
print(sens_calib, row.names = FALSE, right = FALSE)

# ==============================================================================
# STEP 16 — ANALYSIS 4: CONSISTENCY ACROSS CENTERS
# ------------------------------------------------------------------------------
# ==============================================================================

ratio_pool <- function(pool, sigma = SIGMA_EPS, or = OR_REF) {
  pool <- pool[is.finite(pool)]
  lp <- pool_logit(pool); pb <- mean(pool)
  d  <- solve_delta(lp, sigma, pb)
  nf <- nsn_formule(pb, or)
  nu <- nsn_unadj(lp, d, sigma, function(p) apply_OR(p, or))
  c(risque_moyen = pb, ecart_type = sd(pool), NSN_formule = nf,
    NSN_non_ajuste = nu, ratio = nu / nf)
}

boot_ratio <- function(ids, B, seed = SEED) {
  set.seed(seed)
  vapply(seq_len(B), function(b) {
    s <- sample(ids, length(ids), replace = TRUE)
    unname(ratio_pool(as.vector(risk_mat[s, , drop = FALSE]))["ratio"])
  }, numeric(1))
}

obs_c <- df %>% mutate(risk = 1 - triss2) %>%
  filter(!is.na(risk), DECES_bin %in% c(0, 1))

desc_c <- df %>% group_by(centre) %>%
  summarise(n = n(), deces = sum(DECES_bin, na.rm = TRUE), .groups = "drop")
centres_elig <- as.character(desc_c$centre[desc_c$n >= MIN_N_CENTRE &
                                             desc_c$deces >= MIN_DC_CENTRE])
cat(sprintf("\n[16] Centers analyzed: %d / %d (n >= %d and deaths >= %d)\n",
            length(centres_elig), nlevels(df$centre),
            MIN_N_CENTRE, MIN_DC_CENTRE))

# ---- 16.1 Ratio by center ----------------------------------------------------
res_centres <- do.call(rbind, lapply(centres_elig, function(cc) {
  cat(sprintf("  [centre] %s\n", cc))
  ids <- which(df$centre == cc)
  est <- ratio_pool(as.vector(risk_mat[ids, , drop = FALSE]))
  bt  <- boot_ratio(ids, B_CENTRE)
  o   <- obs_c[obs_c$centre == cc, ]
  pool_c <- as.vector(risk_mat[ids, , drop = FALSE])
  lp_c   <- pool_logit(pool_c)
  ix     <- indice_R2(pool_c, lp_c, solve_delta(lp_c, SIGMA_EPS, mean(pool_c)),
                      SIGMA_EPS)
  data.frame(
    centre         = cc,
    n              = length(ids),
    deces          = sum(df$DECES_bin[ids], na.rm = TRUE),
    mortalite_obs  = round(mean(df$DECES_bin[ids], na.rm = TRUE), 4),
    risque_moyen   = round(unname(est["risque_moyen"]), 4),
    ecart_type     = round(unname(est["ecart_type"]), 4),
    AUC            = round(auc_fast(o$DECES_bin, o$risk), 3),
    OE             = round(mean(o$DECES_bin) / mean(o$risk), 2),
    NSN_formule    = unname(est["NSN_formule"]),
    NSN_non_ajuste = unname(est["NSN_non_ajuste"]),
    ratio          = round(unname(est["ratio"]), 2),
    IC95_bas       = round(unname(quantile(bt, .025)), 2),
    IC95_haut      = round(unname(quantile(bt, .975)), 2),
    R2_score       = round(unname(ix["R2_score"]), 3),
    R2_latent      = round(unname(ix["R2_latent"]), 3),
    ratio_predit   = round(1 / (1 - unname(ix["R2_score"]))^2, 2),
    stringsAsFactors = FALSE)
}))

# ---- 16.2 Adjusted sample size by center -------------------------------------
nsn_adj_extrap <- function(pool, effect, N0, nsim, sigma, delta, seed,
                           power = POWER, alpha = ALPHA) {
  set.seed(seed)
  est <- t(vapply(seq_len(nsim), function(s) {
    d <- sim_trial(pool, effect, N0, sigma, delta)
    fit_trt(d$y, d$trt, d$lp, "libre")[c("b", "se")]
  }, numeric(2)))
  est <- est[is.finite(est[, 1]) & is.finite(est[, 2]), , drop = FALSE]
  b   <- abs(mean(est[, 1])); se0 <- mean(est[, 2]); sd0 <- sd(est[, 1])
  as.integer(ceiling(N0 * ((qnorm(1 - alpha / 2) * se0 +
                              qnorm(power) * sd0) / b)^2))
}

n_extrap_gen <- nsn_adj_extrap(
  pool_gen, function(p) apply_OR(p, OR_REF),
  N0 = res1$NSN_ajuste[res1$sous_groupe == "Overall"],
  nsim = NSIM_ADJ, sigma = SIGMA_EPS, delta = d_gen, seed = SEED)
cat(sprintf(paste0("\n[16.2] Check, whole population: adjusted sample size ",
                   "bisection = %d | extrapolation = %d (difference %+.1f %%)\n"),
            res1$NSN_ajuste[res1$sous_groupe == "Overall"],
            n_extrap_gen,
            100 * (n_extrap_gen /
                     res1$NSN_ajuste[res1$sous_groupe == "Overall"] - 1)))

if (RUN_ADJ_CENTRE) {
  calc_centre <- function(i) {
    cc     <- res_centres$centre[i]
    pool_c <- as.vector(risk_mat[df$centre == cc, , drop = FALSE])
    d_c    <- solve_delta(pool_logit(pool_c), SIGMA_EPS, mean(pool_c))
    gain   <- with(res1[res1$sous_groupe == "Overall", ],
                   NSN_ajuste / NSN_non_ajuste)
    N0     <- max(100L, as.integer(res_centres$NSN_non_ajuste[i] * gain))
    nsn_adj_extrap(pool_c, function(p) apply_OR(p, OR_REF), N0, NSIM_ADJ,
                   SIGMA_EPS, d_c, seed = SEED + i)
  }
  idx <- seq_len(nrow(res_centres))
  t0  <- Sys.time()
  res_centres$NSN_ajuste <- unlist(
    if (.Platform$OS.type == "unix" && N_CORES > 1)
      parallel::mclapply(idx, calc_centre, mc.cores = N_CORES)
    else lapply(idx, calc_centre))
  res_centres$ratio_ajuste <- round(res_centres$NSN_ajuste /
                                      res_centres$NSN_formule, 2)
  cat(sprintf("[16.2] Adjusted sample size by center computed in %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

cat("\n=== TABLE S11: sample size and ratio by center (OR = 0.70) ===\n")
print(res_centres, row.names = FALSE, right = FALSE)

# ---- 16.3 Summary ------------------------------------------------------------
est_global <- ratio_pool(as.vector(risk_mat))
bt_global  <- boot_ratio(seq_len(nrow(df)), B_CENTRE)

synthese_centres <- data.frame(
  indicateur = c("Centers (total)", "Centers analyzed",
                 "Patients in the centers analyzed (%)",
                 "Pooled ratio", "Bootstrap 95% CI of the pooled ratio",
                 "Ratio by center: median [IQR]",
                 "Ratio by center: min – max",
                 "Centers with a lower 95% CI bound > 1",
                 "Spearman correlation, ratio / mean risk of the center",
                 "Spearman correlation, ratio / O:E of the center",
                 "Spearman correlation, ratio / R² of the center",
                 "Median relative difference, predicted (R²) vs observed ratio [min – max]"),
  valeur = c(
    nlevels(df$centre), nrow(res_centres),
    round(100 * sum(res_centres$n) / nrow(df), 1),
    round(unname(est_global["ratio"]), 2),
    sprintf("%.2f – %.2f", quantile(bt_global, .025), quantile(bt_global, .975)),
    sprintf("%.2f [%.2f – %.2f]", median(res_centres$ratio),
            quantile(res_centres$ratio, .25), quantile(res_centres$ratio, .75)),
    sprintf("%.2f – %.2f", min(res_centres$ratio), max(res_centres$ratio)),
    sprintf("%d / %d", sum(res_centres$IC95_bas > 1), nrow(res_centres)),
    round(cor(res_centres$ratio, res_centres$risque_moyen, method = "spearman"), 2),
    round(cor(res_centres$ratio, res_centres$OE, method = "spearman"), 2),
    round(cor(res_centres$ratio, res_centres$R2_score, method = "spearman"), 2),
    with(res_centres, {
      e <- ratio_predit / ratio - 1
      sprintf("%+.1f %% [%+.1f – %+.1f %%]", 100 * median(e),
              100 * min(e), 100 * max(e))
    })))

cat("\n=== TABLE S12: summary of consistency across centers ===\n")
print(synthese_centres, row.names = FALSE, right = FALSE)

cat(sprintf("Check: pooled ratio (step 16) = %.2f | Table 2 = %.2f\n",
            est_global["ratio"],
            res1$ratio_non_ajuste[res1$sous_groupe == "Overall"]))

# ---- 16.4 Leave-one-center-out -----------------------------------------------
loco <- do.call(rbind, lapply(levels(df$centre), function(cc) {
  ids <- which(df$centre != cc)
  r   <- ratio_pool(as.vector(risk_mat[ids, , drop = FALSE]))
  data.frame(centre_exclu = cc, n_restant = length(ids),
             risque_moyen = round(unname(r["risque_moyen"]), 4),
             ratio = round(unname(r["ratio"]), 2))
}))

cat("\n=== TABLE S13: pooled ratio excluding each center in turn ===\n")
print(loco, row.names = FALSE, right = FALSE)
cat(sprintf("Leave-one-center-out range: %.2f to %.2f\n",
            min(loco$ratio), max(loco$ratio)))

# ==============================================================================
# STEP 17 — FIGURES
# ==============================================================================

pal <- c("Overall" = "#1c4e80", "Hemorrhagic shock"   = "#a30000",
         "Penetrating" = "#b8860b", "Blunt"     = "#009e73",
         "Traumatic brain injury" = "#7030a0")
grp_levels <- names(pools)

ECHELLE_RISQUE <- scale_x_log10(
  limits = c(1e-4, 1), breaks = c(1e-4, 1e-3, 1e-2, 0.1, 1),
  labels = c("0.01%", "0.1%", "1%", "10%", "100%"))

# ---- FIGURE 1: distribution of predicted risk --------------------------------
  stats_A <- data.frame(
    grp    = factor(grp_levels, levels = grp_levels),
    moy    = sapply(pools[grp_levels], mean),
    med    = sapply(pools[grp_levels], median))

set.seed(SEED)
fig1 <- do.call(rbind, lapply(grp_levels, function(nm)
  data.frame(risk = clamp(echant(pools[[nm]])), grp = nm))) %>%
  mutate(grp = factor(grp, levels = grp_levels)) %>%
  ggplot(aes(risk, fill = grp, color = grp)) +
  geom_density(alpha = .18, linewidth = .6) +
  geom_vline(data = stats_A, aes(xintercept = moy, color = grp),
             linetype = "dashed", linewidth = .5) +
  geom_vline(data = stats_A, aes(xintercept = med, color = grp),
             linetype = "dotted", linewidth = .6) +
  facet_wrap(~ grp, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = pal) + scale_color_manual(values = pal) +
  scale_x_continuous(trans = scales::logit_trans(),
                     breaks = c(1e-3, 1e-2, 0.1, 0.5, 0.9, 0.99),
                     labels = c("0.1%", "1%", "10%", "50%", "90%", "99%")) +
  coord_cartesian(xlim = c(1e-3, 0.995)) +
  labs(x = "Predicted risk of death (logit scale)", y = "Density") +
  theme_journal() + theme(legend.position = "none")
print(fig1)

# ---- FIGURE 2: absolute risk reduction according to baseline risk ------------
courbe <- data.frame(p = seq(0.001, 0.999, length.out = 600)) %>%
  mutate(arr = p - apply_OR(p, OR_REF))
points_grp <- data.frame(grp = factor(grp_levels, levels = grp_levels),
                         p = unname(pbar[grp_levels])) %>%
  mutate(arr = p - apply_OR(p, OR_REF))

fig2 <- ggplot(courbe, aes(p, arr)) +
  geom_line(linewidth = 1, color = "#4a4a4a") +
  geom_point(data = points_grp, aes(color = grp), size = 2.8) +
  scale_color_manual(values = pal, name = NULL) +
  labs(x = "Predicted risk of death (p)",
       y = "Absolute risk reduction") +
  theme_journal() +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))
print(fig2)

# ---- FIGURE 3: sample size by population and method --------------------------
lab_methode <- c(NSN_formule    = "Standard formula",
                 NSN_non_ajuste = "Simulation, unadjusted analysis",
                 NSN_ajuste     = "Simulation, adjusted analysis")

df_fig3 <- res1 %>%
  select(sous_groupe, NSN_formule, NSN_non_ajuste, NSN_ajuste) %>%
  pivot_longer(-sous_groupe, names_to = "methode", values_to = "NSN") %>%
  mutate(methode     = factor(lab_methode[methode], levels = lab_methode),
         sous_groupe = factor(sous_groupe, levels = grp_levels))

y_floor <- min(df_fig3$NSN, na.rm = TRUE) * 0.75

fig3 <- ggplot(df_fig3, aes(sous_groupe, NSN, color = methode)) +
  geom_linerange(aes(ymin = y_floor, ymax = NSN),
                 position = position_dodge(.65), linewidth = .5) +
  geom_point(position = position_dodge(.65), size = 3) +
  geom_text(aes(label = format(NSN, big.mark = ",", trim = TRUE)),
            position = position_dodge(.65), vjust = -1.1, size = 2.6,
            show.legend = FALSE) +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", trim = TRUE),
                expand = expansion(mult = c(0.02, 0.12))) +
  scale_color_manual(values = c("#6a7683", "#a30000", "#1c4e80"), name = NULL) +
  labs(x = NULL, y = "Patients per arm (log scale)") +
  theme_journal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  guides(color = guide_legend(nrow = 1))
print(fig3)

# ---- FIGURE 4: sample-size inflation factor by method ------------------------
fig4 <- res1 %>%
  transmute(sous_groupe = factor(sous_groupe, levels = grp_levels),
            `Unadjusted` = ratio_non_ajuste, `Adjusted` = ratio_ajuste) %>%
  pivot_longer(-sous_groupe, names_to = "analyse", values_to = "ratio") %>%
  mutate(analyse = factor(analyse, levels = c("Unadjusted", "Adjusted"))) %>%
  ggplot(aes(sous_groupe, ratio, fill = analyse)) +
  geom_col(position = position_dodge(.7), width = .6) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey40") +
  geom_text(aes(label = sprintf("×%.1f", ratio)),
            position = position_dodge(.7), vjust = -.4, size = 3) +
  scale_fill_manual(values = c("Unadjusted" = "#a30000",
                               "Adjusted" = "#1c4e80"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(x = NULL, y = "Required / formula sample size") +
  theme_journal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
print(fig4)

# ---- FIGURE 5: ratio according to the target odds ratio ----------------------
fig5 <- res2 %>%
  mutate(sous_groupe = factor(sous_groupe, levels = grp_levels)) %>%
  ggplot(aes(OR, ratio, color = sous_groupe)) +
  geom_hline(yintercept = 1, linetype = 3, color = "grey60") +
  geom_line(linewidth = .8) + geom_point(size = 2) +
  annotate("text", x = 0.90, y = 1.18, label = "Formula (reference)",
           color = "grey45", size = 3, hjust = 0) +
  scale_color_manual(values = pal, name = NULL) +
  scale_x_reverse(breaks = seq(0.50, 0.90, 0.05)) +
  scale_y_continuous(breaks = seq(1, 7, 1)) +
  coord_cartesian(ylim = c(0.9, max(res2$ratio) * 1.05)) +
  labs(x = "Patient-level odds ratio (stronger effect to the right)",
       y = "Required / formula sample size") +
  theme_journal() +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))
print(fig5)

# ---- FIGURE 6: sample-size inflation factor by effect scale ------------------
fig6 <- res3 %>%
  transmute(sous_groupe = factor(sous_groupe, levels = grp_levels),
            `Constant OR` = ratio_OR, `Constant ARR` = ratio_ARR,
            `Constant RR` = ratio_RR) %>%
  pivot_longer(-sous_groupe, names_to = "echelle", values_to = "ratio") %>%
  mutate(echelle = factor(echelle,
                          levels = c("Constant OR", "Constant ARR", "Constant RR"))) %>%
  ggplot(aes(sous_groupe, ratio, fill = echelle)) +
  geom_col(position = position_dodge(.75), width = .68) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey40") +
  geom_text(aes(label = sprintf("×%.1f", ratio)),
            position = position_dodge(.75), vjust = -.35, size = 2.8) +
  scale_fill_manual(values = c("Constant OR" = "#a30000",
                               "Constant ARR" = "#b8860b",
                               "Constant RR" = "#1c4e80"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(x = NULL, y = "Required / formula sample size") +
  theme_journal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
print(fig6)

# ---- FIGURE 7: sample size according to the target effect, whole population ---
set.seed(SEED)
df_fig7 <- nsn_master %>%
  filter(sous_groupe == "Overall") %>%
  transmute(OR, Formule = NSN_formule, Non_ajuste = NSN_non_ajuste)
df_fig7$Ajuste <- vapply(df_fig7$OR, function(or)
  find_nsn(function(N) power_adj(pool_gen, function(p) apply_OR(p, or), N,
                                 NSIM_ADJ, SIGMA_EPS, d_gen, "libre"),
           hi = 1e5), numeric(1))

lab7 <- c(Formule    = "Formula (homogeneous risk)",
          Non_ajuste = "Required, unadjusted analysis",
          Ajuste     = "Required, adjusted analysis")

fig7 <- ggplot(df_fig7, aes(x = OR)) +
  geom_ribbon(aes(ymin = Formule, ymax = Non_ajuste),
              fill = "#a30000", alpha = 0.08) +
  geom_line(aes(y = Formule, color = lab7["Formule"]), linewidth = 1) +
  geom_point(aes(y = Formule, color = lab7["Formule"]), size = 2) +
  geom_line(aes(y = Non_ajuste, color = lab7["Non_ajuste"]), linewidth = 1) +
  geom_point(aes(y = Non_ajuste, color = lab7["Non_ajuste"]), size = 2) +
  geom_line(aes(y = Ajuste, color = lab7["Ajuste"]), linewidth = 1,
            linetype = "longdash") +
  geom_point(aes(y = Ajuste, color = lab7["Ajuste"]), size = 2) +
  scale_y_log10(breaks = c(300, 1000, 3000, 10000, 30000, 100000),
                labels = function(x) format(x, big.mark = ",", trim = TRUE,
                                            scientific = FALSE)) +
  scale_x_reverse(breaks = seq(0.50, 0.90, 0.05)) +
  scale_color_manual(breaks = unname(lab7),
                     values = setNames(c("#6a7683", "#a30000", "#1c4e80"),
                                       unname(lab7)), name = NULL) +
  labs(x = "Patient-level odds ratio (stronger effect to the right)",
       y = "Patients per arm (log scale)") +
  theme_journal() + guides(color = guide_legend(nrow = 3))
print(fig7)

# ---- FIGURE 8: three effect scales, one anchor point -------------------------
p1_gen  <- pbar["Overall"]
g_ref   <- apply_OR(p1_gen, OR_REF)
RR_ref  <- g_ref / p1_gen
ARR_ref <- p1_gen - g_ref

set.seed(SEED)
fig8 <- tibble(p = seq(0, 1, length.out = 600)) %>%
  mutate(`Constant OR`  = apply_OR(p, OR_REF),
         `Constant RR`  = apply_RR(p, RR_ref),
         `Constant ARR` = apply_ARR(p, ARR_ref)) %>%
  pivot_longer(-p, names_to = "echelle", values_to = "traite") %>%
  mutate(echelle = factor(echelle,
                          levels = c("Constant OR", "Constant RR", "Constant ARR"))) %>%
  ggplot(aes(p, traite, color = echelle)) +
  geom_density(data = tibble(p = echant(pool_gen)),
               aes(x = p, y = after_stat(scaled) * 0.9),
               inherit.aes = FALSE, fill = "grey88", color = NA) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey55") +
  geom_line(linewidth = 1, alpha = 0.95) +
  annotate("point", x = p1_gen, y = g_ref, size = 3, color = "black") +
  annotate("text", x = p1_gen + 0.04, y = g_ref - 0.05,
           label = sprintf("anchor: mean risk = %.1f%%", 100 * p1_gen),
           hjust = 0, size = 3) +
  scale_color_manual(values = c("Constant OR" = "#a30000",
                                "Constant RR" = "#1c4e80",
                                "Constant ARR" = "#b8860b"), name = NULL) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Baseline risk of death (p)",
       y = "Risk of death under treatment") +
  theme_journal()
print(fig8)

# ---- FIGURE C1: ratio by center (forest plot) --------------------------------
ratio_poole <- unname(est_global["ratio"])

fig_C1 <- res_centres %>%
  mutate(centre = fct_reorder(centre, ratio)) %>%
  ggplot(aes(ratio, centre)) +
  geom_vline(xintercept = 1, linetype = 3, color = "grey55") +
  geom_vline(xintercept = ratio_poole, linetype = 2, color = "#a30000") +
  geom_linerange(aes(xmin = IC95_bas, xmax = IC95_haut), color = "grey35") +
  geom_point(aes(size = n), shape = 15, color = "#1c4e80") +
  scale_size_area(max_size = 4, guide = "none") +
  labs(x = "Unadjusted / formula sample size (bootstrap 95% CI)", y = NULL) +
  theme_journal()
print(fig_C1)

# ---- FIGURE C2: ratio according to the mean risk of the center ---------------
fig_C2 <- ggplot(res_centres, aes(risque_moyen, ratio)) +
  geom_hline(yintercept = 1, linetype = 3, color = "grey55") +
  geom_linerange(aes(ymin = IC95_bas, ymax = IC95_haut), color = "grey70") +
  geom_point(aes(size = n), color = "#1c4e80", alpha = .8) +
  annotate("point", x = unname(est_global["risque_moyen"]), y = ratio_poole,
           shape = 18, size = 5, color = "#a30000") +
  scale_size_area(max_size = 5, guide = "none") +
  labs(x = "Mean risk of death in the center",
       y = "Unadjusted / formula sample size") +
  theme_journal()
print(fig_C2)

# ---- FIGURE C3: shape of the distribution by center --------------------------
set.seed(SEED)
dens_c <- do.call(rbind, lapply(centres_elig, function(cc)
  data.frame(risk = clamp(echant(as.vector(risk_mat[df$centre == cc, , drop = FALSE]),
                                 5e4)), centre = cc)))

fig_C3 <- ggplot() +
  geom_density(data = dens_c, aes(risk, group = centre),
               color = "grey65", linewidth = .3) +
  geom_density(data = data.frame(risk = clamp(echant(pool_gen))), aes(risk),
               color = "#1c4e80", linewidth = 1) +
  ECHELLE_RISQUE +
  labs(x = "Predicted risk of death (log scale)", y = "Density") +
  theme_journal()
print(fig_C3)

# ---- FIGURE R2 (Fig. 2 of the manuscript): factor according to R² ------------
r2_max <- max(c(table_R2$R2_score, res_centres$R2_score)) + 0.05
courbes_R2 <- tibble(R2 = seq(0, r2_max, length.out = 300)) %>%
  mutate(ratio = 1 / (1 - R2)^2)

pop_order  <- c("Overall", "Hemorrhagic shock", "Penetrating", "Blunt",
                "Traumatic brain injury")
pop_labels <- c("Whole population", "Hemorrhagic shock", "Penetrating", "Blunt",
                "Traumatic brain injury")
pal_R2 <- pal[pop_order]

x_lab <- r2_max - 0.02
fig_R2 <- ggplot() +
  geom_hline(yintercept = 1, linetype = 3, color = "grey60") +
  geom_line(data = courbes_R2, aes(R2, ratio), color = "grey35", linewidth = .7) +
  annotate("text", x = x_lab, y = 1 / (1 - x_lab)^2, hjust = 1.1, vjust = -0.4,
           label = "frac(1, (1 - R^2)^2)", parse = TRUE, size = 3.5,
           color = "grey25") +
  geom_point(data = res_centres, aes(R2_score, ratio), color = "grey60", size = 2) +
  geom_point(data = table_R2,
             aes(R2_score, ratio_observe, color = sous_groupe,
                 shape = "Unadjusted analysis"), size = 3) +
  geom_point(data = table_R2,
             aes(R2_score, ratio_ajuste_observe, color = sous_groupe,
                 shape = "TRISS-adjusted analysis"), size = 3, stroke = 1) +
  scale_color_manual(values = pal_R2, breaks = pop_order, labels = pop_labels,
                     name = NULL) +
  scale_shape_manual(values = c("Unadjusted analysis" = 16,
                                "TRISS-adjusted analysis" = 5),
                     breaks = c("Unadjusted analysis", "TRISS-adjusted analysis"),
                     name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  coord_cartesian(clip = "off") +
  labs(x = expression(R^2 ~ "of baseline risk" ==
                        frac(Var(p), bar(p) ~ (1 - bar(p)))),
       y = "Required / formula sample size") +
  theme_journal() +
  theme(legend.box = "vertical", legend.box.just = "left",
        legend.spacing.y = unit(1, "mm")) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE, order = 1,
                              override.aes = list(shape = 16)),
         shape = guide_legend(order = 2))
print(fig_R2)

# ---- FIGURE PANEL (Fig. 1 of the manuscript): A–D ----------------------------
set.seed(SEED)
pA <- fig1 +
  labs(y = "Density") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_text()) +
  guides(fill = "none", color = "none")

pB <- fig2
pB$layers <- c(geom_density(data = tibble(p = echant(pool_gen)),
                            aes(x = p, y = after_stat(scaled) * 0.085),
                            bounds = c(0, 1), fill = "grey85", color = NA,
                            inherit.aes = FALSE),
               pB$layers)
pB <- pB +
  scale_x_continuous(labels = function(x) paste0(100 * x, "%")) +
  scale_y_continuous(labels = function(x) 100 * x) +
  labs(x = "Predicted risk of death", y = "Absolute risk reduction (points)") +
  guides(color = "none")

lab_panel <- c("Standard formula", "Required, unadjusted analysis",
               "Required, TRISS-adjusted analysis")
col_meth  <- setNames(c("#6a7683", "#a30000", "#1c4e80"), lab_panel)
pC <- df_fig3 %>%
  mutate(methode = factor(lab_panel[as.integer(methode)], levels = lab_panel),
         sous_groupe = fct_rev(sous_groupe)) %>%
  ggplot(aes(sous_groupe, NSN, color = methode, group = fct_rev(methode))) +
  geom_linerange(aes(ymin = y_floor, ymax = NSN), position = position_dodge(.7),
                 linewidth = .5, show.legend = FALSE) +
  geom_point(position = position_dodge(.7), size = 2.5) +
  geom_text(aes(label = format(NSN, big.mark = ",", trim = TRUE)),
            position = position_dodge(.7), hjust = -0.3, size = 2.4,
            show.legend = FALSE) +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", trim = TRUE),
                expand = expansion(mult = c(0.02, 0.20))) +
  scale_color_manual(values = col_meth, breaks = lab_panel, name = NULL) +
  coord_flip() +
  labs(x = NULL, y = "Patients per arm (log scale)") +
  theme_journal()

pD <- fig7 +
  geom_vline(xintercept = OR_REF, linetype = 3, color = "grey50") +
  scale_x_reverse(breaks = seq(0.50, 0.90, 0.10)) +
  labs(x = "Patient-level odds ratio") +
  guides(color = "none")

fig_panel <- ((pA | pB) / (pC | pD)) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom", legend.text = element_text(size = 7.5),
        legend.key.width = unit(3, "mm"),
        plot.tag = element_text(face = "bold", size = 12))
print(fig_panel)

# ==============================================================================
# STEP 18 — EXPORT OF AGGREGATED OUTPUTS (shareable folder)
# ==============================================================================

tables_export <- list(
  S0a_flux = flow_table, S0b_statut_variables = table_statut,
  S0c_manquants_centres = as.data.frame(manq_centres),
  S0d_centres = desc_centres_global,
  T1_caracteristiques = table1, T2_NSN = res1, T2bis_indice_R2 = table_R2,
  T3_echelles = res3,
  S1_validation = table_validation, S2_pools = table_pools,
  S3_dgp = table_dgp, S3bis_auc = table_auc_sg,
  S4_validation_analytique = table_validation_analytique,
  S5_typeI = table_typeI, S6_efficience = table_efficience,
  S7_grille_OR = as.data.frame(res2), S8_sens_auc = sens_auc,
  S9_sens_calib = sens_calib, S10_sigma_pop = sens_sigma_pop,
  S11_centres = res_centres, S12_synthese_centres = synthese_centres,
  S13_loco = loco, S14_fig1D = df_fig7, S15_puissance = table_puissance,
  S16_activite_centres = desc_activite)

for (nm in names(tables_export))
  write.csv2(tables_export[[nm]], file.path(OUTPUT_DIR, paste0(nm, ".csv")),
             row.names = FALSE, fileEncoding = "UTF-8")

figs_export <- list(
  fig_S1 = list(fig_S1, 5, 5), fig1 = list(fig1, 5, 8),
  fig2 = list(fig2, 6, 4.5), fig3 = list(fig3, 8, 5),
  fig4 = list(fig4, 7, 4.5), fig5 = list(fig5, 7, 5),
  fig6 = list(fig6, 8, 5), fig7 = list(fig7, 7, 5),
  fig8 = list(fig8, 6, 6), fig_C1 = list(fig_C1, 6, 7),
  fig_C2 = list(fig_C2, 6, 4.5), fig_C3 = list(fig_C3, 6, 4.5))

DEVICE_PDF <- if (Sys.info()[["sysname"]] == "Darwin")
  function(filename, ...) grDevices::quartz(file = filename, type = "pdf", ...) else "pdf"

for (nm in names(figs_export))
  ggsave(file.path(OUTPUT_DIR, paste0(nm, ".pdf")), figs_export[[nm]][[1]],
         width = figs_export[[nm]][[2]], height = figs_export[[nm]][[3]],
         device = DEVICE_PDF)

ggsave(file.path(OUTPUT_DIR, "fig_panel.pdf"), fig_panel,
       width = 174, height = 190, units = "mm", device = DEVICE_PDF)
ggsave(file.path(OUTPUT_DIR, "fig_R2.pdf"), fig_R2,
       width = 174, height = 130, units = "mm", device = DEVICE_PDF)

cat("\n\n=== OUTPUTS ===\n")
cat("Shareable folder:", normalizePath(OUTPUT_DIR), "\n")
cat("Confidential folder (DO NOT share):", normalizePath(LOCAL_DIR), "\n")
cat(sprintf("Mode = %s | M = %d | SIGMA_EPS = %.3f | target AUC = %.2f | O:E = %.2f\n",
            MODE, M, SIGMA_EPS, AUC_CIBLE, OE_RATIO))
sink()
