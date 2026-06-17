# =============================================================================
# ANNOTATION PIPELINE — Drug waste LC-MS (MS-DIAL output)
#
# Description:
#   This script performs suspect-level annotation of LC-HRMS features from
#   MS-DIAL alignment output. It matches observed m/z values against a
#   curated reference database of drug-related compounds, applies median
#   normalisation and blank correction, and visualises compound distributions
#   per sample group as heatmaps.
#
# Input:
#   Area_2_2026_04_22_08_38_08.csv — MS-DIAL peak area export
#
# Output:
#   annotaties.csv              — all suspect annotations (Level 2 and 3)
#   annotaties_ms2.mgf          — MS2 spectra for GNPS/MassBank validation
#   annotaties_met_areas.csv    — Level 2 annotations with normalised areas
#   groep_per_categorie.csv     — mean signal per drug category per group
#   groep_per_verbinding.csv    — mean signal per compound per group
#   aanwezig_per_groep.csv      — detection percentage per compound per group
#   heatmap_categorie.pdf       — heatmap by drug category
#   heatmap_compound.pdf        — heatmap by individual compound
#
# References:
#   Rainer et al. (2022), Metabolites 12(2):173
#     https://doi.org/10.3390/metabo12020173
#   Schymanski et al. (2014), Environ. Sci. Technol. 48(4):2097
#     https://doi.org/10.1021/es5002105
#   Ritchie et al. (2015), Nucleic Acids Res. 43:e47
#     https://doi.org/10.1093/nar/gkv007
#   Exact masses verified via PubChem (https://pubchem.ncbi.nlm.nih.gov)
# =============================================================================

# Installation (run once):
# BiocManager::install(c("MetaboAnnotation", "MetaboCoreUtils"))
# install.packages(c("tidyverse", "limma", "ggplot2", "viridis"))

library(MetaboAnnotation)   # Rainer et al. 2022
library(tidyverse)
library(limma)              # Ritchie et al. 2015
library(ggplot2)
library(viridis)

# Helper: convert European decimal notation (comma) to numeric
eu2num <- function(x) as.numeric(gsub(",", ".", as.character(x)))

# =============================================================================
# 1. Load MS-DIAL alignment data
# =============================================================================
# The MS-DIAL export contains 4 metadata rows above the actual header (skip=4).
# Columns include feature metadata and one peak area column per sample.

raw <- read.csv("Area_2_2026_04_22_08_38_08.csv",
                header = TRUE, skip = 4,
                check.names = FALSE, stringsAsFactors = FALSE)

# Query object: one row per alignment feature
query <- data.frame(
  feature_id    = as.integer(raw[["Alignment ID"]]),
  mz            = eu2num(raw[["Average Mz"]]),
  rt_min        = eu2num(raw[["Average Rt(min)"]]),
  adduct_msdial = raw[["Adduct type"]],
  fill_pct      = eu2num(raw[["Fill %"]]),
  msms_assigned = raw[["MS/MS assigned"]] == "True",
  msms_spectrum = raw[["MS/MS spectrum"]]
)

# Identify sample columns using strict regex (prefix + digit)
# This prevents metadata columns (e.g. "Metabolite name") from being included
blanco_cols <- grep("^Bl_", names(raw), value = TRUE)
real_cols   <- grep("^(BH|G[0-9]|H[0-9]|M[0-9]|T[0-9])", names(raw), value = TRUE)

# Helper: assign sample group label from column name
get_group <- function(nm) case_when(
  grepl("^BH", nm) ~ "BH",
  grepl("^G",  nm) ~ "G",
  grepl("^H",  nm) ~ "H",
  grepl("^M",  nm) ~ "M",
  grepl("^T",  nm) ~ "T",
  TRUE ~ NA_character_
)

cat("Features loaded:", nrow(query), "\n")
cat("Real sample columns:", length(real_cols), "\n")
cat("Blank columns:", length(blanco_cols), "\n")

# =============================================================================
# 2. Reference database
# =============================================================================
# Curated database of 41 drug-related compounds.
# Exact monoisotopic masses [M] (neutral molecule) verified via PubChem.
# PubChem Compound IDs (CIDs) are listed in comments for traceability.
# Categories: cocaine, adulterant, amphetamine, precursor-amphetamine,
#             MDMA precursor, 2CB, opioid, cathinone, solvent

target <- data.frame(
  compound_id = 1:41,
  name = c(
    "Cocaine",               "Benzoylecgonine",        # PubChem: 446220, 2723724
    "Ecgonine",              "Ecgonine methyl ester",  # 5758, 5770
    "Norcocaine",            "Cocaethylene",            # 65034, 107738
    "Levamisole",            "Phenacetin",              # 26879, 4754
    "Lidocaine",             "Caffeine",                # 3676, 2519
    "Paracetamol",           "Diltiazem",               # 1983, 39186
    "Amphetamine",           "Methamphetamine",         # 3007, 10836
    "MDMA",                  "MDA",          "MDEA",    # 1615, 1816, 26455
    "Ephedrine",             "Pseudoephedrine",         # 9294, 7028
    "Phenylacetone (P2P)",   "Phenylacetaldehyde",      # 7410, 998
    "BMK glycidate",                                    # 12673
    "PMK (3,4-MDP2P)",       "PMK-glycidate",           # 69591, 2723877
    "Safrole",               "Piperonal",               # 5152, 995
    "2C-B",                  "2C-I",         "2C-H",    # 7038, 10325, 16606
    "Fentanyl",              "Heroin",                  # 3676, 5462
    "Morphine",              "Tramadol",   "Methadone", # 5288826, 33741, 4095
    "Mephedrone",            "MDPV",        "Alpha-PVP",# 45266475, 20429, 11979
    "Acetone",               "Toluene",                 # 180, 1140
    "Benzaldehyde",          "Formic acid"              # 240, 284
  ),
  exactmass = c(
    303.14712, 289.13147, 185.10519, 199.12084, 289.13147, 317.16277,
    204.06688, 179.09463, 234.17325, 194.08038, 151.06333, 414.15564,
    135.10480, 149.12045, 193.11028, 179.09463, 207.12593,
    165.11535, 165.11535, 134.07317, 120.05752, 250.08904,
    192.07858, 280.09452, 162.06808, 150.03226,
    259.03950, 305.02510, 181.07898,
    336.20689, 369.15819, 285.10011, 263.18853, 309.20925,
    177.11535, 275.15212, 231.16232,
    58.04187,  92.06260, 106.04187,  46.00548
  ),
  categorie = c(
    rep("cocaine",              6),
    rep("adulterant",           6),
    rep("amphetamine",          5),
    rep("precursor-amphetamine",5),
    rep("MDMA precursor",       4),
    rep("2CB",                  3),
    rep("opioid",               5),
    rep("cathinone",            3),
    rep("solvent",              4)
  ),
  stringsAsFactors = FALSE
)

# =============================================================================
# 3. MS1 annotation via matchMz() (Rainer et al. 2022)
# =============================================================================
# Mass2MzParam converts exact masses to theoretical adduct m/z values and
# matches them against observed query m/z values within the given tolerances.
#
# Adducts considered: [M+H]+, [M+2H]2+, [M+NH4]+, [M+CH3CN+H]+
# These correspond to the adduct types assigned by MS-DIAL in this dataset.
#
# Tolerance settings (Schymanski et al. 2014):
#   tolerance = 0.005 Da  (absolute lower bound)
#   ppm = 10              (relative; gives Level 2 < 5 ppm, Level 3 5-10 ppm)

param <- Mass2MzParam(
  adducts   = c("[M+H]+", "[M+2H]2+", "[M+NH4]+", "[M+C2H3N+H]+"),
  tolerance = 0.005,
  ppm       = 10
)

mtch <- matchMz(query, target, param = param)
mtch  # print summary of Matched object

# =============================================================================
# 4. Extract results via matchedData() (Rainer et al. 2022)
# =============================================================================
# matchedData() returns one row per query-target match.
# Features without a match receive NA for all target columns.
# Annotation confidence assigned per Schymanski et al. (2014):
#   Level 2: putative annotation (mass error < 5 ppm)
#   Level 3: tentative candidate (mass error 5-10 ppm)

res <- matchedData(mtch, c("feature_id", "mz", "rt_min", "adduct_msdial",
                            "fill_pct", "msms_assigned",
                            "target_name", "target_exactmass",
                            "target_categorie", "adduct",
                            "score", "ppm_error")) |>
  as.data.frame() |>
  filter(!is.na(target_name)) |>
  mutate(
    schymanski = case_when(
      ppm_error < 5  ~ "Level 2 (< 5 ppm)",
      ppm_error < 10 ~ "Level 3 (5-10 ppm)"
    )
  ) |>
  arrange(ppm_error)

cat("\nMatches < 10 ppm:", nrow(res), "\n")
cat("Level 2 (< 5 ppm):", sum(res$schymanski == "Level 2 (< 5 ppm)"), "\n")
cat("Level 3 (5-10 ppm):", sum(res$schymanski == "Level 3 (5-10 ppm)"), "\n\n")
print(res)

write.csv(res, "annotaties.csv", row.names = FALSE)

# =============================================================================
# 5. MGF export for MS2 validation
# =============================================================================
# Level 1 identification (confirmed structure) requires MS2 spectral matching
# against a reference library (Schymanski et al. 2014).
# Export annotated features as MGF for upload to:
#   GNPS:     https://gnps.ucsd.edu
#   MassBank: https://massbank.eu

mgf_q <- filter(query,
                feature_id %in% unique(res$feature_id),
                msms_assigned,
                !is.na(msms_spectrum))

mgf <- unlist(lapply(seq_len(nrow(mgf_q)), function(i) {
  r     <- mgf_q[i, ]
  peaks <- sapply(strsplit(strsplit(trimws(r$msms_spectrum), " ")[[1]], ":"),
                  function(p) if (length(p) == 2)
                    paste(gsub(",", ".", p[1]), gsub(",", ".", p[2])) else NA)
  peaks <- peaks[!is.na(peaks)]
  if (!length(peaks)) return(NULL)
  c("BEGIN IONS",
    paste0("FEATURE_ID=", r$feature_id),
    paste0("PEPMASS=",     round(r$mz, 5)),
    paste0("RTINSECONDS=", round(r$rt_min * 60, 1)),
    peaks, "END IONS", "")
}))

writeLines(mgf, "annotaties_ms2.mgf")
cat("Saved: annotaties.csv | annotaties_ms2.mgf\n")

# =============================================================================
# 6. Median normalisation
# =============================================================================
# normalizeMedianValues() (limma; Ritchie et al. 2015) scales each sample
# column so that all samples share the same median, correcting for
# sample-to-sample loading differences.
# Applied to real sample columns only; blanks are kept separate.

area_raw <- raw[, c(real_cols, blanco_cols)]
area_raw[] <- lapply(area_raw, eu2num)
area_raw[is.na(area_raw)] <- 0
rownames(area_raw) <- as.integer(raw[["Alignment ID"]])

area_norm <- as.data.frame(
  normalizeMedianValues(as.matrix(area_raw[, real_cols]))
)
rownames(area_norm) <- rownames(area_raw)

# =============================================================================
# 7. Blank correction
# =============================================================================
# The mean signal of the 13 procedural blanks (Bl_miliQ) is subtracted
# from each normalised sample value on a per-feature basis (row-wise sweep).
# Values below zero are set to zero, indicating signal not exceeding background.

blanco_mean <- rowMeans(area_raw[, blanco_cols], na.rm = TRUE)

area_clean <- as.data.frame(
  sweep(as.matrix(area_norm), 1, blanco_mean, "-")
)
area_clean[area_clean < 0] <- 0
rownames(area_clean)  <- rownames(area_norm)
area_clean$feature_id <- as.integer(rownames(area_clean))

cat("Features with signal after blank correction:",
    sum(rowSums(area_clean[, real_cols]) > 0), "\n")
cat("Features fully zero after blank correction: ",
    sum(rowSums(area_clean[, real_cols]) == 0), "\n")

# =============================================================================
# 8. Link Level 2 annotations to normalised peak areas
# =============================================================================
# Only Level 2 annotations (< 5 ppm) are used here for highest confidence.
# Results are summarised per drug category and per compound per sample group.

res_l2 <- res |> filter(schymanski == "Level 2 (< 5 ppm)")

annotated <- res_l2 |>
  select(feature_id,
         compound   = target_name,
         category   = target_categorie,
         ppm_error,
         schymanski) |>
  inner_join(area_clean, by = "feature_id")

cat("Annotated features (Level 2) after blank correction:", nrow(annotated), "\n")

annotated_long <- annotated |>
  pivot_longer(cols      = all_of(real_cols),
               names_to  = "sample",
               values_to = "area_norm") |>
  mutate(group = get_group(sample)) |>
  filter(!is.na(group), area_norm > 0)

# Mean signal per drug category per sample group
group_category <- annotated_long |>
  group_by(category, group) |>
  summarise(mean_area  = mean(area_norm, na.rm = TRUE),
            n_features = n_distinct(feature_id),
            .groups    = "drop") |>
  arrange(category, group)

# Mean signal per compound per sample group
group_compound <- annotated_long |>
  group_by(compound, category, group) |>
  summarise(mean_area = mean(area_norm, na.rm = TRUE), .groups = "drop") |>
  arrange(category, compound, group)

cat("\n=== MEAN SIGNAL PER CATEGORY PER SAMPLE GROUP ===\n")
print(as.data.frame(group_category), row.names = FALSE)

write.csv(annotated,       "annotaties_met_areas.csv", row.names = FALSE)
write.csv(group_category,  "groep_per_categorie.csv",  row.names = FALSE)
write.csv(group_compound,  "groep_per_verbinding.csv", row.names = FALSE)

# =============================================================================
# 9. Heatmaps
# =============================================================================
# Colour scale: log1p-transformed viridis (direction = -1: darker = higher signal)
# Saved as PDF at high resolution for publication.

# --- Heatmap 1: mean signal per drug category per sample group ---------------
ggplot(group_category, aes(x = group, y = category, fill = mean_area)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(
    option    = "viridis",
    direction = -1,
    trans     = "log1p",
    breaks    = c(100, 1000, 10000, 100000),
    labels    = scales::label_scientific(),
    name      = "Mean peak area\n(log scale)"
  ) +
  labs(
    title = "Mean signal per drug category per sample group",
    x     = "Sample group",
    y     = "Drug category"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x       = element_text(angle = 0),
    panel.grid        = element_blank(),
    legend.position   = "right",
    legend.key.height = unit(1.5, "cm"),
    legend.key.width  = unit(0.5, "cm"),
    legend.title      = element_text(size = 10),
    legend.text       = element_text(size = 9),
    plot.margin       = margin(10, 20, 10, 10)
  )

ggsave("heatmap_categorie.pdf", width = 9, height = 5)

# --- Heatmap 2: mean signal per compound per sample group --------------------
ggplot(group_compound, aes(x = group, y = compound, fill = mean_area)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(
    option    = "viridis",
    direction = -1,
    trans     = "log1p",
    breaks    = c(100, 1000, 10000, 100000),
    labels    = scales::label_scientific(),
    name      = "Mean peak area\n(log scale)"
  ) +
  labs(
    title = "Drug signal fingerprint per sample group",
    x     = "Sample group",
    y     = "Compound"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid        = element_blank(),
    legend.position   = "right",
    legend.key.height = unit(1.5, "cm"),
    legend.key.width  = unit(0.5, "cm"),
    legend.title      = element_text(size = 10),
    legend.text       = element_text(size = 9),
    plot.margin       = margin(10, 20, 10, 10)
  )

ggsave("heatmap_compound.pdf", width = 9, height = 8)

# =============================================================================
# 10. Presence table: detection per compound per sample group
# =============================================================================
# Includes all annotation levels (Level 2 and 3).
# Reports mean area, number of detections, total samples, and detection %.

presence <- res |>
  select(feature_id,
         compound  = target_name,
         category  = target_categorie,
         adduct,
         ppm_error,
         schymanski) |>
  inner_join(area_clean, by = "feature_id") |>
  pivot_longer(cols      = all_of(real_cols),
               names_to  = "sample",
               values_to = "area") |>
  mutate(group = get_group(sample)) |>
  filter(!is.na(group))

presence_table <- presence |>
  group_by(compound, category, schymanski, group) |>
  summarise(
    mean_area    = mean(area, na.rm = TRUE),
    n_detections = sum(area > 0),
    n_samples    = n(),
    pct_detected = round(sum(area > 0) / n() * 100),
    .groups      = "drop"
  ) |>
  filter(mean_area > 0) |>
  arrange(category, compound, group)

cat("\n=== DETECTED COMPOUNDS PER SAMPLE GROUP ===\n")
print(as.data.frame(presence_table), row.names = FALSE)
write.csv(presence_table, "aanwezig_per_groep.csv", row.names = FALSE)

cat("\nSaved:\n")
cat("  annotaties.csv\n")
cat("  annotaties_ms2.mgf\n")
cat("  annotaties_met_areas.csv\n")
cat("  groep_per_categorie.csv\n")
cat("  groep_per_verbinding.csv\n")
cat("  aanwezig_per_groep.csv\n")
cat("  heatmap_categorie.pdf\n")
cat("  heatmap_compound.pdf\n")
