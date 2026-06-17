# =============================================================================
# MS-DIAL NORMALISATION AND MDS ANALYSIS
# (all samples + blanks included)
#
# References:
#   Ritchie et al. (2015), Nucleic Acids Res. 43:e47          (limma)
#   Robinson et al. (2010), Bioinformatics 26:139             (edgeR)
#   Su et al. (2017), Bioinformatics 33:2050                  (Glimma)
# =============================================================================

# Installation (run once):
# BiocManager::install(c("edgeR", "Glimma", "limma"))
# install.packages("dplyr")

library(edgeR)
library(Glimma)
library(limma)
library(dplyr)

# =============================================================================
# 1. Prepare count table from MS-DIAL alignment output
# =============================================================================
# Input: MS-DIAL area CSV export (e.g. Area_2_2026_04_22_08_38_08.csv)

projectcount <- read.csv(
  "Area_2_2026_04_22_08_38_08.csv",
  skip        = 4,
  sep         = ",",
  check.names = FALSE
)

# Create feature IDs from m/z and retention time
projectcount$FeatureID <- paste0(
  "MZ_",  round(as.numeric(gsub(",", ".", projectcount[["Average Mz"]])),       4),
  "_RT_", round(as.numeric(gsub(",", ".", projectcount[["Average Rt(min)"]])), 2)
)

# Identify sample columns (everything after MS/MS spectrum, excluding last 2)
start_col <- which(colnames(projectcount) == "MS/MS spectrum") + 1
end_col   <- ncol(projectcount) - 2

sample_data <- projectcount[, start_col:end_col]

# Build and export count table
count_table            <- cbind(FeatureID = projectcount$FeatureID, sample_data)
rownames(count_table)  <- count_table$FeatureID
count_table$FeatureID  <- NULL

write.csv(count_table, "MSDIAL_count_table.csv")

# =============================================================================
# 2. Prepare factor data from count table column names
# =============================================================================
# Factor data is derived from the sample column names in the count table.
# Expected naming convention: Place_Layer_Number (e.g. BH_boven_1, G_NB_3)
# Blanks are excluded from the factor table (rows 12–24 in the original).

low           <- read.csv("MSDIAL_count_table.csv", header = TRUE)
row.names(low) <- low[, 1]

# Sample name vector (excluding first column = feature ID)
all_sample_names <- colnames(low)[-1]

# Parse Place and Layer from sample names
# Adjust indices below if column order in your export differs
# BH samples: columns 2-12 (indices 1-11), blanks: 13-25, rest: 26-61
Fac_all <- data.frame(
  Sample = all_sample_names,
  Place  = sub("_.*", "", all_sample_names),          # first part before underscore
  Layer  = sub("^[^_]+_([^_]+)_.*", "\\1", all_sample_names)  # second part
)

# Save factor file for reuse
write.table(Fac_all, "Factordata_MSDIAL.txt", row.names = FALSE)

# Remove blanks from factor table (keep only real samples)
Fac <- Fac_all[!grepl("^blanko|^Bl_", Fac_all$Sample, ignore.case = TRUE), ]
write.table(Fac, "Fac_MSDAIL_rest.txt", row.names = FALSE)

# Layer subset (BH + selected G and H samples with depth information)
Fac_layers <- Fac[grepl("boven|midden|onder", Fac$Layer, ignore.case = TRUE), ]
write.table(Fac_layers, "Fac_MSDAIL_layers.txt", sep = ",",
            col.names = TRUE, row.names = FALSE)

# =============================================================================
# 3. Load and transform data
# =============================================================================
# All samples + blanks, log2-transformed
Rest <- log2(low[, -1] + 1)

# Blank columns (columns 13–25 in original data)
blank_cols  <- grep("^blanko|^Bl_", colnames(Rest), ignore.case = TRUE)
row_means     <- rowMeans(low[, blank_cols + 1], na.rm = TRUE)  # +1 for ID col
row_means_log <- log2(row_means + 1)

# =============================================================================
# 4. Blank correction
# =============================================================================
# Log2 blank mean subtracted per feature; negative values set to zero
Rest_corrected              <- Rest - row_means_log
Rest_corrected[Rest_corrected < 0] <- 0
write.table(Rest_corrected, "Data_MSDAIL_corrected.txt")

# =============================================================================
# 5. Median normalisation
# =============================================================================
# Not blank-corrected
sample_idx   <- which(colnames(Rest) %in% Fac$Sample)
Rest_samples <- Rest[, sample_idx]
write.table(Rest_samples, "Data_MSDAIL.txt")

dat.trans2           <- as.data.frame(normalizeMedianValues(as.matrix(Rest_samples)))
row.names(dat.trans2) <- row.names(Rest)
colnames(dat.trans2)  <- Fac$Sample

# Blank-corrected
Rest_corr_samples   <- Rest_corrected[, sample_idx]
dat.trans2_corrected <- as.data.frame(
  normalizeMedianValues(as.matrix(Rest_corr_samples))
)
row.names(dat.trans2_corrected) <- row.names(Rest)
colnames(dat.trans2_corrected)  <- Fac$Sample

# =============================================================================
# 6. Boxplots: normalisation check
# =============================================================================
boxplot(Rest_samples,        main = "Raw data (not blank corrected)",           ylab = "log2 intensity", las = 2, cex.axis = 0.6)
boxplot(dat.trans2,          main = "Normalised data (not blank corrected)",     ylab = "log2 intensity", las = 2, cex.axis = 0.6)
boxplot(Rest_corr_samples,   main = "Raw data (blank corrected)",                ylab = "log2 intensity", las = 2, cex.axis = 0.6)
boxplot(dat.trans2_corrected,main = "Normalised data (blank corrected)",         ylab = "log2 intensity", las = 2, cex.axis = 0.6)

# =============================================================================
# 7. Experimental design (all samples)
# =============================================================================
experimentFactors <- lapply(apply(Fac, 2, split, ""), unlist)
experimentFactors <- as.data.frame(lapply(experimentFactors, as.factor))

Groups <- as.factor(paste0(experimentFactors$Place, "_", experimentFactors$Layer))
design <- model.matrix(~0 + Groups)

fixCols         <- paste(c("Groups", "experimentFactors", "\\$", "\\:", "\\-",
                            colnames(experimentFactors)), sep = "", collapse = "|")
colnames(design) <- gsub(fixCols, "", colnames(design))

samples       <- as.character(Groups)
names(samples) <- Fac$Sample

# =============================================================================
# 8. MDS: all samples (raw and normalised)
# =============================================================================
# Raw
countDF               <- Rest_samples[, names(samples)]
countDF[is.na(countDF)] <- 0
f <- DGEList(counts = countDF, group = as.character(samples))

glMDSPlot(f,
          labels = rownames(f$samples),
          groups = Fac,
          folder = "Glimma_plots_MSDAIL",
          launch = TRUE)

# Normalised
countDF.trans2               <- dat.trans2[, names(samples)]
countDF.trans2[is.na(countDF.trans2)] <- 0
f.trans2 <- DGEList(counts = countDF.trans2, group = as.character(samples))

glMDSPlot(f.trans2,
          labels = rownames(f.trans2$samples),
          groups = Fac,
          folder = "Glimma_plots_MSDAIL_trans",
          launch = TRUE)

# =============================================================================
# 9. Layer subset analysis
# =============================================================================
Fac_layers2 <- read.table("Fac_MSDAIL_layers.txt",
                           header = TRUE, sep = ",")

experimentFactors_layers <- lapply(apply(Fac_layers2, 2, split, ""), unlist)
experimentFactors_layers <- as.data.frame(lapply(experimentFactors_layers, as.factor))

Groups_layers  <- as.factor(paste0(experimentFactors_layers$Place, "_",
                                   experimentFactors_layers$Layer))
samples_layers <- as.character(Groups_layers)
names(samples_layers) <- Fac_layers2$Sample

# Layer subset: raw
countDF_layers               <- Rest_samples[, names(samples_layers)]
countDF_layers[is.na(countDF_layers)] <- 0
countDF_layers_clean <- countDF_layers %>%
  filter(!if_any(everything(), ~ is.infinite(.x) & .x < 0))

f_layers <- DGEList(counts = countDF_layers_clean,
                    group  = as.character(samples_layers))

glMDSPlot(f_layers,
          labels = rownames(f_layers$samples),
          groups = Fac_layers2,
          folder = "Glimma_plots_MSDAIL_layers",
          launch = TRUE)

# Layer subset: normalised
dat.trans2_layers <- as.data.frame(
  normalizeMedianValues(as.matrix(Rest_samples[, names(samples_layers)]))
)
colnames(dat.trans2_layers) <- Fac_layers2$Sample

countDF.trans2_layers               <- dat.trans2_layers[, names(samples_layers)]
countDF.trans2_layers[is.na(countDF.trans2_layers)] <- 0
countDF.trans2_layers_clean <- countDF.trans2_layers %>%
  filter(!if_any(everything(), ~ is.infinite(.x) & .x < 0))

f.trans2_layers <- DGEList(counts = countDF.trans2_layers_clean,
                            group  = as.character(samples_layers))

glMDSPlot(f.trans2_layers,
          labels = rownames(f.trans2_layers$samples),
          groups = Fac_layers2,
          folder = "Glimma_plots_MSDAIL_trans_layers",
          launch = TRUE)

# =============================================================================
# 10. MDS subset: up to 3 random samples per location group
# =============================================================================
set.seed(42)
groups_unique    <- unique(Fac$Place)
selected_samples <- c()

for (grp in groups_unique) {
  samples_in_group <- Fac$Sample[Fac$Place == grp]
  n <- min(3, length(samples_in_group))
  selected_samples <- c(selected_samples, sample(samples_in_group, n))
}

subset_counts               <- dat.trans2[, selected_samples]
subset_counts[is.na(subset_counts)] <- 0

subset_Fac    <- Fac[Fac$Sample %in% selected_samples, ]
subset_groups <- as.factor(paste0(subset_Fac$Place, "_", subset_Fac$Layer))

f_subset <- DGEList(counts = subset_counts, group = subset_groups)

glMDSPlot(f_subset,
          labels = subset_Fac$Sample,
          groups = subset_Fac,
          folder = "Glimma_MDS_subset_3perGroup",
          launch = TRUE)
