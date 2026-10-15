#!/usr/bin/env Rscript
# Main analysis script. Run from the repository root.
# All models adjust for age, sex and years of education.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop('Usage: Rscript run_analysis.R INPUT.csv NEW_OUTPUT_DIR')
}
required <- c('lmerTest', 'emmeans', 'ggplot2')
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop('Install required packages: ', paste(missing, collapse = ', '))
suppressPackageStartupMessages({library(lmerTest); library(emmeans)})
options(stringsAsFactors = FALSE)
input <- normalizePath(args[1], mustWork = TRUE)
out <- args[2]
if (file.exists(out)) stop('Output already exists. Choose a new directory.')

# Check the private analysis table and calculate the participant-level ROI.
source('R/functions.R')
d <- prepare_data(input)
p <- subset(d, Electrode != 'ROI')
p$Electrode <- factor(p$Electrode, levels = sites)
roi <- subset(d, Electrode == 'ROI')
u <- subset(d, Electrode == 'FCz')
dir.create(out, recursive = TRUE)

# Group differences, ERP-phenotype associations and medication models.
source('R/models.R')
tables <- list(Clinical = cl, ERP_groups = gr, Group_electrode = gi,
               ERP_phenotype = assoc, ERP_group_slopes = si, Medication = med)
for (name in names(tables)) {
  write.csv(tables[[name]], file.path(out, paste0(name, '.csv')), row.names = FALSE, na = '')
}

# Table 2 uses the FCz-specific FDR family rather than the five-electrode family.
fcz <- subset(assoc, Site == 'FCz' & Tier == 'Primary')
fcz$q <- fcz$q_within_site
fcz$family_n <- fcz$local_family_n
write.csv(fcz, file.path(out, 'Table2_FCz_associations.csv'), row.names = FALSE)

# Pairwise differences between electrode-specific group contrasts.
post <- data.frame()
for (m in gi$ERP[gi$q < .05]) {
  model <- models[[paste('LMM', 'Main', m, sep = '__')]]
  marginal <- emmeans(model, ~ GROUP | Electrode, lmer.df = 'satterthwaite', weights = 'equal')
  contrasts <- contrast(marginal, list(FEP_minus_HC = c(-1, 1)))
  result <- as.data.frame(summary(pairs(contrasts, by = NULL, adjust = 'none'), infer = c(TRUE, TRUE)))
  result$p_Holm <- p.adjust(result$p.value, 'holm')
  result$ERP <- m
  post <- rbind(post, result)
}
write.csv(post, file.path(out, 'Posthoc_electrode_contrasts.csv'), row.names = FALSE)

# Figure 2 and software versions. No participant-level model objects are exported.
source('R/figure2.R')
capture.output(sessionInfo(), file = file.path(out, 'sessionInfo.txt'))
cat('\nAnalysis complete: ', normalizePath(out), '\n')
