# =============================================================================
# Title: Figure 3: Side-by-side forest plot: Observational vs Mendelian Randomisation
#
# Goal:  Join site-level pooled estimates from observational meta-analysis and MR,
#        format display columns, and draw a two-panel forest plot with I² labels.
# =============================================================================

# ---- Packages ----------------------------------------------------------------
pacman::p_load("data.table", "tidyverse", "data.table", "forestploter", "devEMF", "readxl")

# ---- Inputs ------------------------------------------------------------------
# - MR meta-analysis results: "MR_metaanalysis.csv"
# - Observational meta-analysis results : "meta_analysis_res.csv"


# ---- MR meta-analysis results ------------------------------------------------
MR <- fread("MR_metaanalysis.csv") %>% 
  dplyr::select(Site, N_cases, OR, LCI, UCI, N_studies, I2) %>% 
  rename(N_cases_MR = N_cases, 
         N_studies_MR = N_studies,
         I2_MR = I2,
         OR_MR = OR, 
         LCI_MR = LCI,
         UCI_MR = UCI) %>%
  
  # Format text for plot
  # Harmonise site names to match observational results
  mutate(I2_MR = as.character(I2_MR),  # Convert I2_MR to character
         I2_MR = case_when(is.na(I2_MR) ~ "-", I2_MR == "0" ~ "0.00", TRUE ~ I2_MR),
         I2_MR = paste("   ", I2_MR), 
         N_cases_MR = format(N_cases_MR, big.mark = ",", scientific = FALSE),
         Site = case_when(Site == "Pancreatic" ~ "Pancreas", # Rename to match observational
                          Site == "Oesophageal (SQ)" ~ "Oesophageal (SQ, never-smokers)",
                          Site == "Lung (adj for smoking)" ~ "Lung (never-smokers)",
                          Site == "Lung" ~ "Lung (univariate MR)",
                          Site == "Head and neck (adj for smoking)" ~ "Head and neck (never-smokers)",
                          Site == "Head and neck" ~ "Head and neck (univariate MR)",
                          Site == "Bladder" ~ "Bladder (never-smokers)",
                          TRUE ~ Site )) %>% 
  
  # Drop univariate MR for lung and head & neck
  filter(Site!="Lung (univariate MR)" & 
           Site != "Head and neck (univariate MR)")

# ---- observational meta-analysis results ------------------------------------------------
# And combine with MR 

combined <- fread("meta_analysis_res.csv") %>% 
  dplyr::select(Site, N_cases, N_cohorts, OR, LCI, UCI, I2, IQR) %>%
  rename(`N cases` = N_cases,
         `N cohorts` = N_cohorts) %>%
  mutate(I2 = round(I2, 2)) %>%
  full_join(MR, by="Site") 

# ---- Display columns for forestploter ----------------------------------------
forest_dat <- combined %>%  
  # Two display columns of effect ± CI as strings (Obs vs MR)
  mutate(`Observational`=ifelse(is.na(OR), "",
                                sprintf("%.2f (%.2f, %.2f)   ",
                                        OR, 
                                        LCI, 
                                        UCI)),
         `Mendelian randomisation`= ifelse(is.na(OR_MR), "",
                      sprintf("%.2f (%.2f, %.2f)   ",
                              OR_MR, 
                              LCI_MR, 
                              UCI_MR)),
         # Special override for a site with missing MR
         `Mendelian randomisation`= ifelse(is.na(OR_MR) & (Site == "Gastric (cardia)"), 
                                           "NA",
                      `Mendelian randomisation`),
         I2_MR= ifelse(is.na(OR_MR) & (Site == "Gastric (cardia)"), 
                                           "   NA",
                       I2_MR),
         # Spacer columns to visually separate table groups
         ` `=paste(rep(" ", 50), collapse = " "),
         `  `=paste(rep(" ", 2), collapse = " ") ) %>%
  dplyr:: select(c(Site, `N cases`, `N cohorts`, N_cases_MR, `N_studies_MR`,
                   ` `, 
                   `Observational`, I2, 
                   `Mendelian randomisation`,  I2_MR,
                   OR, LCI, UCI, OR_MR, LCI_MR, UCI_MR,  
                   )) %>%
  
  # Compute SEs from 95% CIs for bubble sizes
  arrange(desc(OR)) %>% 
  mutate(se = (log(UCI) -log(LCI)) / 3.92,
         se_MR = (log(UCI_MR) -log(LCI_MR)) / 3.92,
         se_inv = (1/se)/120, 
         se_MRinv = (1/se_MR)/120,
         
         # Footnote marks (Unicode): dagger \u2020, asterisk *
         Site = recode(Site, "Bladder (never-smokers)" = "Bladder\u2020",
                       "Oesophageal (SQ, never-smokers)" = "Oesophageal (SQ)\u2020",
                       "Head and neck (never-smokers)"= "Head and neck*",
                       "Lung (never-smokers)" = "Lung*" )) %>%
  
  # Rename printed headers (newlines to stack headers nicely)
  rename("\nRR (95% CI)" = "Observational",
         "\nRR (95% CI) " = "Mendelian randomisation",
         "\n   I\U00B2"="I2", 
         "\nN studies"=`N cohorts`, 
         "\nN cases"=`N cases`,
         "MR\nN cases"=`N_cases_MR`, 
         "MR\nN studies"=`N_studies_MR`,
         "MR\n   I\U00B2"=`I2_MR`) 
  
# ---- Theme for forestploter ---------------------------------------------------

tm <- forest_theme(refline_lty = c("solid"),
                   refline_col = c("#bababa"),
                   base_size = 8,
                   legend_name = "Study type",
                   legend_value = c(" Observational", " Mendelian randomisation"),
                   ci_col = c("#0072B2", "#F5C710"),
                   core=list(bg_params=list(fill = c("white"))))

# ---- Output ------------------------------------------------------------------

pdf(file = paste0("Fig3.pdf"), width = 20, height = 10)

forest(forest_dat[,c(1:10)],
            est = list(forest_dat$OR,
                       forest_dat$OR_MR),
            lower = list(forest_dat$LCI,
                         forest_dat$LCI_MR),
            upper = list(forest_dat$UCI,
                         forest_dat$UCI_MR),
            sizes = list(forest_dat$se_inv,
                         forest_dat$se_MRinv),
            arrow_lab = c("Reduced risk", "Elevated risk"),
            ci_column = c(6),
            xlim = c(0.25, 2.0),
            ticks_at = c( 0.5, 1, 1.5, 2),
            ref_line = 1,
            nudge_y = 0.3,
            xlab = expression("RR per 5 kg/m"^2~"increase"),
            theme = tm)


dev.off()

