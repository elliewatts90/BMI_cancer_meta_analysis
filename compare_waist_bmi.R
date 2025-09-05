# =============================================================================
# Title: BMI vs Waist — subgroup meta-analysis and Wald test of difference
# 
#   - For each cancer site, fit random-effects meta-analyses by Exposure (BMI vs Waist)
#   - Extract subgroup pooled RRs and 95% CIs
#   - Test BMI vs Waist difference using a correlated Wald test (r = 0.87)
#
# =============================================================================


# ---- Packages ----------------------------------------------------------------
pacman::p_load("tidyverse", "data.table", "janitor", "meta", "readxl", "scales")


# ---- Helpers -----------------------------------------------------------------
# Pretty p-value formatter
pvalr <- function(pvals, sig.limit = .001, digits = 3, html = FALSE) {
  
  roundr <- function(x, digits = 1) {
    res <- sprintf(paste0('%.', digits, 'f'), x)
    zzz <- paste0('0.', paste(rep('0', digits), collapse = ''))
    res[res == paste0('-', zzz)] <- zzz
    res
  }
  
  sapply(pvals, function(x, sig.limit) {
    if (x < sig.limit)
      if (html)
        return(sprintf('&lt; %s', format(sig.limit))) else
          return(sprintf('< %s', format(sig.limit)))
    if (x > .1)
      return(roundr(x, digits = 2)) else
        return(roundr(x, digits = digits))
  }, sig.limit = sig.limit)
}

# Round to 2 decimals unless value < 0.01; then keep to first non-zero place
my_round = function(x, n=3) {
  max(abs(round(x, n)), abs(signif(x, 1))) * sign(x)
} 


# ---- Data --------------------------------------------------------------------
cancer_allsites <- read_excel("SupplementaryData.xlsx", 
                              sheet = "BMI_waist_comparison") 

# Move first row to column names
colnames(cancer_allsites) <- as.character(unlist(cancer_allsites[1, ]))
cancer_allsites <- cancer_allsites[-1, ] # Remove the old first row

# Format data
cancer_allsites <- cancer_allsites %>%
  mutate(OR = as.numeric(Risk))%>%
  mutate(Lci = as.numeric(LCI)) %>%
  mutate(Uci = as.numeric(UCI)) %>% 
  filter(Forest_cancer!= "Prostate (overall)" )

# Sites to analyze
cancer_list <- unique(cancer_allsites$Forest_cancer)
cancer_list <- cancer_list[!is.na(cancer_list)]

### Investigate het by bmi/waist

# ---- BMI–Waist correlation (from PMID: 34654381) -----------------------------
# Used for the correlated Wald test comparing subgroup pooled effects.

r <- 0.87

# ---- Per-site loop: subgroup meta + Wald test --------------------------------

het_res <- data.frame()

for (c in cancer_list) {
  
  cancer_dta <- cancer_allsites %>%
    subset(Forest_cancer == c) 

  # Random-effects meta-analysis with subgroup = Exposure (BMI vs Waist)
  metares <- metagen(TE = log(OR),  
                     lower=log(Lci),  
                     upper = log(Uci),  
                     level.ci = 0.95,  
                     data = cancer_dta,  
                     subgroup = Exposure,
                     method.tau = "PM") 
  
  # Back-transform subgroup pooled effects to RR scale
  Risk_sum <- exp(metares[["TE.random.w"]])
  LCI_sum <- exp(metares[["lower.random.w"]]) 
  UCI_sum <-  exp(metares[["upper.random.w"]]) 
  
  se_vec <- metares[["seTE.random.w"]]
  
  # By-subgroup p-values (random-effects) — as vector
  PValue_sum <- pvalr(metares[["pval.random.w"]])
  
  # Subgroup labels
  Measurement <- metares[["bylevs"]]
  
  # Case counts by Exposure
  N_cases_bmi <- cancer_dta %>%
    filter(Exposure == "BMI") %>%
    mutate(`N cases` = suppressWarnings(as.numeric(`N cases`))) %>%
    summarise(N_cases_bmi = sum(`N cases`, na.rm = TRUE), .groups = "drop") %>%
    dplyr::pull(N_cases_bmi)

  N_cases_waist <- cancer_dta %>%
    filter(Exposure == "Waist") %>%
    mutate(`N cases` = suppressWarnings(as.numeric(`N cases`))) %>%
    summarise(N_cases_waist = sum(`N cases`, na.rm = TRUE), .groups = "drop") %>%
    dplyr::pull(N_cases_waist)
  
  # ---- Wald test: BMI vs Waist (correlated) ----------------------------------
  idx_BMI   <- match("BMI",   Measurement)
  idx_Waist <- match("Waist", Measurement)
  

  # Variances from 95% CIs : (CI width / (2*1.96))^2
  var_BMI <- se_vec[idx_BMI]^2
  var_waist <- se_vec[idx_Waist]^2
  
  # Wald statistic for difference of two correlated estimates:
  # z = (log(RR_BMI) - log(RR_Waist)) / sqrt(var_BMI + var_Waist - 2*r*sqrt(var_BMI*var_Waist))
  
  numerator <- log(Risk_sum[idx_BMI])- log(Risk_sum[idx_Waist])
  
  denominator <- var_BMI + var_waist - (2 * r * sqrt(var_BMI * var_waist))
 
  # Wald statistic 
  wald <- numerator / sqrt(denominator)
  
  Het_grps <- pvalr(2 * pnorm(-abs(wald)))
  Het_grps_nf <- 2 * pnorm(-abs(wald))
  
  temp <- data.frame(Site = paste(c), 
                     Measurement = Measurement, 
                     N_cases_waist = N_cases_waist,  
                     N_cases_bmi = N_cases_bmi,
                     OR = Risk_sum, 
                     LCI = LCI_sum, 
                     UCI = UCI_sum, 
                     `RR (95% CI)` = sprintf("%.2f (%.2f, %.2f)",
                                             as.numeric(Risk_sum),
                                             as.numeric(LCI_sum),
                                             as.numeric(UCI_sum)),
                     `P-value` = PValue_sum, 
                     Phet = Het_grps, 
                     Phet_nf = Het_grps_nf)
  
  het_res <- plyr::rbind.fill(temp, het_res)
  
}

fwrite(het_res, "wc_bmi_compare_HRs_wald.csv")
