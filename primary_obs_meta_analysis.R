# =============================================================================
# Meta-analysis of Observational Studies: BMI & Cancer
# =============================================================================
# PURPOSE
#   Run the primary random-effects meta-analyses by cancer site using study-level
#   risk ratios (RR) and 95% CIs. Perform sensitivity analyses (Egger’s test and
#   trim-and-fill, plus planned exclusions) and test heterogeneity by geographic
#   region (RegionG) and output forest plots for each cancer by study
#
# INPUT
#   File:  "SupplementaryData.xlsx", sheet = "BMI_cancer_data"
#
# METHODS
#   - Model: meta::metagen on log(RR); random effects only; tau^2 via Paule–Mandel.
#   - Report: pooled RR (95% CI), p-value, I²; subgroup heterogeneity p_het by RegionG.
#   - Small-study effects: Egger’s test (meta::metabias(method.bias = "linreg")).
#   - Publication bias adjustment: trim-and-fill (meta::trimfill) with funnel plots.
#
# OUTPUT
#   - pooled estimates by cancer site.
#   - regional subgroup estimates + heterogeneity tests.
#   - sensitivity_results: Egger & trim-and-fill summary per cancer.
#   - Figures: funnel plots with Egger & trim-and-fill results (Supplementary Figure 27) 
#              forest plots showing individual study results for each cancer type (Supplementary Figures 2-26)
#
# =============================================================================

## Load packages 
pacman::p_load("tidyverse", "data.table", "janitor", "meta", "readxl", "dmetar")

# Functions for p-value rounding

# round to 2 decimal places unless value is < 0.01, else round to the last non zero decimal 
my_round = function(x, n=3) {
  max(abs(round(x, n)), abs(signif(x, 1))) * sign(x)
} 

# Alternative p-value function
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

# ---- Data formatting ----
cancer_allsites <- read_excel("SupplementaryData.xlsx", 
                              sheet = "BMI_cancer_data") 

# Move first row to column names
colnames(cancer_allsites) <- as.character(unlist(cancer_allsites[1, ]))
cancer_allsites <- cancer_allsites[-1, ] # Remove the old first row

# Variable formatting
cancer_allsites <- cancer_allsites %>%
   mutate(RR = as.numeric(`Risk estimate`),
          LCI = as.numeric(`95% lower confidence interval`),
          UCI = as.numeric(`95% upper confidence interval`),
          `N analytic set` = as.numeric(`N analytic set`),
          `N cases` = as.numeric(`N cases`),
          `N cases` = round(`N cases`, 0),
          `N cohorts` = as.numeric(`N cohorts`),
          RegionG = factor(Region , levels=c("Australia", 
                                             "East Asia",
                                             "Europe", 
                                             "North America", 
                                             "Middle East", 
                                             "South Asia"), 
          labels=c("Australia", "East Asia", "Europe", "North America", "Middle East", "South Asia"))) %>%
  rename(smok_indicator = "Smoking adjustment",
         health_dat = "Large health database")

# Create a list of cancers for loop
cancer_list <- unique(cancer_allsites$`Top 25 cancer`)


# ---- Meta-analysis per site ----

results <- data.frame()

for (c in cancer_list){

cancer_dta <- cancer_allsites %>%
  subset(`Top 25 cancer` == c) 

# Primary random-effects meta-analysis on log RR
metares <- metagen(TE = log(RR),  
                   lower=log(LCI),  
                   upper = log(UCI),  
                   level.ci = 0.95,  
                   data = cancer_dta,  
                   sm = "RR",
                   random = TRUE, 
                   common = FALSE, 
                   method.tau = "PM"  # Paule-Mandel
                   ) 

# Back-transform pooled estimate & CI
Risk_sum <- round(exp(metares[["TE.random"]]), 2) 
LCI_sum <- round(exp(metares[["lower.random"]]), 2) 
UCI_sum <-  round(exp(metares[["upper.random"]]), 2) 

# Study counts
N_cohort <- sum(cancer_dta$`N cohorts`)

# Case counts (special adjustment for head & neck cancers)
N_cases <- sum(cancer_dta$`N cases`, na.rm = TRUE)

if (c == "Head and neck (never-smokers)") {
  N_cases <- N_cases + 796

  }

N_cases_formatted <- format(N_cases, big.mark = ",", scientific = FALSE)

# IQR of individual-study RRs (Q1, median, Q3)
iqr_hr <- round(quantile(cancer_dta$RR, 
                         probs=c(0.25, 0.5, 0.75)), 2) 

  # IQR of individual-study RRs (Q1, median, Q3)
temp <- data.frame(Site = paste(c), 
                   N_cases = N_cases_formatted, 
                   N_cohorts = N_cohort, 
                   OR = Risk_sum, 
                   LCI = LCI_sum, 
                   UCI = UCI_sum, 
                   `RR (95% CI)` = paste0(formatC(Risk_sum, digits = 2, format = "f")," (",
                                          formatC(LCI_sum, digits = 2, format = "f"),", ", 
                                          formatC(UCI_sum, digits = 2, format = "f"),")"), 
                   `P-value` = metares[["pval.random"]], 
                   I2 = metares[["I2"]],
                   IQR = paste0(formatC(iqr_hr[[1]], digits = 2, format = "f"),", ", 
                                formatC(iqr_hr[[3]], digits = 2, format = "f")))

# merge together
results <- rbind(temp,results)

}

# Save-out
fwrite(results, "meta_analysis_res.csv")


# --------- Run sensitivity analyses ----------

# ---- Egger and trim-and-fill  ----

## Sensitivity analyses - Egger test, trim-and-fill, and Funnel plots (Supplementary Figure 27) ##

for (c in cancer_list){
  
  cancer_dta <- cancer_allsites %>%
    subset(`Top 25 cancer` == c) 

  # Cancer type to sentence case for axis labels
  cancer_string <- paste0(tolower(substr(c, 1, 1)), substr(c, 2, nchar(c)))
  
  # Primary meta-analysis (random effects)
  metares <- metagen(TE = log(RR),  
                     lower=log(LCI),  
                     upper = log(UCI),  
                     level.ci = 0.95,  
                     data = cancer_dta,  
                     sm="RR",
                     random = TRUE, 
                     common = FALSE, 
                     method.tau = "PM") 
  
  # Egger's test
  res.et <- eggers.test(metares) 
  selected_values <- res.et[c("intercept", "llci", "ulci", "t", "p")]
  
  # Trim-and-fill
  tf_result <- trimfill(metares)
  
  summary(tf_result)
  
  # Store statistics for plot annotations
  n_filled_studies <- sum(tf_result$trimfill)
  
  # Store fill-and-trim summary estimate
  random_or <- exp(tf_result$TE.random)
  random_CI <- paste0(round(random_or,2)," (", round(exp(tf_result$lower.random),2),
                      " to ", 
                      round(exp(tf_result$upper.random),2),
                      ")")
  
  # Primary random-effects OR and CI
  original_or <- exp(metares$TE.random)
  original_CI <- paste0(round(original_or,2)," (", round(exp(metares$lower.random),2),
                        " to ", 
                        round(exp(metares$upper.random),2),
                        ")")
  
  # Format Egger results for plot annotations
  egger_intercept <- round(selected_values$intercept, 2)
  egger_ci <- paste0(round(selected_values$llci, 2), " to ", 
                     round(selected_values$ulci, 2))
  egger_p <- selected_values$p
  
  # Difference in RR from trim-and-fill to RR from primary meta-analysis  
  change_in_or <- (exp(tf_result$TE.random) - exp(metares$TE.random) )
  
  # --- Plot: funnel with contours & filled studies ---
  png(filename = paste0( c, "_fill_trim_plot.png"), width = 8, height = 8, units = "in", res = 300)
  
  par(mar = c(10, 4, 4, 2))  # Increase the bottom margin to 10 for text annotations
  
  # Create funnel plot with filled studies
  funnel(tf_result,
         xlab = paste("Risk ratio", cancer_string),
         contour = c(0.9, 0.95, 0.99),
         xlim = c(0.2, 10), 
         pch=ifelse(tf_result$trimfill, 1, 16),
         col = ifelse(tf_result$trimfill, "red", "black"),
         col.contour = c("darkgray", "gray", "lightgray"))
  
  # Add a legend for both contours and points
  legend(x = 3, y = 0.01,
         legend = c("p < 0.1", "p < 0.05", "p < 0.01", 
                    "Observed Studies", "Filled Studies"),
         col = c("darkgray", "gray", "lightgray", 
                 "black", "red"),
         pch = c(NA, NA, NA, 16, 1),  # NA for contour lines, symbols for points
         fill = c("darkgray", "gray", "lightgray", 
                  NA, NA),  # Fill for contours, NA for points
         border = c("darkgray", "gray", "lightgray", 
                    NA, NA))  
  
  # Add annotations at the bottom of the plot
  mtext(sprintf("N filled studies: %d", n_filled_studies), 
        side = 1, line = 5, adj = 0)
  mtext(paste0("Meta-analysis RR (95% CI) = ", original_CI), 
        side = 1, line = 7, adj = 0)
  mtext(paste0("Trim-and-fill RR (95% CI) = ", random_CI), 
        side = 1, line = 6, adj = 0)
  if(n_filled_studies != 0) {
    mtext(sprintf("Difference in RRs = %.2f", change_in_or), 
          side = 1, line = 8, adj = 0)
  }
  else{
    mtext(paste0("Difference in OR = N/A"), 
          side = 1, line = 8, adj = 0)
  }
  
  mtext(paste0("Egger's test: Intercept (95% CI) = ", egger_intercept, 
               " (", egger_ci, "), p = ", round(egger_p, 3)), 
        side = 1, line = 9, adj = 0)
  
  dev.off()
  
}

# ---- Exclusions for sensitivity analyses  ----

# Inclusions are based on :
#   1) Studies that adjusted for smoking
#   2) Traditional cohort studies
#   3) Large health databases
#   4) High quality studies only (bias score == 2)


cancer_list <- unique(cancer_allsites$`Top 25 cancer`)

# Function to perform meta-analysis and format results
perform_meta_analysis <- function(data, cancer_type) {
  
  metares <- metagen(
    TE = log(RR),
    lower = log(LCI),
    upper = log(UCI),
    level.ci = 0.95,
    data = data,
    sm = "RR",
    random = TRUE,
    common = FALSE,
    method.tau = "PM"
  )
  
  # Calculate summary statistics
  risk_sum <- round(exp(metares[["TE.random"]]), 2)
  lci_sum <- round(exp(metares[["lower.random"]]), 2)
  uci_sum <- round(exp(metares[["upper.random"]]), 2)
  n_cohort <- sum(data$`N cohorts`)
  
  # Handle special case for head and neck cancers
  n_cases <- sum(data$`N cases`, na.rm = TRUE)
  if (cancer_type == "Head and neck (never-smokers)") {
    n_cases <- n_cases + 796
  }
  
  # Create formatted results table
  data.frame(
    Site = cancer_type,
    Analysis_type = "Main analysis", 
    N_cases = n_cases,
    N_cohorts = n_cohort,
    OR = risk_sum,
    LCI = lci_sum,
    UCI = uci_sum,
    `RR (95% CI)` = sprintf("%.2f (%.2f, %.2f)", risk_sum, lci_sum, uci_sum),
    `P-value` = metares[["pval.random"]],
    stringsAsFactors = FALSE
  )
}

# ---- Build sensitivity results ----
sensitivity_results <- data.frame()

for (cancer_type in cancer_list) {
 
  # Get data for current cancer type
  cancer_data <- cancer_allsites %>%
    subset(`Top 25 cancer` == cancer_type) %>% 
    mutate(`N cases` = as.numeric(`N cases`))
  
  # 1. Main analysis
  main_result <- perform_meta_analysis(cancer_data, cancer_type)
  main_result$Analysis_type <- "Main analysis"
  sensitivity_results <- rbind(sensitivity_results, main_result)
  
  # 2. Smoking-adjusted analysis
  smok_adj <- subset(cancer_data, smok_indicator == 1)
  if (nrow(smok_adj) > 0) {
    smok_result <- perform_meta_analysis(smok_adj, cancer_type)
    smok_result$Analysis_type <- "Smoking-adjusted only"
    sensitivity_results <- rbind(sensitivity_results, smok_result)
  }
  
  # 3. Traditional cohorts analysis
  trad_cohorts <- subset(cancer_data, health_dat == 0)
  if (nrow(trad_cohorts) > 0) {
    cohort_result <- perform_meta_analysis(trad_cohorts, cancer_type)
    cohort_result$Analysis_type <- "Traditional cohorts only"
    sensitivity_results <- rbind(sensitivity_results, cohort_result)
  }
  
  # 4. health care datasets analysis
  health_dat_cohorts <- subset(cancer_data, health_dat == 1)
  if (nrow(health_dat_cohorts) > 0) {
    health_dat_result <- perform_meta_analysis(health_dat_cohorts, cancer_type)
    health_dat_result$Analysis_type <- "Health cohorts only"
    sensitivity_results <- rbind(sensitivity_results, health_dat_result)
  }
  
  # 5. High quality studies
  hi_quality <- subset(cancer_data, `Final bias grade (0-2)` == 2)
  if (nrow(hi_quality) > 0) {
    high_result <- perform_meta_analysis(hi_quality, cancer_type)
    high_result$Analysis_type <- "Low bias studies only"
    sensitivity_results <- rbind(sensitivity_results, high_result)
  }
}

# ---- Format final sensitivity table ----
results_formatted <- sensitivity_results %>%
  mutate(Site = factor(Site, levels = cancer_list),
         PValue_sum = ifelse(P.value < 0.001, "< 0.001", 
                             round(P.value, 3)),
         N_cases = case_when(Site == "Head and neck (never-smokers)" &
                             (Analysis_type == "Health cohorts only" | 
                             Analysis_type == "Low bias studies only") ~ N_cases - 796,   # Remove additional case counts for head and neck cancer where the pooling study was not included
                             TRUE ~ N_cases),
         N_cases = format(N_cases, big.mark = ",", scientific = FALSE)) %>%
  group_by(Site) %>%
  mutate(
    # Get the main analysis RR
    main_RR = OR[Analysis_type == "Main analysis"],
    # Calculate delta for each other analysis
    delta_RR = case_when(
      Analysis_type == "Main analysis" ~ NA,
      Analysis_type == "Health cohorts only" ~ OR - main_RR,
      Analysis_type == "Low bias studies only" ~ OR - main_RR,
      Analysis_type == "Smoking-adjusted only" ~ OR - main_RR,
      Analysis_type == "Traditional cohorts only" ~ OR - main_RR,
      TRUE ~ NA_real_
    )
  ) %>%
  arrange(Site) %>%
  select(Site, Analysis_type, N_cases, N_cohorts, "RR..95..CI.", PValue_sum, delta_RR) 



# ---- Test for generalisability by region ----
# Helper: consistent rounding & p-value formatting
fmt_p <- function(p) {
  ifelse(
    is.na(p), NA,
    ifelse(p < 0.001, "< 0.001",
           ifelse(p < 0.01,
                  # show 3 decimals for very small values (but ≥0.001)
                  sprintf("%.3f", round(p, 3)),
                  # show 2 decimals otherwise
                  sprintf("%.2f", round(p, 2))
           )
    )
  )
}

## create list with new order of cancer sites based on pooled RRs
cancer_list_ordered <- results %>% arrange(desc(OR)) %>% pull(Site)

# ---- Prepare data ----
cancer_region <- cancer_allsites %>% 
  mutate(`N cases` = case_when(is.na(`N cases`) ~ 0, TRUE ~ `N cases` )) # Set N cases to 0 if missing

# ---- Regional heterogeneity results ----
het_res <- data.frame()

for (c in cancer_list_ordered) {
  
  cancer_dta <- cancer_region %>%
    subset(`Top 25 cancer` == c) 
  
  metares <- metagen(TE = log(RR),  
                     lower=log(LCI),  
                     upper = log(UCI),  
                     level.ci = 0.95,  
                     data = cancer_dta,  
                     subgroup = RegionG,
                     method.tau = "PM"
  ) 
  
  # Back-transform subgroup pooled results
  Risk_sum <- round(exp(metares[["TE.random.w"]]),2) 
  LCI_sum <- round(exp(metares[["lower.random.w"]]),2) 
  UCI_sum <-  round(exp(metares[["upper.random.w"]]),2) 
  
  # Format P-values
  PValue_sum <- fmt_p(as.numeric(metares[["pval.random.w"]]))

  # Format P-het
  Het_grps <- metares[["pval.Q.b.random"]]
  Het_grps <- pvalr(Het_grps)
  
  # Subgroup names (match ordering of the pooled vectors)
  Region <- metares[["bylevs"]]
  
  # Case counts by subgroup (align by RegionG order)
  n_cases_by_region <- cancer_dta %>%
    group_by(RegionG) %>%
    summarise(N_cases_region = sum(`N cases`, na.rm = TRUE), .groups = "drop") %>%
    arrange(RegionG)
  
  # Results dataframe
  temp <- data.frame(Site = paste(c), 
                     Region = Region, 
                     N_cases = n_cases_by_region$N_cases_region[match(Region, n_cases_by_region$RegionG)], 
                     OR = Risk_sum, 
                     LCI = LCI_sum, 
                     UCI = UCI_sum, 
                     `RR (95% CI)` = paste0(Risk_sum," (",LCI_sum,", ", UCI_sum,")"),
                     `P-value` = PValue_sum, 
                     Phet = Het_grps)
  
  het_res <- plyr::rbind.fill(temp,het_res)
  
}

fwrite(het_res, "het_region_per_cancer.csv")


# ---- Forest plots for individual studies by cancer ----
# Supplementary Figures 2-26

# format N cases with commas to improve readability 
cancer_allsites <- cancer_allsites %>% 
  mutate(`N cases` = format(`N cases`, big.mark=",", scientific=FALSE))

# Loop over each cancer site in pre-defined order
for (c in cancer_list_ordered){
  
  # ---- Prepare study-level data ----
  cancer_dta <- cancer_allsites %>%
    subset(Forest_cancer == c) %>% 
    arrange(Date) %>%
    mutate(studylab= paste0(Author, " (", Date,")"),
           "HRCI"= paste0(round(OR, 2), " (", round(Lci, 2),", ",round(Uci, 2),")")
    )
  
  # ---- Run meta-analysis ----
  metares <- metagen(TE = log(RR),  
                     lower = log(LCI),  
                     upper = log(UCI),  
                     level.ci = 0.95,  
                     data = cancer_dta,  
                     sm = "RR",
                     random = TRUE, 
                     fixed = FALSE, 
                     comb.random = TRUE, 
                     studlab = studylab,
                     method.tau = "PM") 

  # ---- Add custom study-level metadata for plotting and improve spacing ----
  metares$study <- paste("    ", cancer_dta$Study,  sep = "")
  metares$site <- paste("    ", cancer_dta$`Cancer site`,"    ",  sep = "")
  metares$PMID <- paste("    ", cancer_dta$`PMID`,"    ",  sep = "")
  metares$Ncases <- paste("    ", cancer_dta$`N cases`,"    ",  sep = "")
  
  
  # Compute IQR of study-specific ORs for annotation line
  iqr_hr <- round(quantile(cancer_dta$OR, probs=c(0.25, 0.5, 0.75)), 2) 
  

  # Create figure  
    png(
      filename = paste0("forestplot_",c,".png"),
      width = 14, 
      height = max(15, nrow(cancer_dta)* -0.8),
      units = "in",
      res = 300, 
      type = "cairo"  
    )

  meta::forest(metares, 
               xlim = c(0.5, 2.5), 
               at=c(0.5, 0.75, 1, 1.5, 2, 2.5),
               xlab = expression("RR per 5 kg/m"^2~"increase"),
               
               leftcols = c("studylab", "site", "Ncases", "Country"),
               leftlabs = c("Author", "    Site", "N cases", "Country/Region"),
               
               rightlabs=c("RR", "[95% CI]", "   Cohort", "PMID"),
               rightcols=c("effect", "ci", "study", "PMID"),
               
               comb.random = TRUE,
               just="left",   
               just.addcols="left", 
               addrow.overall = TRUE,
               addrows.below.overall = 3,
               
               colgap = "0.2cm",
               col.diamond.random	= "lightblue",
               col.square = "darkblue",
               col.square.lines = "darkblue",
               # Text annotation
               text.addline1 = paste0("Median RR (25th, 75th percentile) = ", 
                                      formatC(iqr_hr[[2]], digits = 2, format = "f"), 
                                      " (", formatC(iqr_hr[[1]], digits = 2, format = "f"),", ", 
                                      formatC(iqr_hr[[3]], digits = 2, format = "f"),")")
  )
  
  dev.off()
}


