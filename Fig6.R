############################################################################################
# Sex-stratified meta-analysis (men vs women) and forest plot
#
# Purpose:
#   1) Read sex-specific risk estimates by cancer site from the Supplementary Data.
#   2) Reshape to long format (one row per cancer site × sex).
#   3) Run subgroup meta-analysis by sex for each cancer site and extract:
#        - pooled RR (random effects)
#        - p-value for pooled estimate
#        - p-value for between-sex heterogeneity (subgroup difference; P-het)
#        - total cases per sex
#   4) Export results to CSV.
#   5) Build a two-column (men/women) forest plot of pooled estimates by cancer (Fig 6).
#
############################################################################################

# # Function to load packages 
pacman::p_load("forestploter", "tidyverse", "janitor", "data.table", "meta", "readxl")

# Read in data and recode 
cancer_allsites <- read_excel("SupplementaryData.xlsx", 
                              sheet = "Sex_strat_working")

## Pivot to a longform table (1 row per sex)
cancer_allsites_long <- cancer_allsites %>%
  # make names consistent so pivoting is easy
  rename_with(~ str_replace_all(.x, "__+", "_")) %>%   # collapses double underscores
  rename_with(~ str_replace(.x, "_women$", "_women")) %>% # harmless, keeps consistent
  pivot_longer(
    cols = matches("_(men|women)$"),
    names_to = c(".value", "sex"),
    names_pattern = "^(.*)_(men|women)$"
  ) %>%
  # nicer column names
  rename(
    cases = cases,
    analytic_set = `N analytic set`,
    RR = `Risk estimate`,
    LCI = LCI,
    UCI = UCI
  ) %>%
  # drop rows where everything sex-specific is missing
  filter(!(is.na(cases) & is.na(analytic_set) & is.na(RR) & is.na(LCI) & is.na(UCI))) %>%
  # Ensure numeric columns are numeric
  mutate(
    across(c(cases, analytic_set, RR, LCI, UCI), ~ suppressWarnings(as.numeric(.x)))
  ) %>%
  # Drop missing rows
  filter(!is.na(RR)) %>% 
  clean_names()

# List of cancer sites available in the data
cancer_list <- unique(cancer_allsites_long$top_25_cancer)

# Drop oesophageal SQ as no data for men were available
cancer_list <- cancer_list[cancer_list != "Oesophageal (SQ, never-smokers)"]

## ---- Prepare cases and run heterogeneity-by-sex meta-analyses ----

# Coerce cases cleanly and set missing cases to 0 (so sums work)
cancer_sex <- cancer_allsites_long %>% 
  mutate(cases = str_replace_all(cases, "\\s", ""), # remove spaces e.g. "1 234"
         cases = as.numeric(cases),
         cases = case_when(is.na(cases) ~ 0, TRUE ~ cases ))

het_res <- data.frame()
sex_order <- c("men", "women")


# Loop over cancer sites and compute sex-specific pooled estimates +
# between-sex heterogeneity test (subgroup difference).

for (c in cancer_list) {
  
  # Subset to a single cancer site and enforce consistent sex ordering
  cancer_dta <- cancer_sex %>%
    filter(top_25_cancer == c) %>%
    mutate(sex = factor(sex, levels = sex_order))
  
  # Random-effects meta-analysis with sex as a subgroup.
  metares <- metagen(
    TE = log(rr),
    lower = log(lci),
    upper = log(uci),
    level.ci = 0.95,
    data = cancer_dta,
    subgroup = sex,
    method.tau = "PM"
  )
  
  # Extract subgroup labels as returned by meta
  sex <- metares$bylevs
  
  # Safety check: ensure we always have both sexes in the expected order.
  stopifnot(identical(as.character(sex), sex_order))
  
  # Extract subgroup-specific random-effects pooled estimates and exponentiate.
  Risk_sum <- sprintf("%.2f", exp(metares$TE.random.w))
  LCI_sum  <- sprintf("%.2f", exp(metares$lower.random.w))   
  UCI_sum  <- sprintf("%.2f", exp(metares$upper.random.w))  
  
  # P-value for pooled effect within subgroup
   PValue_sum <- ifelse(metares$pval.random.w < 0.001, "< 0.001",
                       round(metares$pval.random.w, 3))
  
   # P-value for between-subgroup heterogeneity (sex difference)
    Het_grps <- ifelse(metares$pval.Q.b.random < 0.001, "< 0.001",
                     round(metares$pval.Q.b.random, 2))
  
    # Total cases per sex (summed within cancer site)
    Ncases <- cancer_dta %>%
      group_by(sex) %>%
      summarise(N_cases_sex = sum(cases, na.rm = TRUE), .groups = "drop") %>%
      arrange(factor(sex, levels = sex_order)) %>%
      pull(N_cases_sex)
  
    # Assemble per-site output in a long format (one row per sex)
    temp <- data.frame(
      Site = c,
      Sex = sex,
      N_cases = Ncases,
      OR = Risk_sum,
      LCI = LCI_sum,
      UCI = UCI_sum,
      `RR (95% CI)` = paste0(Risk_sum, " (", LCI_sum, ", ", UCI_sum, ")"),
      `P-value` = PValue_sum,
      Phet = Het_grps
      )
  # Append results
  het_res <- plyr::rbind.fill(het_res, temp)
}


## ---- Forest plot data prep ----

summary_sex <- het_res %>% 
  # Ensure vars are numeric
  mutate(
    across(
      c(N_cases, OR, LCI, UCI),
      ~ suppressWarnings(as.numeric(.x))
    )
  ) %>% 
  rename(`RR (95% CI)`="RR..95..CI.") %>%
  mutate(`N_cases` = round(`N_cases`,0)) %>%
  select("Site" , "Sex", "N_cases", "OR",  "LCI", "UCI", "RR (95% CI)" ,"Phet" ) %>%
  # Wide format so each cancer is one row with men/women columns side-by-side
  tidyr::pivot_wider(names_from = Sex, values_from = c(OR, LCI, UCI, `RR (95% CI)`, `N_cases`)) %>%
  # Spacer columns used by forestploter to control layout
    mutate(` ` = paste(rep(" ", 50), collapse = " "), 
         `  ` = paste(rep(" ", 1), collapse = " "),
         N_cases_men = format(N_cases_men, big.mark=",", scientific=FALSE),
         `N_cases_women` = format(`N_cases_women`, big.mark=",", scientific=FALSE))  %>%
  # set column order
  dplyr::select("Site", ` `, "N_cases_men", "OR_men", "LCI_men", "UCI_men", "RR (95% CI)_men",
                "N_cases_women", "OR_women", "LCI_women", "UCI_women", "RR (95% CI)_women",
                `  `, "Phet") %>%
  ## remove Head and neck as insufficient case numbers 
  filter(Site!="Head and neck (never-smokers)" ) %>%
  
  ### Add spaces to RR 95% CIs to improve spacing
  mutate(`RR (95% CI)_men` = paste(`RR (95% CI)_men`, "           "),
         `RR (95% CI)_women` = paste0(`RR (95% CI)_women`, "           " )) %>%
  
  # Order cancers for plotting (here by men’s pooled RR)
  arrange(desc(OR_men)) %>%
  # Format column headers
  rename("Cancer type" = "Site", 
         "Women \nN cases"="N_cases_women",
         "Men \nN cases"="N_cases_men",
         "Women \nRR (95% CI)"="RR (95% CI)_women",
         "Men \nRR (95% CI)"="RR (95% CI)_men",
         "P-het"="Phet") %>%
  # Approximate standard errors from 95% CI on the log scale:
  mutate(se_women = (log(UCI_women) -log(LCI_women)) / 3.92,
         se_men = (log(UCI_men) -log(LCI_men)) / 3.92) %>%
  # Convert to an inverse-SE proxy for point sizes; scaled down to avoid huge dots
  mutate(se_women_inv = ((1/se_women)/100),
         se_men_inv = ((1/se_men)/100)) 

## ---- Add blank spacer rows to improve readability ----

# Create a blank row with the same columns; keep numeric plotting columns as NA
blank_row <- data.frame(matrix("", ncol = ncol(summary_sex)))
colnames(blank_row) <- colnames(summary_sex)

blank_row <- blank_row %>%
  mutate(across(contains("OR_"), ~NA),
         across(contains("LCI_"), ~NA),
         across(contains("UCI_"), ~NA),
         across(contains("se_"), ~NA))

# Interleave a blank row after each data row
new_data <- list()

# Loop through the original data frame and add a blank row after every row
for (i in 1:nrow(summary_sex)) {
  new_data[[i * 2 - 1]] <- summary_sex[i, ]
  new_data[[i * 2]] <- blank_row
}

final_data <- do.call(rbind, new_data)
final_data <- final_data[-nrow(final_data), ]

## ---- Forest plot ----

# Set-up theme
tm <- forest_theme(refline_lty = c("solid"),
                   refline_col = c("#bababa"),
                   base_size = 15, 
                   legend_name = "Sex",
                   legend_value = c(" Men", "Women"),
                   ci_col = c("#0072B2", "#D55E00"),
                   core=list(bg_params=list(fill = c("white"))))

# Save to pdf
pdf(file = "Forrest_plot_sex_test.pdf", width = 20, height = 30)

forestploter::forest(final_data[,c(1:3,7:8,12:14)],
                     est = list(
                       final_data$`OR_men`, final_data$`OR_women`),
                     lower = list(
                       final_data$`LCI_men`,final_data$`LCI_women`),
                     upper = list(
                       final_data$`UCI_men`, final_data$`UCI_women`),
                     sizes = list(
                       final_data$se_men_inv, final_data$se_women_inv),
                     ci_column = c(2),
                     xlim = c(0.5, 2.0),
                     ticks_at = c( 0.5, 1, 1.5, 2),
                     ref_line = 1,
                     nudge_y = 0.4,
                     xlab = expression("RR per 5 kg/m"^2~"increase"),
                     theme = tm)

dev.off()
