# =============================================================================
# Title: Figure 5: Regional summary forest plot (Europe vs North America vs East Asia)
#
# Requires observational meta-analysis to already be run
# =============================================================================

# ---- Packages ----------------------------------------------------------------

pacman::p_load("data.table", "tidyverse", "forestploter", "janitor")

# ---- Set minimum case numbers per cancer to display in figure ----------------

case_minimum <- 500


# ---- Build wide table: Site × {Europe, North America, East Asia} -------------

summary_region <-fread("het_region_per_cancer.csv") %>% 
  # Keep reasonably powered rows first
  filter(N_cases > case_minimum) %>%
  # Fix over-rounding of RRs
  mutate(
    `RR (95% CI)` = dplyr::case_when(
      is.na(OR) | is.na(LCI) | is.na(UCI) ~ NA_character_,
      TRUE ~ sprintf("%.2f (%.2f, %.2f)", as.numeric(OR), as.numeric(LCI), as.numeric(UCI))
    ),
    `N_cases` = round(`N_cases`, 0)
  ) %>%
  dplyr::select(-c(P.value, "RR..95..CI."))   %>%

  # Pivot to one row per Site with per-region columns
  tidyr::pivot_wider(
    names_from = Region, 
    values_from = c(OR, LCI, UCI, `RR (95% CI)`, `N_cases`)
    ) %>%
  
  # Spacer columns improve readable gaps in forestploter tables
  mutate(
    ` ` = paste(rep(" ", 50), collapse = " "), 
         `  ` = paste(rep(" ", 1), collapse = " "),
    
  # Pretty-print counts with commas
         N_cases_Europe = format(N_cases_Europe, big.mark=",", scientific=FALSE),
         `N_cases_North America` = format(`N_cases_North America`, big.mark=",", scientific=FALSE),
         `N_cases_East Asia` = format(`N_cases_East Asia`, big.mark=",", scientific=FALSE)
  )  %>%
  
  dplyr::select("Site", ` `, 
                "N_cases_Europe", "OR_Europe", "LCI_Europe", "UCI_Europe", "RR (95% CI)_Europe",
                "N_cases_North America", "OR_North America", "LCI_North America", "UCI_North America", "RR (95% CI)_North America",
                "N_cases_East Asia", "OR_East Asia", "LCI_East Asia", "UCI_East Asia","RR (95% CI)_East Asia",
                `  `, "Phet"
                ) %>%
  
  # Keep only sites with all three regional estimates
  filter(!is.na(`OR_North America`) & !is.na(`OR_East Asia`) & !is.na(`OR_Europe`)) %>%
  
  # Add a few trailing spaces to RR strings to avoid crowding in the table
  mutate(
    `RR (95% CI)_Europe` = paste(`RR (95% CI)_Europe`, "           "),
    `RR (95% CI)_North America` = paste0(`RR (95% CI)_North America`, "           " ),
    `RR (95% CI)_East Asia` = paste0(`RR (95% CI)_East Asia`, "            ")
    ) %>%
  arrange(desc(OR_Europe)) %>%
  # Friendly headers for printed columns (with newlines to stack)
  
  rename(
    "N Am \nN cases"="N_cases_North America",
    "Europe \nN cases"="N_cases_Europe",
    "East Asia \nN cases"="N_cases_East Asia",
    "North America \nRR (95% CI)"="RR (95% CI)_North America",
    "Europe \nRR (95% CI)"="RR (95% CI)_Europe",
    "East Asia \nRR (95% CI)"="RR (95% CI)_East Asia",
    "LCI_EA"="LCI_East Asia", 
    "UCI_EA"="UCI_East Asia" ,
    "LCI_NA"="LCI_North America", 
    "UCI_NA"="UCI_North America",
    "P-het"="Phet"
    ) %>%
  
  # Compute SEs from CIs on log scale 
   mutate(se_NA = (log(UCI_NA)      -log(LCI_NA))     / 3.92,
          se_EU = (log(UCI_Europe)  -log(LCI_Europe)) / 3.92,
          se_EA = (log(UCI_EA)      -log(LCI_EA))     / 3.92, 
          
  # set inverse SE for point sizing
          se_NA_inv = ((1/se_NA)/100),
          se_EU_inv = ((1/se_EU)/100),
          se_EA_inv = ((1/se_EA)/100)) 


# ---- Insert blank spacer rows for vertical breathing room --------------------
# This doubles rows and inserts a blank after each real row.
blank_row <- data.frame(matrix("", ncol = ncol(summary_region)))
colnames(blank_row) <- colnames(summary_region)
blank_row <- blank_row %>%
  mutate(across(contains("OR_"), ~ NA),
         across(contains("LCI_"), ~ NA),
         across(contains("UCI_"), ~ NA),
         across(contains("se_"), ~ NA))

# Initialize a list to store the modified rows
new_data <- list()
# Loop through the original data frame and add a blank row after every row
for (i in 1:nrow(summary_region)) {
  new_data[[i * 2 - 1]] <- summary_region[i, ]
  new_data[[i * 2]] <- blank_row
}
final_data <- do.call(rbind, new_data)
final_data <- final_data[-nrow(final_data), ]


# ---- Forest plot theme --------------------------------------------------------
tm <- forest_theme(refline_lty = c("solid"),
                   refline_col = c("#bababa"),
                   base_size = 15, 
                   legend_name = "Region",
                   legend_value = c(" Europe", " North America",  " East Asia"),
                   ci_col = c("darkblue", "darkgreen", "purple"),
                   core=list(bg_params=list(fill = c("white"))))


# ---- Draw forest plot (log-scaled x-axis) ------------------------------------

pdf(file = paste0("Fig5.pdf"), width = 50, height = 30)

forest(final_data[,c(1:3,7:8,12:13, 17:19)],
       est = list(
         final_data$`OR_Europe`, final_data$`OR_North America`,
         final_data$`OR_East Asia`),
       lower = list(
         final_data$`LCI_Europe`,final_data$`LCI_NA`,
         final_data$`LCI_EA`), 
       upper = list(
         final_data$`UCI_Europe`, final_data$`UCI_NA`,
         final_data$`UCI_EA`),
       sizes = list(
         final_data$se_EU_inv, final_data$se_NA_inv,
         final_data$se_EA_inv),
       ci_column = c(2),
       xlim = c(0.47, 2.0),
       x_trans = c("log"),
       ticks_at = c( 0.5, 1, 1.5, 2),
       ref_line = 1,
       nudge_y = 0.4, 
       xlab = expression("RR per 5 kg/m"^2~"increase"),
       theme = tm)

dev.off()
