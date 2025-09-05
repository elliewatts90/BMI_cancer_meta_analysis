# =============================================================================
# Figure 6: Forest plot comparing BMI vs Waist circumference pooled RRs
#
# Inputs:
#   - wc_bmi_compare_HRs_wald.csv  (from compare_waist_bmi.R)
#
# Output:
#   - Fig6.pdf
# =============================================================================


# ---- Setup -------------------------------------------------------------------
pacman::p_load("data.table", "tidyverse",  "forestploter")



# ---- Read & split BMI vs Waist rows ------------------------------------------
res_bmi<-fread("wc_bmi_compare_HRs_wald.csv") %>% 
  filter(Measurement == "BMI") %>% 
  rename(
    OR_bmi = OR, 
    LCI_bmi = LCI, 
    UCI_bmi = UCI, 
    `BMI\nN cases`= N_cases_bmi
    ) %>%
  dplyr::select(
    Site, 
    `BMI\nN cases`, 
    OR_bmi, 
    LCI_bmi, 
    UCI_bmi, 
    Phet
    )

res_full<-fread("wc_bmi_compare_HRs_wald.csv") %>% 
  filter(Measurement == "Waist")  %>% 
  rename(
    OR_wc = OR, 
    LCI_wc = LCI, 
    UCI_wc = UCI, 
    `Waist circumference\nN cases`= N_cases_waist) %>%
  dplyr::select(
    c(
      Site, 
      `Waist circumference\nN cases`, 
      OR_wc, 
      LCI_wc,
      UCI_wc,
      Phet
      )
    ) %>%
  full_join(res_bmi, by="Site") # Join BMI and Waist on Site; resolve duplicated Phet to one column


# ---- Build display table for forestploter ------------------------------------

forest_dat <- res_full %>%  
  mutate(`BMI\nRR (95% CI)`=ifelse(is.na(res_full$OR_bmi), "",
                                sprintf("%.2f (%.2f, %.2f)   ",
                                        res_full$OR_bmi, 
                                        res_full$LCI_bmi, 
                                        res_full$UCI_bmi))) %>% 
  mutate(`Waist circumference\nRR (95% CI)`= ifelse(is.na(res_full$OR_wc), "",
                                           sprintf("%.2f (%.2f, %.2f)            ",
                                                   res_full$OR_wc, 
                                                   res_full$LCI_wc, 
                                                   res_full$UCI_wc))) %>% 
  mutate(` `=paste(rep(" ", 50), collapse = " ") ) %>% 
  dplyr:: select(c(Site, `BMI\nN cases`, `Waist circumference\nN cases`, OR_bmi, LCI_bmi, UCI_bmi, OR_wc, LCI_wc, UCI_wc, ` `, 
                   `BMI\nRR (95% CI)`, `Waist circumference\nRR (95% CI)`, Phet.x)) %>%
  rename(`P-het` =  Phet.x) %>%
  # Pretty-print case counts (keep underlying numerics intact elsewhere)
  
  mutate(`BMI\nN cases` = format(`BMI\nN cases`, big.mark=",", scientific=FALSE),
         `Waist circumference\nN cases` = format(`Waist circumference\nN cases`, big.mark=",", scientific=FALSE)) %>%

  arrange(desc(OR_bmi)) %>% 
  
  # Inverse-SE size proxies (95% CI width = 3.92 * SE on log scale)
  mutate(se_bmi = (log(UCI_bmi) -log(LCI_bmi)) / 3.92,
         se_wc = (log(UCI_wc) -log(LCI_wc)) / 3.92,
         se_inv_bmi = (1/se_bmi)/35,
         se_inv_wc = (1/se_wc)/35) 

# ---- Insert blank spacer rows to add vertical breathing room -----------------

blank_row <- data.frame(matrix("", ncol = ncol(forest_dat)))
colnames(blank_row) <- colnames(forest_dat)

blank_row <- blank_row %>%
  mutate(across(contains("OR_"), ~ NA_real_),
         across(contains("LCI_"), ~ NA_real_),
         across(contains("UCI_"), ~ NA_real_),
         across(contains("se_"), ~ NA_real_))

# Initialize a list to store the modified rows
new_data <- list()

# Loop through the original data frame and add a blank row after every row
for (i in 1:nrow(forest_dat)) {
  new_data[[i * 2 - 1]] <- forest_dat[i, ]
  new_data[[i * 2]] <- blank_row
}

final_data <- do.call(rbind, new_data)
final_data <- final_data[-nrow(final_data), ]

# ---- Theme for forestploter ---------------------------------------------------

tm <- forest_theme(refline_lty = c("solid"),
                   refline_col = c("#bababa"),
                   base_size = 16, 
                   ci_lwd= 3,
                   legend_name = "Measurement",
                   legend_value = c(" BMI", " Waist circumference"),
                   ci_col = c("darkblue", "darkgreen"),
                   core=list(bg_params=list(fill = c("white"))),
                   core=list(fg_params=list(hjust=0, x=0)))

# ---- Draw and export ----------------------------------------------------------

pdf(file = "Fig6.pdf", width = 20, height = 20)

forest(final_data[,c(1:3, 10:13)],
       est = list(final_data$OR_bmi,
                  final_data$OR_wc),
       lower = list(final_data$LCI_bmi,
                    final_data$LCI_wc), 
       upper = list(final_data$UCI_bmi,
                    final_data$UCI_wc),
       sizes = list(final_data$se_inv_bmi,
                    final_data$se_inv_wc),
       ci_column = c(4),
       xlim = c(0.5, 2.0),
       ticks_at = c( 0.5, 1, 1.5, 2),
       ref_line = 1,
       nudge_y = 0.6, 
       xlab = "RR per 1 SD increase",
       theme = tm)


dev.off()



