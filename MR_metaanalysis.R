# =============================================================================
# Title: Meta-analysis of BMI–Cancer MR Studies
# =============================================================================
# Overview
# - Reads MR results from an Excel sheet and runs per-cancer random-effects
#   meta-analyses on the log(RR) scale using {meta}.
# - When ≥2 studies are available for a cancer site: fits Paule–Mandel tau^2
# - When exactly 1 study is available: reports that study’s RR and 95% CI
#
# Inputs
# - File: "SupplementaryData.xlsx"
# - Sheet: "BMI_cancer_MR"
# =============================================================================


# Load packages 
pacman::p_load("tidyverse", "data.table", "janitor", "meta", "readxl", "dmetar")

# p-value formatting function
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

# round to 2 decimal places unless value is < 0.01, else round to the last non zero decimal 
# I use this function if p-values aren't v. small
my_round = function(x, n=3) {
  max(abs(round(x, n)), abs(signif(x, 1))) * sign(x)
} 

# Compute a two-sided p-value from a study's RR and 95% CI on the log scale
p_from_ci <- function(or, lci, uci) {
  if (any(!is.finite(c(or, lci, uci))) || or <= 0 || lci <= 0 || uci <= 0) return(NA_real_)
  se <- (log(uci) - log(or)) / qnorm(0.975)
  if (!is.finite(se) || se <= 0) return(NA_real_)
  z  <- abs(log(or) / se)
  2 * (1 - pnorm(z))
}

# Read in data 
cancer_allsites <- read_excel("SupplementaryData.xlsx", 
                              sheet = "BMI_cancer_MR") 

# Promote first row to header
colnames(cancer_allsites) <- as.character(unlist(cancer_allsites[1, ]))
cancer_allsites <- cancer_allsites[-1, ]

# Make sure values are numeric
cancer_allsites <- cancer_allsites %>%
  mutate(OR = as.numeric(Risk),
         Lci = as.numeric(LCI),
         Uci = as.numeric(UCI),
         N_cases = as.numeric(`N cases`))

# Derive analysis list of cancer sites
cancer_list <- unique(cancer_allsites$"Top 25 cancer")
cancer_list <- cancer_list[!is.na(cancer_list)]


# ---- Primary meta-analysis of MR studies -------------------------------------

results <- data.frame()

for (c in cancer_list){
  
  # Subset to one site
  cancer_dta <- cancer_allsites %>%
    subset(`Top 25 cancer` == c) 
  
  n_studies_avail <- nrow(cancer_dta)
  
  if (n_studies_avail >= 2L) {
    
  # ---- standard random-effects meta-analysis ----
  metares <- metagen(TE = log(OR),  
                     lower=log(Lci),  
                     upper = log(Uci),  
                     level.ci = 0.95,  
                     data = cancer_dta,  
                     sm="RR",
                     random = TRUE, 
                     common = FALSE, 
                     method.tau = "PM"
                     ) 
  
  # Back-transform pooled effect & CI to RR scale for reporting
  Risk_sum <- round(exp(metares[["TE.random"]]), 2) 
  LCI_sum <- round(exp(metares[["lower.random"]]), 2) 
  UCI_sum <-  round(exp(metares[["upper.random"]]), 2) 
  
  # Study count: unique PMIDs (robust to multiple rows per study)
  N_studies <- table(cancer_dta$`Study PMID`) %>% nrow()
  
  # Total cases
  N_cases <-  sum(cancer_dta$`N_cases`, na.rm = TRUE)
  
  # P-value formatting
  PValue_sum <- pvalr(metares[["pval.random"]])
  
  # I^2 (% heterogeneity) rounded when available
  i2 <- metares[["I2"]]

    if (!is.na(i2)) {
    i2 <- round(i2, 2)
    }
  } else if (n_studies_avail == 1L) {
   
   # ---- if single-study ----
    s <- cancer_dta[1, ]
    
    Risk_sum <- round(as.numeric(s$OR),  2)
    LCI_sum  <- round(as.numeric(s$Lci), 2)
    UCI_sum  <- round(as.numeric(s$Uci), 2)
    
    N_studies <- 1L
    N_cases   <- sum(s$N_cases, na.rm = TRUE)
    
    # Compute p-value from the CI (log scale)
    p_single <- p_from_ci(or = as.numeric(s$OR),
                          lci = as.numeric(s$Lci),
                          uci = as.numeric(s$Uci))
    PValue_sum <- pvalr(p_single)
    
    i2 <- NA_real_  # not defined for a single study
  }
  
  # Total cases pretty-print with commas
  N_cases_formatted <- format(N_cases, big.mark = ",", scientific = FALSE)
  
  
  # Assemble one row for this site
  temp <- data.frame(Site = paste(c), 
                     N_cases = N_cases_formatted, 
                     N_studies = N_studies, 
                     OR = Risk_sum, 
                     LCI = LCI_sum, 
                     UCI = UCI_sum, 
                     `RR (95% CI)` = paste0(formatC(Risk_sum, digits = 2, format = "f"),
                                            " (",formatC(LCI_sum, digits = 2, format = "f"),", ", 
                                            formatC(UCI_sum, digits = 2, format = "f"),")"), 
                     `P-value` = PValue_sum, 
                     I2 = i2)
  
  # Bind together
  results <- rbind(temp,results)
  
}

fwrite(results, "MR_metaanalysis.csv")
