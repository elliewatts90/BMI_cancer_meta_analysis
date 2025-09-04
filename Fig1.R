# =============================================================================
# Title: Figure 1: Global map of cohorts
# Study population based on the maximum analytic size for each distinct cohort
# N cohorts based on the number of distinct cohorts for each country
# =============================================================================

# load required packages
pacman::p_load("tidyverse", "data.table", "ggplot2", "readxl", "janitor", 
               "hrbrthemes", "viridis", "ggrepel", "ggthemes","devEMF" , 
               "maps", "mapproj")

cancer_allsites <- read_excel("SupplementaryData.xlsx", 
                              sheet = "BMI_cancer_data") 

# Move first row to column names
colnames(cancer_allsites) <- as.character(unlist(cancer_allsites[1, ]))
cancer_allsites <- cancer_allsites[-1, ] # Remove the old first row


# ---- Identify pooled cohort projects ----------------------------------------
# These labels represent pooled consortia/projects (not single cohorts).
pooled_cohorts <-c("Japan cohorts", "ODDS", "CONOR", "ACC", "ANZDCC", "KRIS, MIHDPS", "US cohorts", "NHS, HPFS", "NHS, NHS-II", "LCCIPP")

# Standardize Study names by stripping any parenthetical suffixes (e.g., (KNHIS (women)))
cancer_allsites$Study <- sub("\\s*\\(.*", "", cancer_allsites$Study)


# -----------------------------------------------------------------------------
# In several pooled analyses the underlying cohort-level rows are not present.
# To count them correctly on the map, we manually add rows for constituent
# cohorts with approximate analytic N when available from publications.
# -----------------------------------------------------------------------------

cohort_sample_size <- cancer_allsites %>%
  clean_names() %>%
  mutate(n_analytic_set = str_replace_all(n_analytic_set, "\\s", "")) %>% # remove whitespace
  mutate(n_analytic_set = as.numeric(n_analytic_set),
         n_cases = as.numeric(n_cases),
         date = as.character(year)) %>%

  # --- Manually expand pooled projects into (approximate) constituent cohorts where not already present ---
  
  # 1) Japan cohorts
  add_row(study = "OSAKA", n_analytic_set = 33893, date = "", pmid="37013939", country="Japan", region="East Asia" ) %>%
  add_row(study = "MIYAGI-I", n_analytic_set = 100610, date = "", pmid="", country="Japan", region="East Asia" ) %>%
  add_row(study = "LSS", n_analytic_set = 17059, date = "", pmid="37013939", country="Japan", region="East Asia" ) %>%

  # 2) ODDS (Sweden)
  add_row(study = "SMCR", n_analytic_set = 1754896 , date = "", pmid="", country="Sweden", region="Europe" ) %>%
  add_row(study = "SMBR", n_analytic_set = 1728111 , date = "", pmid="", country="Sweden", region="Europe" ) %>%
  add_row(study = "SCWC", n_analytic_set = 262622 , date = "", pmid="", country="Sweden", region="Europe" ) %>%
  add_row(study = "SMC", n_analytic_set = 56827 , date = "", pmid="", country="Sweden", region="Europe" ) %>%
  add_row(study = "COSM", n_analytic_set = 33802 , date = "", pmid="", country="Sweden", region="Europe" ) %>%
  add_row(study = "MSP", n_analytic_set = 25700 , date = "", pmid="", country="Sweden", region="Europe" ) %>%
  add_row(study = "MONICA", n_analytic_set = 7500 , date = "", pmid="", country="Sweden", region="Europe" ) %>%
  add_row(study = "VIP", n_analytic_set = 67000 , date = "", pmid="", country="Sweden", region="Europe" ) %>%
  add_row(study = "MPP", n_analytic_set = 20971 , date = "", pmid="", country="Sweden", region="Europe" ) %>%
  add_row(study = "MDCS", n_analytic_set = 11622 , date = "", pmid="", country="Sweden", region="Europe" ) %>%
  add_row(study = "MOS", n_analytic_set = 4700 , date = "", pmid="", country="Sweden", region="Europe" ) %>%
  add_row(study = "WLHS", n_analytic_set = 27257 , date = "", pmid="", country="Sweden", region="Europe" ) %>%
  add_row(study = "SNMC", n_analytic_set = 24866 , date = "", pmid="", country="Sweden", region="Europe" ) %>%
  add_row(study = "WSAS", n_analytic_set = 18561 , date = "", pmid="", country="Sweden", region="Europe" ) %>%
  add_row(study = "WICTORY", n_analytic_set = 14864 , date = "", pmid="", country="Sweden", region="Europe" ) %>%
  add_row(study = "MSS", n_analytic_set = 12743 , date = "", pmid="", country="Sweden", region="Europe" ) %>%
  add_row(study = "Lifegene", n_analytic_set = 11219 , date = "", pmid="", country="Sweden", region="Europe" ) %>%
  add_row(study = "Epihealth", n_analytic_set = 7844 , date = "", pmid="", country="Sweden", region="Europe" ) %>%
  add_row(study = "OLDN", n_analytic_set = 6576 , date = "", pmid="", country="Sweden", region="Europe" ) %>%
  
  # 7 Swedish Twin Registry components in ODDS; split total N ~67,142 equally across 7
  add_row(study = "STR1", n_analytic_set = 9591.714 , date = "", pmid="", country="Sweden", region="Europe" ) %>%
  add_row(study = "STR2", n_analytic_set = 9591.714 , date = "", pmid="", country="Sweden", region="Europe" ) %>%
  add_row(study = "STR3", n_analytic_set = 9591.714 , date = "", pmid="", country="Sweden", region="Europe" ) %>%
  add_row(study = "STR4", n_analytic_set = 9591.714 , date = "", pmid="", country="Sweden", region="Europe" ) %>%
  add_row(study = "STR5", n_analytic_set = 9591.714 , date = "", pmid="", country="Sweden", region="Europe" ) %>%
  add_row(study = "STR6", n_analytic_set = 9591.714 , date = "", pmid="", country="Sweden", region="Europe" ) %>%
  add_row(study = "STR7", n_analytic_set = 9591.714 , date = "", pmid="", country="Sweden", region="Europe" ) %>%
  
  # 3) CONOR (Norway)
  add_row(study = "Tromso4", n_analytic_set = 26925 , date = "", pmid="17984119", country="Norway", region="Europe" ) %>%
  add_row(study = "HUSK", n_analytic_set = 25529 , date = "", pmid="17984119", country="Norway", region="Europe" ) %>%
  add_row(study = "Oslo2", n_analytic_set = 6919 , date = "", pmid="17984119", country="Norway", region="Europe" ) %>%
  add_row(study = "HUBRO", n_analytic_set = 21361 , date = "", pmid="17984119", country="Norway", region="Europe" ) %>%
  add_row(study = "OPPHED", n_analytic_set = 12263 , date = "", pmid="17984119", country="Norway", region="Europe" ) %>%
  add_row(study = "Tromso5", n_analytic_set = 7897 , date = "", pmid="17984119", country="Norway", region="Europe" ) %>%
  add_row(study = "IHUBRO", n_analytic_set = 3614 , date = "", pmid="17984119", country="Norway", region="Europe" ) %>%
  add_row(study = "TROFINN", n_analytic_set = 9032 , date = "", pmid="17984119", country="Norway", region="Europe" ) %>%
  add_row(study = "MoRo2", n_analytic_set = 1989 , date = "", pmid="17984119", country="Norway", region="Europe" ) %>%
  
  # 4) ACC (selected non-Japan components)
  add_row(study = "Mumbai", n_analytic_set = 87397, date ="" , pmid="26243736", country="India", region="South Asia" ) %>%
  add_row(study = "Trivandrum", n_analytic_set = 49598, date = "", pmid="26243736", country="India", region="South Asia" ) %>%
  add_row(study = "CHEFS", n_analytic_set = 75379, date = "", pmid="26243736", country="China", region="East Asia" ) %>%
  add_row(study = "CBCSP", n_analytic_set = 11943, date = "", pmid="26243736", country="China", region="East Asia" ) %>%
  add_row(study = "CVDFACTS", n_analytic_set = 2262, date = "", pmid='26243736', country="China", region="East Asia" ) %>%
  add_row(study = "Seoul", n_analytic_set = 13768, date = "", pmid="26243736", country="South Korea", region="East Asia" ) %>%
  
  # 5) ANZDCC (Australia)
  add_row(study = "ANBP2", n_analytic_set = 5473, date = "", pmid='25810218', country="Australia", region="Australia" ) %>%
  add_row(study = "ALSA", n_analytic_set = 603, date = "", pmid="25810218", country="Australia", region="Australia" ) %>%
  add_row(study = "ADOLS", n_analytic_set = 10224, date = "", pmid="25810218", country="Australia", region="Australia" ) %>%
  add_row(study = "CUS", n_analytic_set = 1350, date = "", pmid="25810218", country="Australia", region="Australia" ) %>%
  add_row(study = "FDS", n_analytic_set = 1258, date = "", pmid="25810218", country="Australia", region="Australia" ) %>%
  add_row(study = "GOS", n_analytic_set = 2881, date = "", pmid="25810218", country="Australia", region="Australia" ) %>%
  add_row(study = "HIM", n_analytic_set = 10482, date = "", pmid='25810218', country="Australia", region="Australia" ) %>%
  add_row(study = "MCCS", n_analytic_set = 39893, date = "", pmid="25810218", country="Australia", region="Australia" ) %>%
  add_row(study = "NWAHS", n_analytic_set = 3681, date = "", pmid="25810218", country="Australia", region="Australia" ) %>%
  add_row(study = "PRFPCS", n_analytic_set = 3613, date = "", pmid="25810218", country="Australia", region="Australia") %>%
  
  # 6) KRIS, MIHDPS (Lithuania)
  add_row(study = "KRIS", n_analytic_set = 2447, date = "", pmid="", country="Lithuania", region="Europe" ) %>%
  add_row(study = "MIHDPS", n_analytic_set = 5933, date = "", pmid="", country="Lithuania", region="Europe") %>%
  
  # --- Standardize cohort names and remove pooled labels ----------------------

  mutate(study= case_when(study == "JPHC-II" ~ "JPHC", 
                          study == "JPHC-I" ~ "JPHC" ,
                          study == "NHS-II" ~ "NHS" , 
                          study == "NHS-I" ~ "NHS",
                          TRUE ~ study )) %>% 
  filter(!(study %in% pooled_cohorts)) %>%  # drop pooled rows
  mutate(country = case_when(country == "S. Korea" ~ "South Korea", TRUE ~ country)) %>%
  
  # --- Collapse duplicates / multi-country entries ----------------------------
  
  group_by(study) %>%
  mutate(n_analytic_set = max(n_analytic_set, na.rm=TRUE)) %>%  # keep max N
  ungroup() %>% 
  distinct(study, .keep_all = TRUE) %>%

  # --- Add EPIC sub-country splits --------------------------------------------

  add_row(study = "EPIC Greece", n_analytic_set = 28572, country="Greece", region="Europe" ) %>%
  add_row(study = "EPIC Spain", n_analytic_set = 41440, country="Spain", region="Europe" ) %>%
  add_row(study = "EPIC Italy", n_analytic_set = 47749, country="Italy", region="Europe" ) %>%
  add_row(study = "EPIC France", n_analytic_set = 72996, country="France", region="Europe" ) %>%
  add_row(study = "EPIC Germany", n_analytic_set = 53094, country="Germany", region="Europe" ) %>%
  add_row(study = "EPIC Netherlands", n_analytic_set = 40072, country="Netherlands", region="Europe" ) %>%
  add_row(study = "EPIC UK", n_analytic_set = 87940, country="UK", region="Europe" ) %>%
  add_row(study = "EPIC Denmark", n_analytic_set = 57054, country="Denmark", region="Europe" ) %>%
  add_row(study = "EPIC Sweden", n_analytic_set = 53830, country="Sweden", region="Europe" ) %>%
  add_row(study = "EPIC Norway", n_analytic_set = 37231, country="Norway", region="Europe" ) %>%
  add_row(study = "Adventist Health Study 2", n_analytic_set = 4059, country="Canada", region = "North America" ) %>%
  add_row(study = "Adventist Health Study 2", n_analytic_set = 86097, country="USA", region = "North America" )

# ---- World basemap -----------------------------------------------------------
world_map <- map_data("world") %>% 
  # Drop wrap-around polygons & Antarctica for cleaner cartogram
  filter(! long > 180) %>% filter(region != "Antarctica") 

# Country list & rough centroids (for placing labels)
countries <- world_map %>% 
  distinct(region) %>% 
  rowid_to_column() 

cnames <- aggregate(cbind(long, lat) ~ region, data=world_map, 
                    FUN=function(x)mean(range(x)))


# ---- Aggregate per-country totals & merge to map -----------------------------
country_count <- cohort_sample_size %>%
  group_by(country) %>%
  summarise(`Study population` = sum(n_analytic_set, na.rm = TRUE),
            study_count = n_distinct(study)) %>%
  ungroup() %>%
  right_join(countries, by=c("country" = "region")) %>%    # keep map countries
  right_join(cnames, by=c("country" = "region")) %>% 
  # Adjust some label anchor positions so they sit more centrally
  mutate(long = case_when(country == "USA" ~ -105,
                          country == "Canada" ~ -115,
                          country == "Australia" ~ 135,
                          TRUE ~ long)) %>%
  mutate(lat = case_when(country == "USA" ~ 40,
                         country == "Canada" ~ 57,
                         country == "Australia" ~ -25,
                         TRUE ~ lat)) 

# ---- Plot: choropleth + study-count labels -----------------------------------
country_count %>% 
  ggplot(aes(fill = `Study population`, map_id = country)) +
  geom_map(map = world_map) +
  expand_limits(x = world_map$long, y = world_map$lat) +
  scale_fill_distiller(palette = "Reds", direction = 1, na.value = "gray90", 
                       labels = scales::comma, 
                       breaks = scales::pretty_breaks(n = 5)) +
  coord_map("moll") +
  theme_map()  + 
  geom_label(data=country_count , aes(long, lat, label = study_count), fill="white", alpha=0.5,  label.size = NA, size = 6)

ggsave("global_figure.png", width = 20, height = 15, device='png',bg="white")


