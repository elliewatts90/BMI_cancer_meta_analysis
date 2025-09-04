# =============================================================================
# Figure 3A and 3B: Percentage of cancer cases by region (per site + overall)
# - Stacked 100% bar chart (per cancer site)
# - Overall pie chart (all sites combined)
# =============================================================================


# ---- Packages ----------------------------------------------------------------

pacman::p_load("tidyverse", "data.table", "janitor", "meta", "readxl", "devEMF", "ggrepel")

# ---- Inputs ------------------------------------------------------------------
cancer_allsites <- read_excel("SupplementaryData.xlsx", 
                              sheet = "BMI_cancer_data") 

# Move first row to column names
colnames(cancer_allsites) <- as.character(unlist(cancer_allsites[1, ]))
cancer_allsites <- cancer_allsites[-1, ] # Remove the old first row

cancer_allsites <- cancer_allsites %>%  
  clean_names() %>%
  select(author:pmid) %>%
  mutate(n_cases = as.numeric(n_cases)) %>% 
  filter(region!= "Global") %>% 
  mutate(region = case_when(region != "Europe" & region != "North America" & region != "East Asia"  ~ "Other",
                            TRUE ~ region))


# ---- Site ordering ------------------------------------------------------------
cancer_list_ordered <- fread("meta_analysis_res.csv") %>% arrange(desc(OR)) %>% dplyr::select(Site) %>% unlist()

# Apply factor levels to ensure bars are ordered top→bottom by the supplied order
cancer_allsites$cancer_site <- factor(cancer_allsites$top_25_cancer, levels = cancer_list_ordered)
cancer_allsites <- cancer_allsites[order(cancer_allsites$cancer_site), ]


# ---- Color palettes -----------------------------------------------------------
cbp1 <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")



# ---- Prepare data for stacked 100% bars --------------------------------------
# Sum n_cases at (site, region); fill missing region categories with 0 so the
# proportions stack to 100% for every site (important for position = "fill").
cases_by_site_region <- cancer_allsites %>%
  group_by(cancer_site, region) %>%
  summarise(n_cases = sum(n_cases, na.rm = TRUE), .groups = "drop") %>%
  tidyr::complete(cancer_site, region, fill = list(n_cases = 0))


# ---- Plot 1: % of cases by region within each cancer site --------------------
ggplot(
  data = cases_by_site_region,
  aes(x = reorder(cancer_site, desc(cancer_site)), y = n_cases, fill = region)) + 
  geom_col(position = "fill") +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = cbp1, name = "Region") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Cancer type", y = "Percentage of cases by region") +
  theme(
    panel.border     = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black")
  )


ggsave("Fig3B.pdf",
       width = 8, height = 8, 
       device = cairo_pdf, bg = "white")

dev.off()

  
# ---- Overall proportions (all sites) -----------------------------------------
# Aggregate *total* cases across sites by region, compute proportions,

overall_proportions <- cancer_allsites %>%
  mutate(region = if_else(region %in% c("Europe", "North America", "East Asia"), region, "Other")) %>%
  group_by(region) %>%
  summarise(n_cases = sum(n_cases, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    proportion = n_cases / sum(n_cases),
    label      = if_else(proportion < 0.0001, NA_real_, proportion * 100),
    label      = if_else(is.na(label), NA_character_, paste0(sprintf("%.1f", label), "%"))
  ) %>%
  arrange(desc(region)) %>%
  # Compute label radial positions
  mutate(ypos = cumsum(proportion) - 0.5 * proportion,
         label = scales::percent(proportion, accuracy = 0.1))


# ---- Plot 2: Overall pie chart ------------------------------------------------

ggplot(overall_proportions, aes(x = "", y = proportion, fill = region)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  blank_theme +
  scale_fill_manual(values = cbp1, name = "Region") +
  theme(axis.text.x=element_blank(), 
        legend.text=element_text(size=15), 
        legend.title=element_text(size=15)) +
  ggrepel::geom_label_repel(
    data = overall_proportions,
    aes(x = 1, y = ypos, label = label),
    nudge_x = 0.7,           
    direction = "y",
    force = 5,
    min.segment.length = 0,   
    segment.size = 0.2,
    size = 5,
    show.legend = FALSE
  )

ggsave("Fig3A.pdf",
       width = 8, height = 8, 
       device = cairo_pdf, bg = "white")

dev.off()



