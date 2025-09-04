# =============================================================================
# Figure 2: Study size and regional totals over time
# Purpose: Clean cohort metadata, lump rare regions, and visualize per-study
#          cancer case counts (log scale) by publication year with a smoothed trend.
# =============================================================================

# ---- Packages ----------------------------------------------------------------
pacman::p_load("tidyverse", "data.table", "ggplot2", "readxl", "janitor", 
               "hrbrthemes",  "ggrepel", "ggthemes")

# ---- Load & clean data -------------------------------------------------------
cancer_allsites <- read_excel("SupplementaryData.xlsx", 
                              sheet = "BMI_cancer_data") 

# Move first row to column names
colnames(cancer_allsites) <- as.character(unlist(cancer_allsites[1, ]))
cancer_allsites <- cancer_allsites[-1, ] # Remove the old first row


cancer_allsites <- cancer_allsites %>%
  clean_names() %>%
  mutate(n_analytic_set = str_replace_all(n_analytic_set, "\\s", ""),
         n_analytic_set = as.numeric(n_analytic_set),
         n_cases = as.numeric(n_cases),
         year = as.numeric(year))

# Remove parenthetical qualifiers from study names, e.g., "KNHIS (women)" -> "KNHIS"
cancer_allsites$study <- sub("\\s*\\(.*", "", cancer_allsites$study)

# ---- Collapse to one row per study within each cancer site -------------------
# Sum duplicated entries (same site+study), drop global pooled rows, remove missing n_cases,
# and lump small/rare regions into "Other" for a cleaner legend.

all_studies <- cancer_allsites %>%
  group_by(top_25_cancer, study) %>%
  mutate(n_analytic_set = sum(n_analytic_set),
         n_cases = sum(n_cases)) %>% 
  ungroup() %>%
  distinct(top_25_cancer, study, .keep_all = TRUE) %>% 
  # # drop global pooled analyses
  filter(region !="Global") %>%
  filter(!is.na(n_cases)) %>%
  mutate(region = case_when(region == "Middle East" |
                              region == "South Asia" |
                              region == "Australia" ~ "Other",
                            TRUE ~ region))

# Select columns for plotting
plot_dat <- all_studies %>% 
  select(year, n_cases, region, n_analytic_set)

# ---- Constraint point to stabilize early years --------------------------------
# Add a tiny "anchor" point at year 1987 to keep the smoother from drifting.
constraint_point <- data.frame(year = 1987, 
                               n_cases = 1,
                               region = "Europe", 
                               n_analytic_set = 1)

plot_dat <- rbind(plot_dat, constraint_point)


# ---- Palette (colorblind-friendly) -------------------------------------------

cbp1 <- c( "#E69F00", "#56B4E9", "#009E73","#CC79A7",
          "#F0E442", "#0072B2", "#D55E00", "#999999")

# ---- Base plot factory --------------------------------------------------------

create_base_plot <- function(data) {
  ggplot(data, aes(x=year, y=n_cases)) +
    geom_jitter(alpha=0.5, shape=21, color="black", width = 1, height = 0, 
                aes(fill=region, size=n_analytic_set)) +
    scale_size(range = c(1, 40), name="Study size", labels = scales::comma) +
    theme_classic() +
    scale_x_continuous(limits =c(1987, 2026), expand=c(0,0)) +
    ylab("N cancer cases (log scale)") +
    xlab("Publication year") +
    scale_fill_manual(values = cbp1) +
    guides(fill = guide_legend(override.aes = list(size = 20), title = "Region")) +
    coord_cartesian(xlim = c(1987, 2027), ylim= c(1, 350000)) +
    scale_y_continuous(limits =c(1, 300000), expand=c(0,0),
                       trans = 'log10',
                       breaks = scales::trans_breaks('log10', function(x) 10^x),
                       labels = scales::comma) +
    theme(axis.text=element_text(size=27), 
          axis.title.x = element_text(size = 30, hjust = 0.5), 
          axis.title.y = element_text(size = 30, hjust = 0.5),
          legend.text = element_text(size = 25), 
          legend.title = element_text(size = 25))
}

# ---- Create figure ------------------------------------------------------------

p3 <- create_base_plot(plot_dat) +
  stat_smooth(method = lm, 
              formula = y ~ splines::ns(x, 3),
              fullrange = TRUE, 
              colour = "darkgrey", 
              se = TRUE, 
              linewidth = 1) 

p3

ggsave("Fig3.pdf",
       plot = p3, width = 16, height = 11, 
       device = cairo_pdf, bg = "white")

