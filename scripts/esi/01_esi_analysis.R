# ============================================================
# CHANGE - UC5 Malaria
#
# ESI analysis
#
# This script reads the annual MAP-ESI table and produces the
# final ESI comparison plot on screen and the concordance tables
#
# ============================================================


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

  library(dplyr)
  library(readr)
  library(ggplot2)
  library(irr)


# ------------------------------------------------------------
# 1. Data file Path
# ------------------------------------------------------------

input_file <- file.path(
  "data",
  "esi",
  "esi_map.rds"
)


# ------------------------------------------------------------
# 2. Parameters
# ------------------------------------------------------------

year_min <- 2000L
year_max <- 2024L
reference_start <- 2000L
reference_end <- 2008L
analysis_start <- 2009L
analysis_end <- 2024L


# ------------------------------------------------------------
# 3. Read annual ESI data
# ------------------------------------------------------------

comparison_data <- readRDS(
  input_file
)

# ------------------------------------------------------------
# 4. Reference range and ESI2 status
# ------------------------------------------------------------
reference_data<-comparison_data|>
  group_by(country)|>
  summarise(
         map_q50_centre=mean(map_q50),
         map_q50_scale=sd(map_q50),
         esi_q50_centre=mean(esi_q50),
         esi_q50_scale=sd(esi_q50),
  )
comparison_data<-comparison_data|>inner_join(reference_data, by=c("country"))

reference_ranges <- comparison_data |>
  filter(
    year >= reference_start,
    year <= reference_end
  ) |>
  group_by(country) |>
  summarise(
    map_reference_min = min(map_q50, na.rm = TRUE),
    map_reference_max = max(map_q50, na.rm = TRUE),
    esi_reference_min = min(esi_q50, na.rm = TRUE),
    esi_reference_max = max(esi_q50, na.rm = TRUE),
    n_reference_years = n_distinct(year),
    map_q50_centre = first(map_q50_centre),
    map_q50_scale = first(map_q50_scale),
    esi_q50_centre = first(esi_q50_centre),
    esi_q50_scale = first(esi_q50_scale),
    .groups = "drop"
  ) |>
  mutate(
    map_reference_min_z =
      (map_reference_min - map_q50_centre) / map_q50_scale,
    map_reference_max_z =
      (map_reference_max - map_q50_centre) / map_q50_scale,
    esi_reference_min_z =
      (esi_reference_min - esi_q50_centre) / esi_q50_scale,
    esi_reference_max_z =
      (esi_reference_max - esi_q50_centre) / esi_q50_scale,
    country = as.factor(country)
  )


classified_q50 <- comparison_data |>
  filter(
    year >= analysis_start,
    year <= analysis_end
  ) |>
  left_join(
    reference_ranges |>
      select(
        country,
        map_reference_min,
        map_reference_max,
        esi_reference_min,
        esi_reference_max
      ),
    by = "country"
  ) |>
  mutate(
    map_within_range =
      map_q50 >= map_reference_min &
      map_q50 <= map_reference_max,
    esi_within_range =
      esi_q50 >= esi_reference_min &
      esi_q50 <= esi_reference_max
  )

q50_points <- bind_rows(
  classified_q50 |>
    transmute(
      country,
      year,
      series = "MAP incidence Q50",
      q50_z = map_q50_z,
      range_status = if_else(
        map_within_range,
        "Within range",
        "Outside range"
      )
    ),
  classified_q50 |>
    transmute(
      country,
      year,
      series = "ESI Q50",
      q50_z = esi_q50_z,
      range_status = if_else(
        esi_within_range,
        "Within range",
        "Outside range"
      )
    )
) |>
  mutate(
    country = as.factor(country),
    range_status = factor(
      range_status,
      levels = c("Within range", "Outside range")
    )
  )


# ------------------------------------------------------------
# 5. Figure: annual MAP incidence vs ESI, 2000-2024
# ------------------------------------------------------------

line_colours <- c(
  "MAP Pf incidence Q50" = "#2166AC",
  "ESI spatial Q50" = "#D6604D"
)

line_types <- c(
  "MAP Pf incidence Q50" = "solid",
  "ESI spatial Q50" = "solid"
)

ribbon_fills <- c(
  "MAP spatial Q05-Q95" = "#2166AC",
  "ESI spatial Q05-Q95" = "#D6604D"
)

point_fills <- c(
  "Within range" = "#70AD47",
  "Outside range" = "#C83E3A"
)

ggplot(
  comparison_data,
  aes(x = year)
) +
  geom_ribbon(
    aes(
      ymin = map_q05_z,
      ymax = map_q95_z,
      fill = "MAP spatial Q05-Q95"
    ),
    alpha = 0.18,
    colour = NA
  ) +
  geom_ribbon(
    aes(
      ymin = esi_q05_z,
      ymax = esi_q95_z,
      fill = "ESI spatial Q05-Q95"
    ),
    alpha = 0.18,
    colour = NA
  ) +
  geom_hline(
    data = reference_ranges,
    aes(yintercept = map_reference_min_z),
    colour = "#2166AC",
    linetype = "dotted",
    linewidth = 0.55,
    alpha = 0.85,
    inherit.aes = FALSE
  ) +
  geom_hline(
    data = reference_ranges,
    aes(yintercept = map_reference_max_z),
    colour = "#2166AC",
    linetype = "dotted",
    linewidth = 0.55,
    alpha = 0.85,
    inherit.aes = FALSE
  ) +
  geom_hline(
    data = reference_ranges,
    aes(yintercept = esi_reference_min_z),
    colour = "#D6604D",
    linetype = "dotted",
    linewidth = 0.55,
    alpha = 0.85,
    inherit.aes = FALSE
  ) +
  geom_hline(
    data = reference_ranges,
    aes(yintercept = esi_reference_max_z),
    colour = "#D6604D",
    linetype = "dotted",
    linewidth = 0.55,
    alpha = 0.85,
    inherit.aes = FALSE
  ) +
  geom_line(
    aes(
      y = map_q50_z,
      colour = "MAP Pf incidence Q50",
      linetype = "MAP Pf incidence Q50"
    ),
    linewidth = 0.95
  ) +
  geom_line(
    aes(
      y = esi_q50_z,
      colour = "ESI spatial Q50",
      linetype = "ESI spatial Q50"
    ),
    linewidth = 0.95
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.35,
    colour = "grey45"
  ) +
  geom_point(
    data = q50_points,
    aes(
      x = year,
      y = q50_z,
      fill = range_status
    ),
    shape = 21,
    size = 2.8,
    stroke = 0.45,
    colour = "white",
    inherit.aes = FALSE
  ) +
  facet_wrap(
    ~country,
    ncol = 2,
    scales = "free_y"
  ) +
  scale_x_continuous(
    limits = c(year_min, year_max),
    breaks = seq(year_min, year_max, 4),
    minor_breaks = seq(year_min, year_max, 1),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_colour_manual(
    values = line_colours,
    name = NULL
  ) +
  scale_linetype_manual(
    values = line_types,
    name = NULL
  ) +
  scale_fill_manual(
    values = c(ribbon_fills, point_fills),
    name = NULL
  ) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = "Standardised annual value",
    caption = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(
      colour = "grey90",
      linewidth = 0.35
    ),
    strip.background = element_rect(
      fill = "grey95",
      colour = "grey75",
      linewidth = 0.4
    ),
    strip.text = element_text(
      face = "bold",
      size = 10.5
    ),
    plot.title = element_text(
      face = "bold",
      size = 14,
      margin = margin(b = 5)
    ),
    plot.subtitle = element_text(
      size = 9.2,
      lineheight = 1.1,
      margin = margin(b = 9)
    ),
    plot.caption = element_text(
      hjust = 0,
      size = 8.2,
      colour = "grey35",
      lineheight = 1.1,
      margin = margin(t = 8)
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 10
    ),
    axis.text = element_text(colour = "grey25"),
    legend.position = "none"
  )


# ------------------------------------------------------------
# 6. Concordance table and Cohen's K
# ------------------------------------------------------------

# ------------------------------------------------------------
# 6.1. Function:
# calculate Cohen's kappa and observed agreement
# ------------------------------------------------------------

calculate_kappa <- function(
    
  classified_data
  
) {
  
  if (nrow(classified_data) == 0) {
    
    return(
      
      tibble(
        
        N = 0L,
        
        agreement_n = 0L,
        
        observed_agreement = NA_real_,
        
        observed_agreement_percent = NA_real_,
        
        kappa = NA_real_,
        
        kappa_p_value = NA_real_
        
      )
      
    )
    
  }
  
  pairs <- classified_data |>
    
    transmute(
      
      ESI = as.factor(
        
        if_else(esi_within_range,"IN","OUT")

      ),
      
      Incidence = as.factor(
        
        if_else(map_within_range,"IN","OUT")
        
      )
      
    )
  
  
  # Observed agreement
  
  agreement_n <- sum(
    
    pairs$ESI ==
      
      pairs$Incidence
    
  )
  
  
  observed_agreement <-
    
    agreement_n /
    
    nrow(pairs)
  
  
  # Cohen's kappa
  
  kappa_result <- tryCatch(
    
    irr::kappa2(
      
      pairs,
      
      weight = "unweighted"
      
    ),
    
    error = function(e) {
      
      NULL
      
    }
    
  )
  
  
  if (is.null(kappa_result)) {
    
    kappa_value <-
      
      NA_real_
    
    kappa_p_value <-
      
      NA_real_
    
  } else {
    
    kappa_value <-
      
      as.numeric(
        
        kappa_result$value
        
      )
    
    kappa_p_value <-
      
      as.numeric(
        
        kappa_result$p.value
        
      )
    
  }
  
  
  tibble(
    
    N =
      
      nrow(pairs),
    
    agreement_n =
      
      agreement_n,
    
    observed_agreement =
      
      observed_agreement,
    
    observed_agreement_percent =
      
      100 *
      
      observed_agreement,
    
    kappa =
      
      kappa_value,
    
    kappa_p_value =
      
      kappa_p_value
    
  )
  
}


sapply(
  unique(classified_q50$country),
  function(x)
    {calculate_kappa(classified_q50%>%
                       filter(country==x)
                     )
  }
)
