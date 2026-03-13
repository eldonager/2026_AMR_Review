##Pooled AMR prevalence by antibiotic compound and species per region

library(tidyverse)
library(ggpubr)
library(patchwork)
library(grid)
library(gridExtra)
library(ggthemes)
library(flextable)
library(officer)

theme_Publication <- function(base_size = 17, base_family = "helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size = base_size, base_family = base_family)
    + theme(
      plot.title       = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
      text             = element_text(),
      panel.background = element_rect(colour = NA),
      plot.background  = element_rect(colour = NA),
      panel.border     = element_rect(colour = NA),
      axis.title       = element_text(face = "bold", size = rel(1)),
      axis.title.y     = element_text(angle = 90, vjust = 2),
      axis.title.x     = element_text(vjust = -0.2),
      axis.text        = element_text(),
      axis.line        = element_line(colour = "black"),
      axis.ticks       = element_line(),
      panel.grid.major = element_line(colour = "#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key       = element_rect(colour = NA),
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.key.size  = unit(0.2, "cm"),
      legend.margin    = margin(0, 0, 0, 0),
      legend.title     = element_text(face = "italic"),
      plot.margin      = unit(c(2, 2, 2, 2), "mm"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text       = element_text(face = "bold")
    ))
}

# Consistent species colour palette across all panels
species_colours <- c(
  "Cattle"      = "#386cb0",
  "Chicken"     = "#fdb462",
  "Environment" = "#7fc97f",
  "Goats"       = "#ef3b2c",
  "Pigs"        = "#662506",
  "Sheep"       = "#a6cee3",
  "Turkey"      = "#984ea3"
)

all_species <- names(species_colours)

###############################################
# Load and clean data
###############################################
AMR_clean <- read.csv("AMR_clean.csv") %>%
  dplyr::select(doi, country, region, iso_3, y_coordinate, x_coordinate,
                sampling_start_year, sampling_end_year, species, pathogen,
                sal_prevalence, no_sample, antimicrobial, antimicrobial_compound,
                antibiotic_class, who_classification, no_isolate,
                no_isolates_resistant, no_isolates_intermediate,
                no_isolates_susceptible, mdr_percentage) %>%
  dplyr::mutate(
    no_isolates_susceptible = no_isolate - no_isolates_resistant,
    species = ifelse(species == "Pig", "Pigs", species)
  ) %>%
  dplyr::select(-no_isolates_intermediate)

# Total isolates per species per region - correct denominator
# Summed before antibiotic aggregation to avoid double counting
isolate_counts <- AMR_clean %>%
  group_by(region, species) %>%
  summarise(NIsolates = sum(no_isolate, na.rm = TRUE), .groups = "drop")

###############################################
# Subtitle functions - dynamic, no hardcoding
###############################################

# Eastern Africa - line break before Environment
build_subtitle_ea <- function() {
  d <- isolate_counts %>%
    filter(region == "Eastern Africa") %>%
    arrange(species) %>%
    mutate(label = paste0(species, " n = ", format(NIsolates, big.mark = ",")))
  line1_labels <- d %>% filter(species < "Environment") %>% pull(label)
  line2_labels <- d %>% filter(species >= "Environment" & species < "Sheep") %>% pull(label)
  line3_labels <- d %>% filter(species >= "Sheep") %>% pull(label)
  line1 <- paste(c("Eastern Africa", line1_labels), collapse = ", ")
  line2 <- paste(line2_labels, collapse = ", ")
  line3 <- paste(line3_labels, collapse = ", ")
  paste0(line1, ",\n", line2, ",\n", line3)
}

# West Africa - line break before Environment
build_subtitle_wa <- function() {
  d <- isolate_counts %>%
    filter(region == "West Africa") %>%
    arrange(species) %>%
    mutate(label = paste0(species, " n = ", format(NIsolates, big.mark = ",")))
  before_env <- d %>% filter(species < "Environment") %>% pull(label)
  env_onward  <- d %>% filter(species >= "Environment") %>% pull(label)
  line1 <- paste(c("West Africa", before_env), collapse = ", ")
  line2 <- paste(env_onward, collapse = ", ")
  paste0(line1, ",\n", line2)
}

# Northern Africa - exclude Environment, line break before Turkey
build_subtitle_na <- function() {
  d <- isolate_counts %>%
    filter(region == "Northern Africa", species != "Environment") %>%
    arrange(species) %>%
    mutate(label = paste0(species, " n = ", format(NIsolates, big.mark = ",")))
  before_turkey <- d %>% filter(species < "Turkey") %>% pull(label)
  turkey_onward  <- d %>% filter(species >= "Turkey") %>% pull(label)
  line1 <- paste(c("Northern Africa", before_turkey), collapse = ", ")
  line2 <- paste(turkey_onward, collapse = ", ")
  paste0(line1, ",\n", line2)
}

# Southern Africa - line break before Goats
build_subtitle_sa <- function() {
  d <- isolate_counts %>%
    filter(region == "Southern Africa") %>%
    arrange(species) %>%
    mutate(label = paste0(species, " n = ", format(NIsolates, big.mark = ",")))
  before_goats <- d %>% filter(species < "Goats") %>% pull(label)
  goats_onward  <- d %>% filter(species >= "Goats") %>% pull(label)
  line1 <- paste(c("Southern Africa", before_goats), collapse = ", ")
  line2 <- paste(goats_onward, collapse = ", ")
  paste0(line1, ",\n", line2)
}

###############################################
# Aggregate and compute pooled prevalence
###############################################
SalRes <- aggregate(no_isolates_resistant ~ antimicrobial_compound + species + region,
                    data = AMR_clean, FUN = sum)
SalAll <- aggregate(no_isolate ~ antimicrobial_compound + species + region,
                    data = AMR_clean, FUN = sum)

SalMean <- round((SalRes$no_isolates_resistant / SalAll$no_isolate) * 100, digits = 0)

SalMeandf <- as.data.frame(cbind(SalRes$antimicrobial_compound, SalRes$species,
                                 SalRes$region, SalMean, SalAll$no_isolate))
colnames(SalMeandf) <- c("Compound", "Species", "Region", "Mean", "NIsolates")
SalMeandf$Mean      <- as.numeric(as.character(SalMeandf$Mean))
SalMeandf$NIsolates <- as.numeric(as.character(SalMeandf$NIsolates))

# Restrict to drug-species pairings with >= 10 isolates
SalMeandf <- SalMeandf[SalMeandf$NIsolates >= 10, ]

hist(SalMeandf$Mean, main = "Distribution of Mean Prevalence", xlab = "Mean Prevalence")

# 95% Wald CI: p +/- 1.96 * sqrt(p*(1-p)/n)
CI.function <- function(x, y) { x + c(-1.96, 1.96) * sqrt(x * (1 - x) / y) }

CIsal     <- as.data.frame(t(mapply(CI.function, SalMeandf$Mean / 100, SalMeandf$NIsolates)))
SalMeandf <- cbind(SalMeandf, round(CIsal * 100, digits = 0))
colnames(SalMeandf) <- c("Compound", "Species", "Region", "Mean", "NIsolates", "CILow", "CIHigh")

SalMeandf$CILow[SalMeandf$CILow   < 0]   <- 0
SalMeandf$CIHigh[SalMeandf$CIHigh > 100] <- 100

# Remove zero-resistance rows
SalMeandf <- SalMeandf[SalMeandf$Mean != 0, ]

###############################################
# Supplementary Tables
###############################################
zebra_header <- "#2E75B6"
zebra_odd    <- "#DEEAF1"
zebra_even   <- "#FFFFFF"

make_zebra_table <- function(df, title) {
  n_rows <- nrow(df)
  flextable(df) %>%
    bg(bg = zebra_header, part = "header") %>%
    color(color = "white", part = "header") %>%
    bold(part = "header") %>%
    bg(i = seq(1, n_rows, by = 2), bg = zebra_odd,  part = "body") %>%
    bg(i = seq(2, n_rows, by = 2), bg = zebra_even, part = "body") %>%
    border_outer(border = fp_border(color = "#2E75B6", width = 1.5)) %>%
    border_inner_h(border = fp_border(color = "#BDD7EE", width = 0.5)) %>%
    font(fontname = "Arial", part = "all") %>%
    fontsize(size = 10, part = "body") %>%
    fontsize(size = 11, part = "header") %>%
    align(align = "center", part = "header") %>%
    align(align = "left",   part = "body") %>%
    set_table_properties(width = 1, layout = "autofit") %>%
    set_caption(caption = title)
}

regions <- unique(SalMeandf$Region)
reg.ca  <- SalMeandf %>% filter(Region == regions[1])
reg.ea  <- SalMeandf %>% filter(Region == regions[2])
reg.na  <- SalMeandf %>% filter(Region == regions[3])
reg.sa  <- SalMeandf %>% filter(Region == regions[4])
reg.wa  <- SalMeandf %>% filter(Region == regions[5])

ft.ca <- make_zebra_table(reg.ca, paste("Table. Pooled AMR prevalence -", regions[1]))
ft.ea <- make_zebra_table(reg.ea, paste("Table. Pooled AMR prevalence -", regions[2]))
ft.na <- make_zebra_table(reg.na, paste("Table. Pooled AMR prevalence -", regions[3]))
ft.sa <- make_zebra_table(reg.sa, paste("Table. Pooled AMR prevalence -", regions[4]))
ft.wa <- make_zebra_table(reg.wa, paste("Table. Pooled AMR prevalence -", regions[5]))

doc <- read_docx() %>%
  body_add_par("Supplementary Tables", style = "heading 1") %>%
  body_add_flextable(ft.ca) %>% body_add_par("") %>%
  body_add_flextable(ft.ea) %>% body_add_par("") %>%
  body_add_flextable(ft.na) %>% body_add_par("") %>%
  body_add_flextable(ft.sa) %>% body_add_par("") %>%
  body_add_flextable(ft.wa)

print(doc, target = "AMRprevalencespeciesregion_tables.docx")

###############################################
# Shared legend - all species
###############################################
dummy_df <- data.frame(Species = all_species, x = 1, y = 1)

dummy_plot <- ggplot(dummy_df, aes(x = x, y = y, fill = Species)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = species_colours) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
  theme(legend.position  = "right",
        legend.direction = "vertical",
        legend.title     = element_text(size = 17, face = "italic"),
        legend.text      = element_text(size = 17),
        legend.key.size  = unit(0.5, "cm"),
        legend.spacing.y = unit(0.4, "cm")) +
  labs(fill = "Species")

get_legend <- function(p) {
  g   <- ggplotGrob(p)
  leg <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
  return(leg)
}

shared_legend <- get_legend(dummy_plot)

###############################################
# Eastern Africa
###############################################
ea_compounds <- c("ERY", "STR", "TET", "AMP", "NIT", "KAN", "NAL", "SMX", "AMC")

sal.ea <- SalMeandf %>%
  filter(Compound %in% ea_compounds, Region == "Eastern Africa")

sal.ea$Compound <- reorder(sal.ea$Compound, -sal.ea$Mean)

sal.ea.plot <- ggbarplot(sal.ea, "Compound", "Mean",
                         fill = "Species",
                         position = position_dodge(width = 0.7),
                         subtitle = build_subtitle_ea(),
                         xlab = FALSE, ylab = FALSE,
                         legend.title = "Species", font.subtitle = c(10)) +
  rotate_x_text(90) +
  geom_linerange(aes(group = Species, ymax = CIHigh, ymin = CILow),
                 position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = species_colours, drop = TRUE) +
  theme_Publication() +
  theme(axis.title.x   = element_blank(),
        axis.title.y   = element_blank(),
        axis.text      = element_text(size = 17),
        legend.position = "none")

sal.ea.plot <- ggpar(sal.ea.plot, ylim = c(0, 100))

###############################################
# West Africa
###############################################
wa_compounds <- c("ERY", "AMX", "AMP", "AMC", "COT", "CTX", "TET", "OXA", "IMP")

sal.wa <- SalMeandf %>%
  filter(Compound %in% wa_compounds, Region == "West Africa")

sal.wa$Compound <- reorder(sal.wa$Compound, -sal.wa$Mean)

sal.wa.plot <- ggbarplot(sal.wa, "Compound", "Mean",
                         fill = "Species",
                         position = position_dodge(width = 0.7),
                         subtitle = build_subtitle_wa(),
                         xlab = FALSE, ylab = FALSE,
                         legend.title = "Species", font.subtitle = c(10)) +
  rotate_x_text(90) +
  geom_linerange(aes(group = Species, ymax = CIHigh, ymin = CILow),
                 position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = species_colours, drop = TRUE) +
  theme_Publication() +
  theme(axis.title.x   = element_blank(),
        axis.title.y   = element_blank(),
        axis.text      = element_text(size = 17),
        legend.position = "none")

sal.wa.plot <- ggpar(sal.wa.plot, ylim = c(0, 100))

###############################################
# Northern Africa
###############################################
na_compounds <- c("AMX", "AMP", "AMC", "COT", "CTX", "STR", "TET", "DOX", "KAN")

sal.na <- SalMeandf %>%
  filter(Compound %in% na_compounds, Region == "Northern Africa")

sal.na$Compound <- reorder(sal.na$Compound, -sal.na$Mean)

sal.na.plot <- ggbarplot(sal.na, "Compound", "Mean",
                         fill = "Species",
                         position = position_dodge(width = 0.7),
                         subtitle = build_subtitle_na(),
                         xlab = FALSE, ylab = FALSE,
                         legend.title = "Species", font.subtitle = c(10)) +
  rotate_x_text(90) +
  geom_linerange(aes(group = Species, ymax = CIHigh, ymin = CILow),
                 position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = species_colours, drop = TRUE) +
  theme_Publication() +
  theme(axis.title.x   = element_blank(),
        axis.title.y   = element_blank(),
        axis.text      = element_text(size = 17),
        legend.position = "none")

sal.na.plot <- ggpar(sal.na.plot, ylim = c(0, 100))

###############################################
# Southern Africa
###############################################
sa_compounds <- c("ERY", "AMX", "TET", "STR", "CTX", "CFX", "AMP", "COT", "TMP")

sal.sa <- SalMeandf %>%
  filter(Compound %in% sa_compounds, Region == "Southern Africa")

sal.sa$Compound <- reorder(sal.sa$Compound, -sal.sa$Mean)

sal.sa.plot <- ggbarplot(sal.sa, "Compound", "Mean",
                         fill = "Species",
                         position = position_dodge(width = 0.7),
                         subtitle = build_subtitle_sa(),
                         xlab = FALSE, ylab = FALSE,
                         legend.title = "Species", font.subtitle = c(10)) +
  rotate_x_text(90) +
  geom_linerange(aes(group = Species, ymax = CIHigh, ymin = CILow),
                 position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = species_colours, drop = TRUE) +
  theme_Publication() +
  theme(axis.title.x   = element_blank(),
        axis.title.y   = element_blank(),
        axis.text      = element_text(size = 17),
        legend.position = "none")

sal.sa.plot <- ggpar(sal.sa.plot, ylim = c(0, 100))


###############################################
# Combine and save Figure 6 - legend on right
###############################################
region_plot <- (sal.ea.plot | sal.wa.plot) / (sal.na.plot | sal.sa.plot) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(plot.tag = element_text(face = "bold"))

yleft  <- textGrob("Percent resistance", rot = 90, gp = gpar(fontsize = 19))
bottom <- textGrob("Antimicrobial",               gp = gpar(fontsize = 19))

plot1 <- patchwork::patchworkGrob(region_plot)

# Legend sized to match figure height
legend_plot <- ggplot(dummy_df, aes(x = x, y = y, fill = Species)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = species_colours) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
  theme_void() +
  theme(legend.position  = "right",
        legend.direction = "vertical",
        legend.title     = element_text(size = 17, face = "italic", vjust = 1),
        legend.text      = element_text(size = 17),
        legend.key.size  = unit(0.6, "cm"),
        legend.spacing.y = unit(0.4, "cm")) +
  labs(fill = "Species")

shared_legend <- get_legend(legend_plot)

final_plot <- gridExtra::arrangeGrob(
  plot1,
  shared_legend,
  ncol   = 2,
  nrow   = 1,
  widths = c(20, 3)
)

tiff("figure6.tiff", width = 5000, height = 4200, res = 300)
grid.arrange(final_plot, left = yleft, bottom = bottom)
dev.off()

###############################################
# Legend figure - standalone for separate use
###############################################
tiff("figure6_legend.tiff", width = 600, height = 800, res = 300)
grid.arrange(shared_legend)
dev.off()