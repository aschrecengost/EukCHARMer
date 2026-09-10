#!/usr/bin/env Rscript
# Generate the taxa barplot and sample map.
#  - barplot: samples merged by habitat type (env_2), relative abundance,
#    filled by Family.
#  - map: the "map of all samples" version (points from the cleaned metadata,
#    not the phyloseq object), colored by env_2, shaped by shallow_deep.
#
# Usage:
#   Rscript generate_figures.R <ps.rds> <metadata_cleaned.csv> \
#       <out_barplot.pdf> <out_barplot.png> <out_map.pdf> <out_map.png>

suppressPackageStartupMessages({
    library(phyloseq)
    library(ggplot2)
    library(dplyr)
    library(sf)
    library(rnaturalearth)
    library(rnaturalearthdata)
    library(viridis)
})

args <- commandArgs(trailingOnly = TRUE)
ps_rds       <- args[1]
metadata_csv <- args[2]
out_bar_pdf  <- args[3]
out_bar_png  <- args[4]
out_map_pdf  <- args[5]
out_map_png  <- args[6]

ps <- readRDS(ps_rds)

# ---- taxa barplot ----

# samples without an env_2 assignment (env_biome empty/unmappable and no
# manual fill) cannot be grouped - drop them before merging
n_na <- sum(is.na(sample_data(ps)$env_2))
if (n_na > 0) {
    message("barplot: dropping ", n_na, " samples with NA env_2")
    ps <- prune_samples(!is.na(sample_data(ps)$env_2), ps)
}

# combine samples by habitat type and convert to relative abundance
# change "env_2" to any metadata variable which you want to group samples by
# add "tax_glom" if you want to agglomerate by a taxonomic rank
ps_rel_habitat <- ps %>%
    merge_samples("env_2") %>%
    transform_sample_counts(function(x) x / sum(x)) # %>% tax_glom(taxrank="Class")

# optional - order habitat types and choose colors
# habitat_names <- c()
# habitat_colors <- c()

p_bar <- plot_bar(ps_rel_habitat, fill = "Family") +
    coord_flip() +
    geom_bar(stat = "identity", color = NA, linewidth = 0) +
    scale_fill_viridis(discrete = TRUE) #+
    #scale_fill_manual(breaks=habitat_names,values=habitat_colors)

ggsave(out_bar_pdf, p_bar, width = 12, height = 6)
ggsave(out_bar_png, p_bar, width = 12, height = 6, dpi = 300)
message("barplot: ", nsamples(ps_rel_habitat), " habitat groups, ",
        ntaxa(ps_rel_habitat), " taxa")

# ---- sample map (all samples) ----

filtered_ps <- prune_samples(sample_sums(ps) >= 1, ps)
                            
dataset_coords<-data.frame(filtered_ps@sam_data)
dataset_coords$lon_converted <- as.numeric(dataset_coords$lon_converted)
dataset_coords$lat_converted <- as.numeric(dataset_coords$lat_converted)
dataset_coords <- dataset_coords[!is.na(dataset_coords$lon_converted) &
                                 !is.na(dataset_coords$lat_converted), ]
dataset_coords_sf <- st_as_sf(dataset_coords,
                              coords = c("lon_converted", "lat_converted"),
                              crs = 4326)

p_map <- ggplot() +
    geom_sf(data = ne_countries(returnclass = "sf"),
            fill = "lightgoldenrod1", color = NA) +  # Base map
    # color corresponds to habitat type, shape corresponds to depth category
    geom_sf(data = dataset_coords_sf,
            aes(fill = env_2, shape = shallow_deep),
            size = 5, show.legend = TRUE) +
    theme_bw() +
    theme(
        text = element_text(size = 20),      # Text size
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.text = element_blank(),         # Remove coordinate labels
        axis.ticks = element_blank(),        # Remove axis ticks
        panel.border = element_blank()       # Remove panel border
    ) +
    scale_fill_viridis(discrete = TRUE, option = "A") +
    scale_shape_manual(values = c(21, 23)) +
    coord_sf(expand = FALSE) +
    ylab("") +
    xlab("") +
    guides(
        fill = guide_legend(title = "Environmental Factor",
                            override.aes = list(shape = 21)),
        shape = guide_legend(title = "Depth Type")
    )

ggsave(out_map_pdf, p_map, width = 14, height = 8)
ggsave(out_map_png, p_map, width = 14, height = 8, dpi = 300)
message("map: ", nrow(dataset_coords), " samples with valid coordinates (of ",
        nrow(read.csv(metadata_csv)), ")")
