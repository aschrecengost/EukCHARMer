#!/usr/bin/env Rscript
# Build the three phyloseq input objects (sample_data, tax_table, otu_table)
# from pipeline outputs, plus the two combined TSVs they are derived from.
#
# Adapted from the import sections of map_and_plot_snakemake.R:
#   meta <- read.table("TableS2.csv", sep=",", header=TRUE)
#   META <- sample_data(data.frame(meta, row.names=meta$'Run'))
#   taxonomy_pp <- read.csv("unassigned/pp_taxonomy.tsv", sep="\t", row.names=1)
#   TAX_pp <- tax_table(as.matrix(taxonomy_pp))
#   otu_all <- read.csv("unassigned/anaerobic_ciliate_otu_table_full.tsv", sep="\t", row.names=1)
#   ASV_all <- otu_table(data.frame(otu_all), taxa_are_rows=TRUE)
#
# Differences from the original:
#  - sample IDs come from `run_accession` (ERR/SRR accessions, matching the
#    count-table sample columns) instead of `Run`.
#  - the taxonomy table is built here from the EDPL/LWR-filtered clade-tree
#    placements (filtered_{APM,Plagio,Scuti}_LWR_EDPL.tsv), so it covers only
#    the high-confidence anaerobic-ciliate placements.
#  - the count table is built here by combining the exported ciliophora and
#    unassigned tables. Their sample columns differ, so this is an outer
#    join on samples (missing combinations filled with 0), not a cat.
#  - env_2 is derived from env_biome (ENVO codes resolved to labels via the
#    BioPortal ENVO export); shallow_deep is derived from depth_m
#    (<= 200 m -> shallow, > 200 m -> deep, missing -> NA).
#
# Usage:
#   Rscript prepare_phyloseq_objects.R \
#     <merged_metadata.csv> \
#     <filtered_APM.tsv> <filtered_Plagio.tsv> <filtered_Scuti.tsv> \
#     <ciliophora_counts.tsv> <unassigned_counts.tsv> \
#     <out_pp_taxonomy.tsv> <out_combined_counts.tsv> <out_metadata_cleaned.csv> \
#     <out_metadata.rds> <out_tax_table.rds> <out_otu_table.rds> <out_ps.rds> \
#     <ENVO_IDs.csv> <out_missing_data_report.txt>
#
# envo_mapping is reference/ENVO_IDs.csv, a BioPortal export of ENVO
# (columns used: ID, Preferred Label).
#
# Public SRA metadata is routinely incomplete or inconsistent across projects,
# so this script does not try to hand-patch specific values for specific
# projects - see report_missing() below, which surfaces exactly what's
# missing (per project, per column) instead, so gaps can be dealt with
# upstream (in the source metadata) rather than papered over here.

suppressPackageStartupMessages({
    library(phyloseq)
})

# Parse a single coordinate token (one half of a split lat_lon pair) into a
# signed decimal-degree value, with no external package dependency (parzer
# has no conda-forge/bioconda build compatible with R 4.4). Handles every
# format observed in this pipeline's metadata plus DMS as a documented
# extension: decimal + hemisphere letter ("24.75 N"), signed decimal with no
# letter ("-69.648632"), and degrees[-minutes[-seconds]] with any separator
# ("40 26 46 N", "40d26m46sN", "40:26:46N").
parse_coord <- function(x) {
    x <- trimws(x)
    if (is.na(x) || x == "") return(NA_real_)

    has_hemi <- grepl("[NSEWnsew]\\s*$", x)
    hemi <- if (has_hemi) toupper(sub(".*?([NSEWnsew])\\s*$", "\\1", x)) else NA
    body <- trimws(if (has_hemi) sub("[NSEWnsew]\\s*$", "", x) else x)

    nums <- suppressWarnings(as.numeric(
        regmatches(body, gregexpr("-?[0-9]+\\.?[0-9]*", body))[[1]]))
    if (length(nums) == 0 || anyNA(nums)) return(NA_real_)

    deg <- nums[1]
    minute <- if (length(nums) >= 2) nums[2] else 0
    second <- if (length(nums) >= 3) nums[3] else 0
    value <- abs(deg) + minute / 60 + second / 3600
    if (deg < 0) value <- -value
    if (has_hemi && hemi %in% c("S", "W")) value <- -value
    value
}

args <- commandArgs(trailingOnly = TRUE)
metadata_csv       <- args[1]
placement_files    <- c(apm = args[2], plagio = args[3], scuti = args[4])
cil_counts_tsv     <- args[5]
una_counts_tsv     <- args[6]
out_taxonomy_tsv   <- args[7]
out_counts_tsv     <- args[8]
out_meta_clean_csv <- args[9]
out_meta_rds       <- args[10]
out_tax_rds        <- args[11]
out_otu_rds        <- args[12]
out_ps_rds         <- args[13]
envo_map_csv       <- args[14]
out_missing_report <- args[15]

# ---- import and clean metadata ----

# Columns kept for figure generation (edit this list to keep more/fewer).
# Dropped: experiment_title, organism_name, library_* (uninformative here),
# sample_accession (redundant with biosample), sample_title_x (always empty).
KEEP_COLS <- c(
    "run_accession", "bioproject", "biosample", "collection_date",
    "geo_loc_name", "lat_lon", "depth", "depth_m", "depth_layer", "elev",
    "env_biome", "env_feature", "env_material", "samp_collect_device",
    "description"
)

meta <- read.table(metadata_csv, sep = ",", header = TRUE)
if (anyDuplicated(meta$run_accession)) {
    stop("Duplicate run_accession values in ", metadata_csv)
}
meta <- meta[, intersect(KEEP_COLS, names(meta))]

# Reformat lat_lon into uniform decimal-degree columns using parse_coord()
# (see above). Observed formats: "24.75 N 122.83 E" (decimal + hemisphere
# letters) and "43.994827 -69.648632" (signed decimal); parse_coord also
# handles DMS should it ever appear. Split the combined field into lat/lon
# halves, then parse each.
split_lat_lon <- function(s) {
    s <- trimws(s)
    if (is.na(s) || s == "") return(c(NA_character_, NA_character_))
    hemi <- regmatches(s, regexec("^(.*[NSns])[ ,;]+(.*[EWew])$", s))[[1]]
    if (length(hemi) == 3) return(hemi[2:3])
    parts <- strsplit(s, "[ ,;]+")[[1]]
    if (length(parts) == 2) return(parts)
    c(NA_character_, NA_character_)
}
halves <- t(vapply(meta$lat_lon, split_lat_lon, character(2)))
meta$lat_converted <- vapply(halves[, 1], parse_coord, numeric(1), USE.NAMES = FALSE)
meta$lon_converted <- vapply(halves[, 2], parse_coord, numeric(1), USE.NAMES = FALSE)

n_fail <- sum(!is.na(meta$lat_lon) & meta$lat_lon != "" &
              (is.na(meta$lat_converted) | is.na(meta$lon_converted)))
message("coordinates: ", sum(!is.na(meta$lat_converted)), "/", nrow(meta),
        " rows converted to decimal degrees (", n_fail,
        " non-empty values failed to parse)")

# shallow_deep from water depth: <= 200 m -> shallow, > 200 m -> deep,
# no depth_m recorded -> NA.
depth_num <- suppressWarnings(as.numeric(meta$depth_m))
derived_sd <- ifelse(is.na(depth_num), NA_character_,
                     ifelse(depth_num <= 200, "shallow", "deep"))
meta$shallow_deep <- derived_sd
message("shallow_deep: ", sum(meta$shallow_deep == "shallow", na.rm = TRUE),
        " shallow, ", sum(meta$shallow_deep == "deep", na.rm = TRUE),
        " deep, ", sum(is.na(meta$shallow_deep)), " NA (no depth_m)")

# env_2 from env_biome: ENVO codes are resolved to their labels via the
# BioPortal ENVO export (reference/ENVO_IDs.csv: ID + Preferred Label);
# plain-text values pass through unchanged; empty values and codes absent
# from the mapping become NA.
envo_map <- read.csv(envo_map_csv, colClasses = "character")
envo_labels <- setNames(envo_map$Preferred.Label, envo_map$ID)
convert_biome <- function(v) {
    v <- trimws(v)
    if (is.na(v) || v == "") return(NA_character_)
    if (grepl("^ENVO[:_][0-9]+$", v)) {
        id <- sub(":", "_", v, fixed = TRUE)
        if (id %in% names(envo_labels)) return(envo_labels[[id]])
        return(NA_character_)
    }
    v
}
meta$env_2 <- vapply(meta$env_biome, convert_biome, character(1),
                     USE.NAMES = FALSE)
n_code <- sum(grepl("^ENVO[:_]", meta$env_biome))
message("env_2: ", sum(!is.na(meta$env_2)), "/", nrow(meta), " assigned (",
        n_code, " env_biome values were ENVO codes; unmapped codes and empty",
        " env_biome give NA)")

# Report any remaining gaps, per project, after derivation (env_2,
# shallow_deep, lat/lon) has already been applied above. Public SRA metadata
# is routinely incomplete and varies project to project, so rather than
# hand-patching specific values for specific projects (which wouldn't
# generalize to a new dataset), this surfaces exactly what's missing - as
# messages in the run log, and as a standalone report file for anyone who
# wants to check completeness without digging through metadata_cleaned.csv.
report_missing <- function(df, out_path) {
    cols <- setdiff(names(df), c("run_accession", "bioproject"))
    lines <- character(0)
    for (col in cols) {
        missing <- is.na(df[[col]]) | df[[col]] == ""
        if (!any(missing)) next
        totals <- table(df$bioproject)
        by_project <- tapply(missing, df$bioproject, sum)
        by_project <- by_project[by_project > 0]
        parts <- sprintf("%s (%d/%d missing)", names(by_project), by_project,
                         totals[names(by_project)])
        lines <- c(lines, paste0(col, " -- ", paste(parts, collapse = ", ")))
    }
    if (length(lines) == 0) lines <- "no gaps found in any column"
    for (l in lines) message("missing data: ", l)

    header <- c(
        paste("# Metadata completeness report -",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
        paste("#", nrow(df), "samples across",
              length(unique(df$bioproject)), "project(s)"),
        "#"
    )
    writeLines(c(header, lines), out_path)
    message("missing data: report written to ", out_path)
}
report_missing(meta, out_missing_report)

write.csv(meta, out_meta_clean_csv, row.names = FALSE)

META <- sample_data(data.frame(meta, row.names = meta$run_accession))
saveRDS(META, out_meta_rds)
message("sample_data: ", nrow(META), " samples x ", ncol(META), " variables")

# ---- build taxonomy from placement results ----
# Sources: the EDPL/LWR-filtered clade-tree placements (apm/plagio/scuti).
# All three trees carry full lineages, up to 8 positional taxopath fields.
# Output columns match the original pp_taxonomy.tsv:
#   name, Kingdom, Phylum, Class, Subclass, Order, Family, Genus, Species, LWR

RANKS <- c("Kingdom", "Phylum", "Class", "Subclass", "Order", "Family",
           "Genus", "Species")

parse_full <- function(fields) {
    tax <- setNames(rep("", length(RANKS)), RANKS)
    n <- min(length(fields), length(RANKS))
    tax[seq_len(n)] <- fields[seq_len(n)]
    tax
}

tax_blocks <- lapply(names(placement_files), function(tree) {
    pq <- read.delim(placement_files[[tree]], colClasses = "character")
    fields <- lapply(strsplit(pq$taxopath, ";", fixed = TRUE),
                     function(x) x[x != ""])
    block <- t(vapply(fields, parse_full, character(length(RANKS))))
    rownames(block) <- pq$name
    cbind(block, LWR = pq$LWR, placement_tree = tree)
})
tax_mat <- do.call(rbind, tax_blocks)
if (anyDuplicated(rownames(tax_mat))) {
    tax_mat <- tax_mat[!duplicated(rownames(tax_mat)), , drop = FALSE]
}

tree_counts <- table(tax_mat[, "placement_tree"])
tax_mat <- tax_mat[, c(RANKS, "LWR")]

write.table(data.frame(name = rownames(tax_mat), tax_mat, check.names = FALSE),
            out_taxonomy_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

TAX_pp <- tax_table(tax_mat)
saveRDS(TAX_pp, out_tax_rds)
message("tax_table: ", nrow(tax_mat), " ASVs (",
        paste(names(tree_counts), tree_counts, sep = "=", collapse = ", "), ")")

# ---- combine count tables ----

read_counts <- function(path) {
    # skip=1 drops the "# Constructed from biom file" comment line
    as.matrix(read.delim(path, skip = 1, row.names = 1, check.names = FALSE))
}
cil_counts <- read_counts(cil_counts_tsv)
una_counts <- read_counts(una_counts_tsv)

all_samples <- union(colnames(cil_counts), colnames(una_counts))
pad_samples <- function(m) {
    missing <- setdiff(all_samples, colnames(m))
    if (length(missing) > 0) {
        m <- cbind(m, matrix(0, nrow = nrow(m), ncol = length(missing),
                             dimnames = list(rownames(m), missing)))
    }
    m[, all_samples, drop = FALSE]
}
otu_mat <- rbind(pad_samples(cil_counts), pad_samples(una_counts))

# ASV sets should be disjoint; sum as a safety net if one ever appears in both
if (anyDuplicated(rownames(otu_mat))) {
    otu_mat <- rowsum(otu_mat, group = rownames(otu_mat))
}

write.table(data.frame("#OTU ID" = rownames(otu_mat), otu_mat,
                       check.names = FALSE),
            out_counts_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

ASV_all <- otu_table(otu_mat, taxa_are_rows = TRUE)
saveRDS(ASV_all, out_otu_rds)
message("otu_table: ", ntaxa(ASV_all), " ASVs x ", nsamples(ASV_all),
        " samples, total counts ", sum(otu_mat))

# ---- generate phyloseq object ----

# (the original script pruned an undefined `ps_pp` object; fixed to `ps` here)
ps <- phyloseq(META, ASV_all, TAX_pp)

# prune (no asvs / samples with 0 counts)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
ps <- prune_samples(sample_sums(ps) > 0, ps)

saveRDS(ps, out_ps_rds)
message("phyloseq: ", ntaxa(ps), " taxa x ", nsamples(ps),
        " samples after pruning zero-count taxa/samples")
