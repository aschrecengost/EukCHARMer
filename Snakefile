# Snakemake file - batch download and process multiple SRA studies
## Modified for batch processing multiple projects
## Original adapted from Sarah Hu https://forum.qiime2.org/t/qiime2-snakemake-workflow-tutorial-18s-16s-tag-sequencing/11334

configfile: "config.yaml"

#-----SET VARIABLES-----#

# Get list of projects to process
PROJECTS = config["projects"]

# Shared paths
RAW_DATA = config["raw_data"]
OUTPUT = config["output"]
VISUALIZATION = config["visualization"]
REFFASTA = config["taxonomy_fasta"]
CLASSIFIER = config["trained_ref_database"]

# Merged paths
MERGED = OUTPUT + "_merged/"

# Placements
PLACEMENT        = OUTPUT + "_placement/"
CIL_PLACEMENT    = OUTPUT + "_cil_placement/"
APM_PLACEMENT    = OUTPUT + "_apm_placement/"
PLAGIO_PLACEMENT = OUTPUT + "_plagio_placement/"
SCUTI_PLACEMENT  = OUTPUT + "_scuti_placement/"

# R figure generation (phyloseq objects, taxa barplots, sample map)
FIGURES = OUTPUT + "_figures/"

# ---- METADATA CLEANING CONFIG ----
# For each target column name, list every raw name it might appear as across projects.
# All aliases will be renamed to the target; duplicates are merged into one column.
METADATA_RENAME = {
    "collection_date":      ["collection date", "collection_date"],
    "geo_loc_name":         ["geographic location (country and/or sea)", "geo_loc_name"],
    "depth":                ["geographic location (depth)", "depth"],
    "elev":                 ["geographic location (elevation)", "elev"],
    "env_biome":            ["environment (biome)", "env_biome"],
    "env_feature":          ["environment (feature)", "env_feature"],
    "env_material":         ["environment (material)", "env_material"],
    "depth_m":              ["depth (m)", "depth_(m)", "depth_m"],
    "depth_layer":          ["depth layer", "Depth_Layer", "depth_layer"],
    "sample_title":         ["sample_title", "sample_title_x"],
}

# Columns to keep in the cleaned output (use the target names from METADATA_RENAME above).
METADATA_KEEP = [
    "run_accession",
    "experiment_title",
    "organism_name",
    "library_strategy",
    "library_source",
    "library_selection",
    "sample_accession",
    "biosample",
    "bioproject",
    "collection_date",
    "geo_loc_name",
    "lat_lon",
    "description",
    "depth",
    "depth_m",
    "depth_layer",
    "elev",
    "env_biome",
    "env_feature",
    "env_material",
    "samp_collect_device",
    "sample_title",
]
# ---- END METADATA CLEANING CONFIG ----

# Helper function to get project-specific or global settings
def get_setting(wildcards, key, default=None):
    project = getattr(wildcards, "project", None)
    project_settings = config.get("project_settings", {})
    if project in project_settings and key in project_settings[project]:
        return project_settings[project][key]
    return config.get(key, default)


rule all:
    input:
        # Download and import - for all projects
        rawreads = expand(RAW_DATA + "{project}/reads", project=PROJECTS),
        q2_import = expand(OUTPUT + "{project}/{project}-PE-demux.qza", project=PROJECTS),

        # Processing - for all projects
        q2_primerRM = expand(OUTPUT + "{project}/{project}-PE-demux-noprimer.qza", project=PROJECTS),
        q2_merged = expand(OUTPUT + "{project}/{project}-PE-demux-noprimer-merged.qza", project=PROJECTS),
        q2_qualfilter = expand(OUTPUT + "{project}/{project}-PE-demux-noprimer-merged-filtered.qza", project=PROJECTS),

        # Core outputs - for all projects
        table = expand(OUTPUT + "{project}/{project}-table.qza", project=PROJECTS),
        rep = expand(OUTPUT + "{project}/{project}-rep-seqs.qza", project=PROJECTS),
        stats = expand(OUTPUT + "{project}/{project}-deblur-stats.qza", project=PROJECTS),
        sklearn = expand(OUTPUT + "{project}/{project}-tax_sklearn.qza", project=PROJECTS),

        # Visualizations - for all projects
        raw = expand(VISUALIZATION + "{project}/{project}-PE-demux.qzv", project=PROJECTS),
        primer = expand(VISUALIZATION + "{project}/{project}-PE-demux-noprimer.qzv", project=PROJECTS),
        merged = expand(VISUALIZATION + "{project}/{project}-PE-demux-noprimer-merged.qzv", project=PROJECTS),
        filtered = expand(VISUALIZATION + "{project}/{project}-PE-demux-noprimer-merged-filtered.qzv", project=PROJECTS),
        deblur_stats = expand(VISUALIZATION + "{project}/{project}-deblur-stats.qzv", project=PROJECTS),

        # Merged outputs
        merged_table = MERGED + "merged-table.qza",
        merged_seqs = MERGED + "merged-seqs.qza",
        merged_taxa = MERGED + "merged-taxa.qza",

        # Ciliophora merged artifacts and exports
        merged_cil_table = MERGED + "merged-ciliophora-table.qza",
        merged_cil_seqs = MERGED + "merged-ciliophora-seqs.qza",
        table_cil_biom = MERGED + "export/table/merged-ciliophora-table.biom",
        table_cil_tsv = MERGED + "export/table/merged-ciliophora-table.tsv",
        rep_seqs_cil_fasta = MERGED + "export/merged-ciliophora-seqs.fasta",

        # Unassigned merged artifacts and exports
        merged_unassigned_table = MERGED + "merged-unassigned-table.qza",
        merged_unassigned_seqs = MERGED + "merged-unassigned-seqs.qza",
        table_unassigned_biom = MERGED + "export/table/merged-unassigned-table.biom",
        table_unassigned_tsv = MERGED + "export/table/merged-unassigned-table.tsv",
        rep_seqs_unassigned_fasta = MERGED + "export/merged-unassigned-seqs.fasta",

        # Taxonomy export
        merged_table_tax = MERGED + "export/merged-taxonomy.tsv",

        # Metadata
        metadata = expand("documents/raw/{project}.csv", project=PROJECTS),
        cleaned = expand("documents/cleaned/{project}.csv", project=PROJECTS),
        merged_metadata = "documents/merged/merged_metadata.csv",

        # R figure generation
        figure_taxonomy = FIGURES + "pp_taxonomy.tsv",
        figure_counts = FIGURES + "combined_count_table.tsv",
        figure_meta_clean = FIGURES + "metadata_cleaned.csv",
        figure_missing_report = FIGURES + "metadata_gaps.txt",
        figure_metadata = FIGURES + "metadata_sample_data.rds",
        figure_tax_rds = FIGURES + "tax_table.rds",
        figure_otu_rds = FIGURES + "otu_table.rds",
        figure_ps = FIGURES + "ps.rds",
        figure_barplot = FIGURES + "taxa_barplot.pdf",
        figure_map = FIGURES + "sample_map.pdf",

        # Phylogenetic placement outputs
        # changed these to _unassigned
        unassigned_epa_jplace = PLACEMENT + "epa_result.jplace",
        unassigned_heat_tree_svg = PLACEMENT + "tree.svg",
        unassigned_per_query = PLACEMENT + "per_query.tsv",
        unassigned_lwr_histogram = PLACEMENT + "lwr-histogram.csv",
        unassigned_edpl_histogram = PLACEMENT + "edpl_histogram.csv",

        # EUK placement filtering
        per_query_ciliates = PLACEMENT + "per_query_ciliates.tsv",
        filtered_euk_placements = PLACEMENT + "filtered_ciliate_LWR_EDPL.tsv",
        unassigned_ciliates_fasta = PLACEMENT + "unassigned_ciliates.fasta",

        # Taxa placement outputs
        taxa_jplace = CIL_PLACEMENT + "epa_result.jplace",
        taxa_heat_tree_svg = CIL_PLACEMENT + "tree.svg",
        taxa_per_query = CIL_PLACEMENT + "per_query.tsv",
        taxa_lwr_histogram = CIL_PLACEMENT + "lwr-histogram.csv",
        taxa_edpl_histogram = CIL_PLACEMENT + "edpl_histogram.csv",

        # Clade filtering
        per_query_Clade1 = CIL_PLACEMENT + "per_query_Clade1.tsv",
        per_query_Clade2 = CIL_PLACEMENT + "per_query_Clade2.tsv",
        per_query_Clade3 = CIL_PLACEMENT + "per_query_Clade3.tsv",
        filtered_Clade1 = CIL_PLACEMENT + "filtered_APM_LWR_EDPL.tsv",
        filtered_Clade2 = CIL_PLACEMENT + "filtered_Plagio_LWR_EDPL.tsv",
        filtered_Clade3 = CIL_PLACEMENT + "filtered_Scuti_LWR_EDPL.tsv",

        # Clade1 placement outputs
        Clade1_jplace = APM_PLACEMENT + "epa_result.jplace",
        Clade1_heat_tree = APM_PLACEMENT + "tree.svg",
        Clade1_per_query = APM_PLACEMENT + "per_query.tsv",
        Clade1_lwr = APM_PLACEMENT + "lwr-histogram.csv",
        Clade1_edpl = APM_PLACEMENT + "edpl_histogram.csv",
        Clade1_filtered = APM_PLACEMENT + "filtered_APM_LWR_EDPL.tsv",
        Clade1_placed_seqs = APM_PLACEMENT + "apm_placed_seqs.fasta",

        # Clade2 placement outputs
        Clade2_jplace = PLAGIO_PLACEMENT + "epa_result.jplace",
        Clade2_heat_tree = PLAGIO_PLACEMENT + "tree.svg",
        Clade2_per_query = PLAGIO_PLACEMENT + "per_query.tsv",
        Clade2_lwr = PLAGIO_PLACEMENT + "lwr-histogram.csv",
        Clade2_edpl = PLAGIO_PLACEMENT + "edpl_histogram.csv",
        Clade2_filtered = PLAGIO_PLACEMENT + "filtered_Plagio_LWR_EDPL.tsv",
        Clade2_placed_seqs = PLAGIO_PLACEMENT + "plagio_placed_seqs.fasta",

        # Clade3 placement outputs
        Clade3_jplace = SCUTI_PLACEMENT + "epa_result.jplace",
        Clade3_heat_tree = SCUTI_PLACEMENT + "tree.svg",
        Clade3_per_query = SCUTI_PLACEMENT + "per_query.tsv",
        Clade3_lwr = SCUTI_PLACEMENT + "lwr-histogram.csv",
        Clade3_edpl = SCUTI_PLACEMENT + "edpl_histogram.csv",
        Clade3_filtered = SCUTI_PLACEMENT + "filtered_Scuti_LWR_EDPL.tsv",
        Clade3_placed_seqs = SCUTI_PLACEMENT + "scuti_placed_seqs.fasta",


rule download_runinfo:
    output:
        runinfo = RAW_DATA + "{project}/runinfo.csv"
    conda:
        "envs/parallelfastqdump.yaml"
    params:
        sra_id = lambda wildcards: get_setting(wildcards, "SRAid", wildcards.project)
    shell:
        """
        mkdir -p $(dirname {output.runinfo})
        esearch -db sra -query {params.sra_id} | efetch -format runinfo > {output.runinfo}
        """


rule get_runs:
    input:
        runinfo = RAW_DATA + "{project}/runinfo.csv"
    output:
        SRRnumbers = RAW_DATA + "{project}/SRR.numbers"
    shell:
        """
        cat {input.runinfo} | cut -d ',' -f 1 > {output.SRRnumbers}
        sed -i '1d' {output.SRRnumbers}
        """


rule fasterq_dump:
    input:
        SRRnumbers = RAW_DATA + "{project}/SRR.numbers"
    output:
        rawreads = directory(RAW_DATA + "{project}/reads")
    conda:
        "envs/parallelfastqdump.yaml"
    threads: 10
    shell:
        """
        mkdir -p {output.rawreads}
        cat {input.SRRnumbers} | parallel -j 1 \
        parallel-fastq-dump --sra-id {{}} \
        --threads {threads} \
        --outdir {output.rawreads} \
        --split-files --gzip
        """


rule rename_files_import:
    input:
        rawreads = RAW_DATA + "{project}/reads"
    output:
        q2_import = OUTPUT + "{project}/{project}-PE-demux.qza"
    log:
        OUTPUT + "{project}/logs/{project}_q2.log"
    conda:
        "envs/qiime2-amplicon-2026.1.yaml"
    shell:
        """
        if ls {input.rawreads}/*_1.fastq.gz 1>/dev/null 2>&1; then
            for f in {input.rawreads}/*_1.fastq.gz; do
                base="${{f%_1.fastq.gz}}"
                mv "$f"                    "${{base}}_S1_L001_R1_001.fastq.gz"
                mv "${{base}}_2.fastq.gz"  "${{base}}_S1_L001_R2_001.fastq.gz"
            done
        fi
        mkdir -p $(dirname {output.q2_import})
        qiime tools import \
            --type 'SampleData[PairedEndSequencesWithQuality]' \
            --input-path {input.rawreads} \
            --output-path {output.q2_import} \
            --input-format CasavaOneEightSingleLanePerSampleDirFmt
        """


rule rm_primers:
    input:
        q2_import = OUTPUT + "{project}/{project}-PE-demux.qza"
    output:
        q2_primerRM = OUTPUT + "{project}/{project}-PE-demux-noprimer.qza"
    log:
        OUTPUT + "{project}/logs/{project}_primer_q2.log"
    conda:
        "envs/qiime2-amplicon-2026.1.yaml"
    params:
        primerF = lambda wildcards: get_setting(wildcards, "primerF"),
        primerR = lambda wildcards: get_setting(wildcards, "primerR"),
        min_length = lambda wildcards: get_setting(wildcards, "p_min_length", 10),
        primer_err = lambda wildcards: get_setting(wildcards, "primer_err", 0.1)
    shell:
        """
        mkdir -p $(dirname {log})
        qiime cutadapt trim-paired \
           --i-demultiplexed-sequences {input.q2_import} \
           --p-front-f {params.primerF} \
           --p-front-r {params.primerR} \
           --p-minimum-length {params.min_length} \
           --p-error-rate {params.primer_err} \
           --p-discard-untrimmed \
           --p-match-adapter-wildcards \
           --p-match-read-wildcards \
           --verbose \
           --o-trimmed-sequences {output.q2_primerRM}
        """


rule merge:
    input:
        q2_primerRM = OUTPUT + "{project}/{project}-PE-demux-noprimer.qza"
    output:
        q2_merged = OUTPUT + "{project}/{project}-PE-demux-noprimer-merged.qza",
        q2_unmerged = OUTPUT + "{project}/{project}-PE-demux-noprimer-unmerged.qza"
    log:
        OUTPUT + "{project}/logs/{project}_merge_q2.log"
    conda:
        "envs/qiime2-amplicon-2026.1.yaml"
    params:
        max_diffs = lambda wildcards: get_setting(wildcards, "max_diffs_merge", 40),
        minovlen = lambda wildcards: get_setting(wildcards, "minovlen_merge", 40)
    threads: 4
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output.q2_unmerged})
        qiime vsearch merge-pairs \
            --i-demultiplexed-seqs {input.q2_primerRM} \
            --o-merged-sequences {output.q2_merged} \
            --o-unmerged-sequences {output.q2_unmerged} \
            --p-maxdiffs {params.max_diffs} \
            --p-allowmergestagger \
            --p-minovlen {params.minovlen} \
            --p-threads {threads} \
            --verbose
        """


rule quality_filter:
    input:
        q2_merged = OUTPUT + "{project}/{project}-PE-demux-noprimer-merged.qza"
    output:
        q2_qualfilter = OUTPUT + "{project}/{project}-PE-demux-noprimer-merged-filtered.qza",
        q2_qualfilter_stats = OUTPUT + "{project}/{project}-PE-demux-filter-stats.qza",
        q2_qualfilter_stats_viz = VISUALIZATION + "{project}/{project}-PE-demux-filter-stats.qzv"
    log:
        OUTPUT + "{project}/logs/{project}_filter_q2.log"
    conda:
        "envs/qiime2-amplicon-2026.1.yaml"
    params:
        min_quality = lambda wildcards: get_setting(wildcards, "filter_minquality", 20)
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output.q2_qualfilter_stats_viz})
        qiime quality-filter q-score \
            --i-demux {input.q2_merged} \
            --p-min-quality {params.min_quality} \
            --o-filtered-sequences {output.q2_qualfilter} \
            --o-filter-stats {output.q2_qualfilter_stats}
        qiime metadata tabulate \
            --m-input-file {output.q2_qualfilter_stats} \
            --o-visualization {output.q2_qualfilter_stats_viz}
        """


# IMPORT REFERENCE DB MANUALLY before running pipeline:
# unset PYTHONPATH
# conda activate qiime2-amplicon-2026.1
# qiime tools import --type 'FeatureData[Sequence]' --input-path reference/reference.fasta --output-path reference/reference_database.qza

rule deblur:
    input:
        q2_qualfilter = OUTPUT + "{project}/{project}-PE-demux-noprimer-merged-filtered.qza",
        ref_database_fasta = REFFASTA
    output:
        table = OUTPUT + "{project}/{project}-table.qza",
        rep = OUTPUT + "{project}/{project}-rep-seqs.qza",
        stats = OUTPUT + "{project}/{project}-deblur-stats.qza"
    log:
        OUTPUT + "{project}/logs/{project}-deblur.log"
    conda:
        "envs/qiime2-amplicon-2026.1.yaml"
    params:
        trim_length = lambda wildcards: get_setting(wildcards, "p_trim_length", 325)
    threads: 16
    shell:
        """
        mkdir -p $(dirname {log})
        qiime deblur denoise-other \
            --i-demultiplexed-seqs {input.q2_qualfilter} \
            --i-reference-seqs {input.ref_database_fasta} \
            --p-trim-length {params.trim_length} \
            --p-no-hashed-feature-ids \
            --o-representative-sequences {output.rep} \
            --o-table {output.table} \
            --p-sample-stats \
            --p-jobs-to-start {threads} \
            --o-stats {output.stats} \
            --p-min-reads 1
        """


rule get_stats:
    input:
        q2_import = OUTPUT + "{project}/{project}-PE-demux.qza",
        q2_primerRM = OUTPUT + "{project}/{project}-PE-demux-noprimer.qza",
        q2_merged = OUTPUT + "{project}/{project}-PE-demux-noprimer-merged.qza",
        q2_qualfilter = OUTPUT + "{project}/{project}-PE-demux-noprimer-merged-filtered.qza",
        stats = OUTPUT + "{project}/{project}-deblur-stats.qza"
    output:
        raw = VISUALIZATION + "{project}/{project}-PE-demux.qzv",
        primer = VISUALIZATION + "{project}/{project}-PE-demux-noprimer.qzv",
        merged = VISUALIZATION + "{project}/{project}-PE-demux-noprimer-merged.qzv",
        filtered = VISUALIZATION + "{project}/{project}-PE-demux-noprimer-merged-filtered.qzv",
        deblur_stats = VISUALIZATION + "{project}/{project}-deblur-stats.qzv"
    log:
        OUTPUT + "{project}/logs/{project}_getviz_q2.log"
    conda:
        "envs/qiime2-amplicon-2026.1.yaml"
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output.raw})
        qiime demux summarize --i-data {input.q2_import} --o-visualization {output.raw}
        qiime demux summarize --i-data {input.q2_primerRM} --o-visualization {output.primer}
        qiime demux summarize --i-data {input.q2_merged} --o-visualization {output.merged}
        qiime demux summarize --i-data {input.q2_qualfilter} --o-visualization {output.filtered}
        qiime deblur visualize-stats --i-deblur-stats {input.stats} --o-visualization {output.deblur_stats}
        """


rule assign_tax:
    input:
        rep = OUTPUT + "{project}/{project}-rep-seqs.qza",
        classifier = CLASSIFIER
    output:
        sklearn = OUTPUT + "{project}/{project}-tax_sklearn.qza"
    log:
        OUTPUT + "{project}/logs/{project}-sklearn_q2.log"
    conda:
        "envs/qiime2-amplicon-2026.1.yaml"
    shell:
        """
        mkdir -p $(dirname {log})
        qiime feature-classifier classify-sklearn \
            --i-classifier {input.classifier} \
            --i-reads {input.rep} \
            --o-classification {output.sklearn}
        """


rule filter_table_taxa:
    input:
        table = OUTPUT + "{project}/{project}-table.qza",
        sklearn = OUTPUT + "{project}/{project}-tax_sklearn.qza"
    output:
        table_taxa = OUTPUT + "{project}/{project}-table-taxa.qza"
    log:
        OUTPUT + "qiime2/logs/{project}/{project}-filter-table-taxa_q2.log"
    params:
        include = config["p_include_taxa"]
    conda:
        "envs/qiime2-amplicon-2026.1.yaml"
    shell:
        """
        mkdir -p $(dirname {log})
        export R_HOME=${{R_HOME:-$CONDA_PREFIX/lib/R}}
        if qiime taxa filter-table \
            --i-table {input.table} \
            --i-taxonomy {input.sklearn} \
            --p-include {params.include} \
            --o-filtered-table {output.table_taxa} 2> {log}; then
            :
        elif grep -qi "all features were filtered" {log}; then
            echo "WARNING: no {params.include} features found for {wildcards.project}; writing empty placeholder table" | tee -a {log}
            python scripts/make_empty_qza.py "FeatureTable[Frequency]" {output.table_taxa}
        else
            cat {log}
            exit 1
        fi
        """


rule filter_table_unassigned:
    input:
        table = OUTPUT + "{project}/{project}-table.qza",
        sklearn = OUTPUT + "{project}/{project}-tax_sklearn.qza"
    output:
        table_unassigned = OUTPUT + "{project}/{project}-table-unassigned.qza"
    log:
        OUTPUT + "qiime2/logs/{project}/{project}-filter-table-unassigned_q2.log"
    params:
        include = config["p_include_unassigned"]
    conda:
        "envs/qiime2-amplicon-2026.1.yaml"
    shell:
        """
        mkdir -p $(dirname {log})
        export R_HOME=${{R_HOME:-$CONDA_PREFIX/lib/R}}
        if qiime taxa filter-table \
            --i-table {input.table} \
            --i-taxonomy {input.sklearn} \
            --p-mode exact \
            --p-include "{params.include}" \
            --o-filtered-table {output.table_unassigned} 2> {log}; then
            :
        elif grep -qi "all features were filtered" {log}; then
            echo "WARNING: no Unassigned/Eukaryota features found for {wildcards.project}; writing empty placeholder table" | tee -a {log}
            python scripts/make_empty_qza.py "FeatureTable[Frequency]" {output.table_unassigned}
        else
            cat {log}
            exit 1
        fi
        """


rule filter_seqs_taxa:
    input:
        rep = OUTPUT + "{project}/{project}-rep-seqs.qza",
        sklearn = OUTPUT + "{project}/{project}-tax_sklearn.qza"
    output:
        seqs_taxa = OUTPUT + "{project}/{project}-rep-seqs-taxa.qza"
    log:
        OUTPUT + "qiime2/logs/{project}/{project}-filter-seqs-taxa_q2.log"
    params:
        include = config["p_include_taxa"]
    conda:
        "envs/qiime2-amplicon-2026.1.yaml"
    shell:
        """
        mkdir -p $(dirname {log})
        export R_HOME=${{R_HOME:-$CONDA_PREFIX/lib/R}}
        if qiime taxa filter-seqs \
            --i-sequences {input.rep} \
            --i-taxonomy {input.sklearn} \
            --p-include {params.include} \
            --o-filtered-sequences {output.seqs_taxa} 2> {log}; then
            :
        elif grep -qi "all features were filtered" {log}; then
            echo "WARNING: no {params.include} features found for {wildcards.project}; writing empty placeholder sequences" | tee -a {log}
            python scripts/make_empty_qza.py "FeatureData[Sequence]" {output.seqs_taxa}
        else
            cat {log}
            exit 1
        fi
        """


rule filter_seqs_unassigned:
    input:
        rep = OUTPUT + "{project}/{project}-rep-seqs.qza",
        sklearn = OUTPUT + "{project}/{project}-tax_sklearn.qza"
    output:
        seqs_unassigned = OUTPUT + "{project}/{project}-rep-seqs-unassigned.qza"
    log:
        OUTPUT + "qiime2/logs/{project}/{project}-filter-seqs-unassigned_q2.log"
    params:
        include = config["p_include_unassigned"]
    conda:
        "envs/qiime2-amplicon-2026.1.yaml"
    shell:
        """
        mkdir -p $(dirname {log})
        export R_HOME=${{R_HOME:-$CONDA_PREFIX/lib/R}}
        if qiime taxa filter-seqs \
            --i-sequences {input.rep} \
            --i-taxonomy {input.sklearn} \
            --p-mode exact \
            --p-include "{params.include}" \
            --o-filtered-sequences {output.seqs_unassigned} 2> {log}; then
            :
        elif grep -qi "all features were filtered" {log}; then
            echo "WARNING: no Unassigned/Eukaryota features found for {wildcards.project}; writing empty placeholder sequences" | tee -a {log}
            python scripts/make_empty_qza.py "FeatureData[Sequence]" {output.seqs_unassigned}
        else
            cat {log}
            exit 1
        fi
        """


rule merge_filtered:
    input:
        tables_cil = expand(OUTPUT + "{project}/{project}-table-taxa.qza", project=PROJECTS),
        seqs_cil = expand(OUTPUT + "{project}/{project}-rep-seqs-taxa.qza", project=PROJECTS),
        tables_unassigned = expand(OUTPUT + "{project}/{project}-table-unassigned.qza", project=PROJECTS),
        seqs_unassigned = expand(OUTPUT + "{project}/{project}-rep-seqs-unassigned.qza", project=PROJECTS)
    output:
        merged_cil_table = MERGED + "merged-ciliophora-table.qza",
        merged_cil_seqs = MERGED + "merged-ciliophora-seqs.qza",
        merged_unassigned_table = MERGED + "merged-unassigned-table.qza",
        merged_unassigned_seqs = MERGED + "merged-unassigned-seqs.qza"
    log:
        MERGED + "logs/merge_filtered_q2.log"
    conda:
        "envs/qiime2-amplicon-2026.1.yaml"
    shell:
        """
        mkdir -p $(dirname {output.merged_cil_table})
        mkdir -p $(dirname {log})
        qiime feature-table merge \
            --i-tables {input.tables_cil} \
            --o-merged-table {output.merged_cil_table}
        qiime feature-table merge-seqs \
            --i-data {input.seqs_cil} \
            --o-merged-data {output.merged_cil_seqs}
        qiime feature-table merge \
            --i-tables {input.tables_unassigned} \
            --o-merged-table {output.merged_unassigned_table}
        qiime feature-table merge-seqs \
            --i-data {input.seqs_unassigned} \
            --o-merged-data {output.merged_unassigned_seqs}
        """


rule merge_projects:
    input:
        tables = expand(OUTPUT + "{project}/{project}-table.qza", project=PROJECTS),
        seqs = expand(OUTPUT + "{project}/{project}-rep-seqs.qza", project=PROJECTS),
        taxa = expand(OUTPUT + "{project}/{project}-tax_sklearn.qza", project=PROJECTS)
    output:
        merged_table = MERGED + "merged-table.qza",
        merged_seqs = MERGED + "merged-seqs.qza",
        merged_taxa = MERGED + "merged-taxa.qza"
    log:
        MERGED + "logs/merge_projects_q2.log"
    conda:
        "envs/qiime2-amplicon-2026.1.yaml"
    shell:
        """
        mkdir -p $(dirname {output.merged_table})
        mkdir -p $(dirname {log})
        qiime feature-table merge \
            $(echo "{input.tables}" | tr ' ' '\\n' | sed 's/^/--i-tables /') \
            --o-merged-table {output.merged_table}
        qiime feature-table merge-seqs \
            $(echo "{input.seqs}" | tr ' ' '\\n' | sed 's/^/--i-data /') \
            --o-merged-data {output.merged_seqs}
        qiime feature-table merge-taxa \
            $(echo "{input.taxa}" | tr ' ' '\\n' | sed 's/^/--i-data /') \
            --o-merged-data {output.merged_taxa}
        """


rule export_ciliophora:
    input:
        table_cil_merged = MERGED + "merged-ciliophora-table.qza",
        rep_cil_merged = MERGED + "merged-ciliophora-seqs.qza"
    output:
        table_cil_biom = MERGED + "export/table/merged-ciliophora-table.biom",
        table_cil_tsv = MERGED + "export/table/merged-ciliophora-table.tsv",
        rep_seqs_cil_fasta = MERGED + "export/merged-ciliophora-seqs.fasta"
    log:
        MERGED + "logs/export_ciliophora_q2.log"
    conda:
        "envs/qiime2-amplicon-2026.1.yaml"
    params:
        table_outdir = MERGED + "export/table/",
        seqs_outdir = MERGED + "export/",
        # private temp dirs so a concurrently-running export_unassigned job can't
        # clobber our intermediate feature-table.biom / dna-sequences.fasta
        table_tmpdir = MERGED + "export/table/.tmp_ciliophora/",
        seqs_tmpdir = MERGED + "export/.tmp_ciliophora_seqs/"
    shell:
        """
        mkdir -p {params.table_outdir}
        mkdir -p {params.seqs_outdir}
        mkdir -p {params.table_tmpdir}
        mkdir -p {params.seqs_tmpdir}
        mkdir -p $(dirname {log})
        qiime tools export --input-path {input.table_cil_merged} --output-path {params.table_tmpdir}
        mv {params.table_tmpdir}feature-table.biom {output.table_cil_biom}
        rmdir {params.table_tmpdir}
        biom convert -i {output.table_cil_biom} -o {output.table_cil_tsv} --to-tsv
        qiime tools export --input-path {input.rep_cil_merged} --output-path {params.seqs_tmpdir}
        mv {params.seqs_tmpdir}dna-sequences.fasta {output.rep_seqs_cil_fasta}
        rmdir {params.seqs_tmpdir}
        """


rule export_unassigned:
    input:
        table_unassigned_merged = MERGED + "merged-unassigned-table.qza",
        rep_unassigned_merged = MERGED + "merged-unassigned-seqs.qza"
    output:
        table_unassigned_biom = MERGED + "export/table/merged-unassigned-table.biom",
        table_unassigned_tsv = MERGED + "export/table/merged-unassigned-table.tsv",
        rep_seqs_unassigned_fasta = MERGED + "export/merged-unassigned-seqs.fasta"
    log:
        MERGED + "logs/export_unassigned_q2.log"
    conda:
        "envs/qiime2-amplicon-2026.1.yaml"
    params:
        table_outdir = MERGED + "export/table/",
        seqs_outdir = MERGED + "export/",
        # private temp dirs so a concurrently-running export_ciliophora job can't
        # clobber our intermediate feature-table.biom / dna-sequences.fasta
        table_tmpdir = MERGED + "export/table/.tmp_unassigned/",
        seqs_tmpdir = MERGED + "export/.tmp_unassigned_seqs/"
    shell:
        """
        mkdir -p {params.table_outdir}
        mkdir -p {params.seqs_outdir}
        mkdir -p {params.table_tmpdir}
        mkdir -p {params.seqs_tmpdir}
        mkdir -p $(dirname {log})
        qiime tools export --input-path {input.table_unassigned_merged} --output-path {params.table_tmpdir}
        mv {params.table_tmpdir}feature-table.biom {output.table_unassigned_biom}
        rmdir {params.table_tmpdir}
        biom convert -i {output.table_unassigned_biom} -o {output.table_unassigned_tsv} --to-tsv
        qiime tools export --input-path {input.rep_unassigned_merged} --output-path {params.seqs_tmpdir}
        mv {params.seqs_tmpdir}dna-sequences.fasta {output.rep_seqs_unassigned_fasta}
        rmdir {params.seqs_tmpdir}
        """


rule export_taxonomy:
    input:
        merged_sklearn = MERGED + "merged-taxa.qza"
    output:
        merged_table_tax = MERGED + "export/merged-taxonomy.tsv"
    log:
        MERGED + "logs/export_taxonomy_q2.log"
    params:
        outdir = MERGED + "export/"
    conda:
        "envs/qiime2-amplicon-2026.1.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})
        qiime tools export --input-path {input.merged_sklearn} --output-path {params.outdir}
        mv {params.outdir}taxonomy.tsv {output.merged_table_tax}
        """

rule get_metadata:
    output:
        metadata = "documents/raw/{project}.csv"
    conda:
        "envs/pysradb.yaml"
    params:
        sra_id = lambda wildcards: get_setting(wildcards, "SRAid", wildcards.project)
    shell:
        """
        mkdir -p documents/raw
        for attempt in 1 2 3 4 5; do
            if pysradb metadata {params.sra_id} --detailed --saveto {output.metadata}; then
                exit 0
            fi
            echo "Attempt $attempt failed (likely rate limit), retrying in $((attempt * 10))s..."
            sleep $((attempt * 10))
        done
        echo "All attempts failed"
        exit 1
        """


rule clean_metadata:
    input:
        metadata = "documents/raw/{project}.csv"
    output:
        cleaned = "documents/cleaned/{project}.csv"
    params:
        rename = METADATA_RENAME,
        keep   = METADATA_KEEP
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/clean_metadata.py"


rule merge_metadata:
    input:
        cleaned = expand("documents/cleaned/{project}.csv", project=PROJECTS)
    output:
        merged = "documents/merged/merged_metadata.csv"
    run:
        import pandas as pd
        import os
        os.makedirs("documents/merged", exist_ok=True)
        df = pd.concat([pd.read_csv(f) for f in input.cleaned], ignore_index=True)
        df.to_csv(output.merged, index=False)
        print(f"Merged {len(input.cleaned)} projects, {len(df)} total rows")

rule papara_unassigned_seqs:
    input:
        rep_seqs_unassigned_fasta = MERGED + "export/merged-unassigned-seqs.fasta",
        euk_tree = config["EUK_TREE"],
        msa_phylip = config["MSA_PHYLIP_EUK"],
        msa_fasta = config["MSA_FASTA_EUK"],
        clades = config["CLADES_EUKS"]
    output:
        papara_alignment = PLACEMENT + "papara_alignment.default"
    params:
        placement_dir = PLACEMENT
    threads: 17
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}

        cp {input.euk_tree} {params.placement_dir}
        cp {input.msa_phylip} {params.placement_dir}
        cp {input.msa_fasta} {params.placement_dir}
        cp {input.clades} {params.placement_dir}
        cp {input.rep_seqs_unassigned_fasta} {params.placement_dir}

        cd {params.placement_dir}

        module load papara_nt/2.5
        papara \
            -t $(basename {input.euk_tree}) \
            -s $(basename {input.msa_phylip}) \
            -q $(basename {input.rep_seqs_unassigned_fasta}) \
            -j {threads} \
            -r

        """


rule epa_split:
    input:
        papara_alignment = PLACEMENT + "papara_alignment.default",
        msa_fasta = config["MSA_FASTA_EUK"]
    output:
        query_fasta = PLACEMENT + "query.fasta",
        ref_fasta = PLACEMENT + "reference.fasta"
    params:
        placement_dir = PLACEMENT
    threads: 17
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}
        cp {input.msa_fasta} {params.placement_dir}

        cd {params.placement_dir}

        epa-ng \
            --split $(basename {input.msa_fasta}) \
            $(basename {input.papara_alignment})

        """


rule raxml_evaluate:
    input:
        ref_fasta = PLACEMENT + "reference.fasta",
        euk_tree = config["EUK_TREE"]
    output:
        best_model = PLACEMENT + "reference.fasta.raxml.bestModel",
        best_tree = PLACEMENT + "reference.fasta.raxml.bestTree",
        log = PLACEMENT + "reference.fasta.raxml.log",
        rba = PLACEMENT + "reference.fasta.raxml.rba",
        start_tree = PLACEMENT + "reference.fasta.raxml.startTree"
    params:
        placement_dir = PLACEMENT
    threads: 17
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}
        cp {input.euk_tree} {params.placement_dir}

        cd {params.placement_dir}

        raxml-ng \
            --evaluate \
            --msa $(basename {input.ref_fasta}) \
            --tree $(basename {input.euk_tree}) \
            --model GTR+G \
            --threads {threads} --force perf_threads

        """


rule epa_placement:
    input:
        euk_tree = config["EUK_TREE"],
        ref_fasta = PLACEMENT + "reference.fasta",
        query_fasta = PLACEMENT + "query.fasta",
        best_model = PLACEMENT + "reference.fasta.raxml.bestModel"
    output:
        jplace = PLACEMENT + "epa_result.jplace",
        info_log = PLACEMENT + "epa_info.log"
    params:
        placement_dir = PLACEMENT,
        filter_acc_lwr = config["epa_filter_acc_lwr"],
        filter_max = config["epa_filter_max"]
    threads: 17
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}
        cp {input.euk_tree} {params.placement_dir}

        # Remove any leftover epa-ng outputs from previous runs to avoid --redo errors
        rm -f {params.placement_dir}epa_result.jplace \
              {params.placement_dir}epa_info.log

        cd {params.placement_dir}

        epa-ng \
            --filter-acc-lwr {params.filter_acc_lwr} \
            --filter-max {params.filter_max} \
            -t $(basename {input.euk_tree}) \
            -s $(basename {input.ref_fasta}) \
            -q $(basename {input.query_fasta}) \
            --model $(basename {input.best_model})

        """


rule gappa_heat_tree:
    input:
        jplace = PLACEMENT + "epa_result.jplace"
    output:
        svg = PLACEMENT + "tree.svg",
        newick = PLACEMENT + "tree.newick",
        nexus = PLACEMENT + "tree.nexus"
    params:
        placement_dir = PLACEMENT,
        mass_norm = config["gappa_mass_norm"]
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}

        cd {params.placement_dir}

        gappa examine heat-tree \
            --jplace-path $(basename {input.jplace}) \
            --mass-norm {params.mass_norm} \
            --allow-file-overwriting \
            --write-svg-tree \
            --write-newick-tree \
            --write-nexus-tree

        """



rule gappa_assign:
    input:
        jplace = PLACEMENT + "epa_result.jplace",
        clades = config["CLADES_EUKS"]
    output:
        per_query = PLACEMENT + "per_query.tsv",
        profile = PLACEMENT + "profile.tsv",
        labelled_tree = PLACEMENT + "labelled_tree.newick"
    params:
        placement_dir = PLACEMENT,
        consensus_thresh = config["gappa_consensus_thresh"]
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}
        cp {input.clades} {params.placement_dir}

        cd {params.placement_dir}

        gappa examine assign \
            --jplace-path $(basename {input.jplace}) \
            --consensus-thresh {params.consensus_thresh} \
            --taxon-file $(basename {input.clades}) \
            --best-hit \
            --per-query-results \
            --allow-file-overwriting

        """


rule gappa_lwr_edpl:
    input:
        jplace = PLACEMENT + "epa_result.jplace"
    output:
        lwr_histogram = PLACEMENT + "lwr-histogram.csv",
        lwr_list = PLACEMENT + "lwr-list.csv",
        edpl_histogram = PLACEMENT + "edpl_histogram.csv",
        edpl_list = PLACEMENT + "edpl_list.csv"
    params:
        placement_dir = PLACEMENT
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}

        cd {params.placement_dir}

        gappa examine lwr-histogram \
            --jplace-path $(basename {input.jplace}) \
            --allow-file-overwriting

        gappa examine lwr-list \
            --jplace-path $(basename {input.jplace}) \
            --allow-file-overwriting

        gappa examine edpl \
            --jplace-path $(basename {input.jplace}) \
            --allow-file-overwriting

        """


# ---- EUK PLACEMENT FILTERING ----

rule filter_ciliophora:
    input:
        tsv = PLACEMENT + "per_query.tsv"
    output:
        tsv = PLACEMENT + "per_query_ciliates.tsv"
    params:
        search_term = config["search_term_taxa"]
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/filter_by_taxopath.py"


rule filter_euk_placements:
    input:
        per_query = PLACEMENT + "per_query_ciliates.tsv",
        lwr_list  = PLACEMENT + "lwr-list.csv",
        edpl_list = PLACEMENT + "edpl_list.csv"
    output:
        tsv = PLACEMENT + "filtered_ciliate_LWR_EDPL.tsv"
    params:
        edpl_threshold = config["edpl_threshold"]
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/filter_placements.py"


rule extract_unassigned_ciliates:
    input:
        tsv   = PLACEMENT + "filtered_ciliate_LWR_EDPL.tsv",
        fasta = MERGED + "export/merged-unassigned-seqs.fasta"
    output:
        fasta = PLACEMENT + "unassigned_ciliates.fasta"
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/extract_fasta.py"


# ---- CILIATE PLACEMENT ----

rule combine_ciliate_seqs:
    input:
        unassigned_ciliates = PLACEMENT + "unassigned_ciliates.fasta",
        asvs_ciliates       = MERGED + "export/merged-ciliophora-seqs.fasta"
    output:
        combined = CIL_PLACEMENT + "ciliate_seqs.fasta"
    shell:
        """
        mkdir -p {CIL_PLACEMENT}
        cat {input.unassigned_ciliates} {input.asvs_ciliates} > {output.combined}
        """


rule papara_ciliate:
    input:
        ciliate_seqs = CIL_PLACEMENT + "ciliate_seqs.fasta",
        cil_tree     = config["CIL_TREE"],
        msa_phylip   = config["MSA_PHYLIP_CIL"],
        msa_fasta    = config["MSA_FASTA_CIL"],
        clades       = config["CLADES_CIL"]
    output:
        papara_alignment = CIL_PLACEMENT + "papara_alignment.default"
    params:
        placement_dir = CIL_PLACEMENT
    threads: 17
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}

        cp {input.cil_tree}   {params.placement_dir}
        cp {input.msa_phylip} {params.placement_dir}
        cp {input.msa_fasta}  {params.placement_dir}
        cp {input.clades}     {params.placement_dir}

        cd {params.placement_dir}

        module load papara_nt/2.5
        papara \
            -t $(basename {input.cil_tree}) \
            -s $(basename {input.msa_phylip}) \
            -q $(basename {input.ciliate_seqs}) \
            -j {threads} \
            -r

        """


rule epa_split_ciliate:
    input:
        papara_alignment = CIL_PLACEMENT + "papara_alignment.default",
        msa_fasta        = config["MSA_FASTA_CIL"]
    output:
        query_fasta = CIL_PLACEMENT + "query.fasta",
        ref_fasta   = CIL_PLACEMENT + "reference.fasta"
    params:
        placement_dir = CIL_PLACEMENT
    threads: 17
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}
        cp {input.msa_fasta}        {params.placement_dir}

        cd {params.placement_dir}

        epa-ng \
            --split $(basename {input.msa_fasta}) \
            $(basename {input.papara_alignment})

        """


rule raxml_evaluate_ciliate:
    input:
        ref_fasta = CIL_PLACEMENT + "reference.fasta",
        cil_tree  = config["CIL_TREE"]
    output:
        best_model = CIL_PLACEMENT + "reference.fasta.raxml.bestModel",
        best_tree  = CIL_PLACEMENT + "reference.fasta.raxml.bestTree",
        log        = CIL_PLACEMENT + "reference.fasta.raxml.log",
        rba        = CIL_PLACEMENT + "reference.fasta.raxml.rba",
        start_tree = CIL_PLACEMENT + "reference.fasta.raxml.startTree"
    params:
        placement_dir = CIL_PLACEMENT
    threads: 17
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}
        cp {input.cil_tree}  {params.placement_dir}

        cd {params.placement_dir}

        raxml-ng \
            --evaluate \
            --msa $(basename {input.ref_fasta}) \
            --tree $(basename {input.cil_tree}) \
            --model GTR+G \
            --threads {threads} --force perf_threads

        """


rule epa_placement_ciliate:
    input:
        cil_tree   = config["CIL_TREE"],
        ref_fasta  = CIL_PLACEMENT + "reference.fasta",
        query_fasta = CIL_PLACEMENT + "query.fasta",
        best_model = CIL_PLACEMENT + "reference.fasta.raxml.bestModel"
    output:
        jplace   = CIL_PLACEMENT + "epa_result.jplace",
        info_log = CIL_PLACEMENT + "epa_info.log"
    params:
        placement_dir   = CIL_PLACEMENT,
        filter_acc_lwr  = config["epa_filter_acc_lwr"],
        filter_max      = config["epa_filter_max"]
    threads: 17
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}
        cp {input.cil_tree}    {params.placement_dir}

        rm -f {params.placement_dir}epa_result.jplace \
              {params.placement_dir}epa_info.log

        cd {params.placement_dir}

        epa-ng \
            --filter-acc-lwr {params.filter_acc_lwr} \
            --filter-max {params.filter_max} \
            -t $(basename {input.cil_tree}) \
            -s $(basename {input.ref_fasta}) \
            -q $(basename {input.query_fasta}) \
            --model $(basename {input.best_model})

        """


rule gappa_heat_tree_ciliate:
    input:
        jplace = CIL_PLACEMENT + "epa_result.jplace"
    output:
        svg    = CIL_PLACEMENT + "tree.svg",
        newick = CIL_PLACEMENT + "tree.newick",
        nexus  = CIL_PLACEMENT + "tree.nexus"
    params:
        placement_dir = CIL_PLACEMENT,
        mass_norm   = config["gappa_mass_norm"]
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}

        cd {params.placement_dir}

        gappa examine heat-tree \
            --jplace-path $(basename {input.jplace}) \
            --mass-norm {params.mass_norm} \
            --allow-file-overwriting \
            --write-svg-tree \
            --write-newick-tree \
            --write-nexus-tree

        """


rule gappa_assign_ciliate:
    input:
        jplace = CIL_PLACEMENT + "epa_result.jplace",
        clades = config["CLADES_CIL"]
    output:
        per_query    = CIL_PLACEMENT + "per_query.tsv",
        profile      = CIL_PLACEMENT + "profile.tsv",
        labelled_tree = CIL_PLACEMENT + "labelled_tree.newick"
    params:
        placement_dir = CIL_PLACEMENT,
        consensus_thresh  = config["gappa_consensus_thresh"]
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}
        cp {input.clades} {params.placement_dir}

        cd {params.placement_dir}

        gappa examine assign \
            --jplace-path $(basename {input.jplace}) \
            --consensus-thresh {params.consensus_thresh} \
            --taxon-file $(basename {input.clades}) \
            --best-hit \
            --per-query-results \
            --allow-file-overwriting

        """


rule gappa_lwr_edpl_ciliate:
    input:
        jplace = CIL_PLACEMENT + "epa_result.jplace"
    output:
        lwr_histogram  = CIL_PLACEMENT + "lwr-histogram.csv",
        lwr_list       = CIL_PLACEMENT + "lwr-list.csv",
        edpl_histogram = CIL_PLACEMENT + "edpl_histogram.csv",
        edpl_list      = CIL_PLACEMENT + "edpl_list.csv"
    params:
        placement_dir = CIL_PLACEMENT,
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}

        cd {params.placement_dir}

        gappa examine lwr-histogram \
            --jplace-path $(basename {input.jplace}) \
            --allow-file-overwriting

        gappa examine lwr-list \
            --jplace-path $(basename {input.jplace}) \
            --allow-file-overwriting

        gappa examine edpl \
            --jplace-path $(basename {input.jplace}) \
            --allow-file-overwriting

        """


# ---- CLADE FILTERING ----

rule filter_Clade1:
    input:
        tsv = CIL_PLACEMENT + "per_query.tsv"
    output:
        tsv = CIL_PLACEMENT + "per_query_Clade1.tsv"
    params:
        search_term = config["search_term_clade1"]
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/filter_by_taxopath.py"


rule filter_Clade2:
    input:
        tsv = CIL_PLACEMENT + "per_query.tsv"
    output:
        tsv = CIL_PLACEMENT + "per_query_Clade2.tsv"
    params:
        search_term = config["search_term_clade2"]
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/filter_by_taxopath.py"


rule filter_Clade3:
    input:
        tsv = CIL_PLACEMENT + "per_query.tsv"
    output:
        tsv = CIL_PLACEMENT + "per_query_Clade3.tsv"
    params:
        search_term = config["search_term_clade3"]
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/filter_by_taxopath.py"


rule filter_APM_placements:
    input:
        per_query = CIL_PLACEMENT + "per_query_Clade1.tsv",
        lwr_list  = CIL_PLACEMENT + "lwr-list.csv",
        edpl_list = CIL_PLACEMENT + "edpl_list.csv"
    output:
        tsv = CIL_PLACEMENT + "filtered_APM_LWR_EDPL.tsv"
    params:
        edpl_threshold = config["edpl_threshold"]
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/filter_placements.py"


rule filter_Plagio_placements:
    input:
        per_query = CIL_PLACEMENT + "per_query_Clade2.tsv",
        lwr_list  = CIL_PLACEMENT + "lwr-list.csv",
        edpl_list = CIL_PLACEMENT + "edpl_list.csv"
    output:
        tsv = CIL_PLACEMENT + "filtered_Plagio_LWR_EDPL.tsv"
    params:
        edpl_threshold = config["edpl_threshold"]
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/filter_placements.py"


rule filter_Scuti_placements:
    input:
        per_query = CIL_PLACEMENT + "per_query_Clade3.tsv",
        lwr_list  = CIL_PLACEMENT + "lwr-list.csv",
        edpl_list = CIL_PLACEMENT + "edpl_list.csv"
    output:
        tsv = CIL_PLACEMENT + "filtered_Scuti_LWR_EDPL.tsv"
    params:
        edpl_threshold = config["edpl_threshold"]
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/filter_placements.py"


# ---- APM PLACEMENT ----

rule extract_apm_seqs:
    input:
        tsv   = CIL_PLACEMENT + "filtered_APM_LWR_EDPL.tsv",
        fasta = CIL_PLACEMENT + "ciliate_seqs.fasta"
    output:
        fasta = APM_PLACEMENT + "apm_seqs.fasta"
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/extract_fasta.py"


rule papara_apm:
    input:
        query_seqs = APM_PLACEMENT + "apm_seqs.fasta",
        tree       = config["APM_TREE"],
        msa_phylip = config["MSA_PHYLIP_APM"],
        msa_fasta  = config["MSA_FASTA_APM"],
        clades     = config["CLADES_APM"]
    output:
        papara_alignment = APM_PLACEMENT + "papara_alignment.default"
    params:
        placement_dir = APM_PLACEMENT
    threads: 17
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}

        cp {input.tree}       {params.placement_dir}
        cp {input.msa_phylip} {params.placement_dir}
        cp {input.msa_fasta}  {params.placement_dir}
        cp {input.clades}     {params.placement_dir}

        cd {params.placement_dir}

        module load papara_nt/2.5
        papara \
            -t $(basename {input.tree}) \
            -s $(basename {input.msa_phylip}) \
            -q $(basename {input.query_seqs}) \
            -j {threads} \
            -r

        """


rule epa_split_apm:
    input:
        papara_alignment = APM_PLACEMENT + "papara_alignment.default",
        msa_fasta        = config["MSA_FASTA_APM"]
    output:
        query_fasta = APM_PLACEMENT + "query.fasta",
        ref_fasta   = APM_PLACEMENT + "reference.fasta"
    params:
        placement_dir = APM_PLACEMENT
    threads: 17
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}
        cp {input.msa_fasta}        {params.placement_dir}

        cd {params.placement_dir}

        epa-ng \
            --split $(basename {input.msa_fasta}) \
            $(basename {input.papara_alignment})

        """


rule raxml_evaluate_apm:
    input:
        ref_fasta = APM_PLACEMENT + "reference.fasta",
        tree      = config["APM_TREE"]
    output:
        best_model = APM_PLACEMENT + "reference.fasta.raxml.bestModel",
        best_tree  = APM_PLACEMENT + "reference.fasta.raxml.bestTree",
        log        = APM_PLACEMENT + "reference.fasta.raxml.log",
        rba        = APM_PLACEMENT + "reference.fasta.raxml.rba",
        start_tree = APM_PLACEMENT + "reference.fasta.raxml.startTree"
    params:
        placement_dir = APM_PLACEMENT
    threads: 17
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}
        cp {input.tree}      {params.placement_dir}

        cd {params.placement_dir}

        raxml-ng \
            --evaluate \
            --msa $(basename {input.ref_fasta}) \
            --tree $(basename {input.tree}) \
            --model GTR+G \
            --threads {threads} --force perf_threads

        """


rule epa_placement_apm:
    input:
        tree        = config["APM_TREE"],
        ref_fasta   = APM_PLACEMENT + "reference.fasta",
        query_fasta = APM_PLACEMENT + "query.fasta",
        best_model  = APM_PLACEMENT + "reference.fasta.raxml.bestModel"
    output:
        jplace   = APM_PLACEMENT + "epa_result.jplace",
        info_log = APM_PLACEMENT + "epa_info.log"
    params:
        placement_dir = APM_PLACEMENT,
        filter_acc_lwr = config["epa_filter_acc_lwr"],
        filter_max     = config["epa_filter_max"]
    threads: 17
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}
        cp {input.tree}        {params.placement_dir}

        rm -f {params.placement_dir}epa_result.jplace {params.placement_dir}epa_info.log

        cd {params.placement_dir}

        epa-ng \
            --filter-acc-lwr {params.filter_acc_lwr} \
            --filter-max {params.filter_max} \
            -t $(basename {input.tree}) \
            -s $(basename {input.ref_fasta}) \
            -q $(basename {input.query_fasta}) \
            --model $(basename {input.best_model})

        """


rule gappa_heat_tree_apm:
    input:
        jplace = APM_PLACEMENT + "epa_result.jplace"
    output:
        svg    = APM_PLACEMENT + "tree.svg",
        newick = APM_PLACEMENT + "tree.newick",
        nexus  = APM_PLACEMENT + "tree.nexus"
    params:
        placement_dir = APM_PLACEMENT,
        mass_norm   = config["gappa_mass_norm"]
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}

        cd {params.placement_dir}

        gappa examine heat-tree \
            --jplace-path $(basename {input.jplace}) \
            --mass-norm {params.mass_norm} \
            --allow-file-overwriting \
            --write-svg-tree \
            --write-newick-tree \
            --write-nexus-tree

        """


rule gappa_assign_apm:
    input:
        jplace = APM_PLACEMENT + "epa_result.jplace",
        clades = config["CLADES_APM"]
    output:
        per_query     = APM_PLACEMENT + "per_query.tsv",
        profile       = APM_PLACEMENT + "profile.tsv",
        labelled_tree = APM_PLACEMENT + "labelled_tree.newick"
    params:
        placement_dir = APM_PLACEMENT,
        consensus_thresh = config["gappa_consensus_thresh"]
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}
        cp {input.clades} {params.placement_dir}

        cd {params.placement_dir}

        gappa examine assign \
            --jplace-path $(basename {input.jplace}) \
            --consensus-thresh {params.consensus_thresh} \
            --taxon-file $(basename {input.clades}) \
            --best-hit \
            --per-query-results \
            --allow-file-overwriting

        """


rule gappa_lwr_edpl_apm:
    input:
        jplace = APM_PLACEMENT + "epa_result.jplace"
    output:
        lwr_histogram  = APM_PLACEMENT + "lwr-histogram.csv",
        lwr_list       = APM_PLACEMENT + "lwr-list.csv",
        edpl_histogram = APM_PLACEMENT + "edpl_histogram.csv",
        edpl_list      = APM_PLACEMENT + "edpl_list.csv"
    params:
        placement_dir = APM_PLACEMENT,
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}

        cd {params.placement_dir}

        gappa examine lwr-histogram \
            --jplace-path $(basename {input.jplace}) \
            --allow-file-overwriting

        gappa examine lwr-list \
            --jplace-path $(basename {input.jplace}) \
            --allow-file-overwriting

        gappa examine edpl \
            --jplace-path $(basename {input.jplace}) \
            --allow-file-overwriting

        """


# ---- APM PLACEMENT FILTERING (replicates filter_euk_placements / extract_unassigned_ciliates
#      on the APM tree's own placement results) ----

rule filter_apm_final:
    input:
        per_query = APM_PLACEMENT + "per_query.tsv",
        lwr_list  = APM_PLACEMENT + "lwr-list.csv",
        edpl_list = APM_PLACEMENT + "edpl_list.csv"
    output:
        tsv = APM_PLACEMENT + "filtered_APM_LWR_EDPL.tsv"
    params:
        edpl_threshold = config["edpl_threshold"]
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/filter_placements.py"


rule extract_apm_final_seqs:
    input:
        tsv   = APM_PLACEMENT + "filtered_APM_LWR_EDPL.tsv",
        fasta = APM_PLACEMENT + "apm_seqs.fasta"
    output:
        fasta = APM_PLACEMENT + "apm_placed_seqs.fasta"
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/extract_fasta.py"


# ---- PLAGIOPYLEA PLACEMENT ----

rule extract_plagio_seqs:
    input:
        tsv   = CIL_PLACEMENT + "filtered_Plagio_LWR_EDPL.tsv",
        fasta = CIL_PLACEMENT + "ciliate_seqs.fasta"
    output:
        fasta = PLAGIO_PLACEMENT + "plagio_seqs.fasta"
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/extract_fasta.py"


rule papara_plagio:
    input:
        query_seqs = PLAGIO_PLACEMENT + "plagio_seqs.fasta",
        tree       = config["PLAGIO_TREE"],
        msa_phylip = config["MSA_PHYLIP_PLAGIO"],
        msa_fasta  = config["MSA_FASTA_PLAGIO"],
        clades     = config["CLADES_PLAGIO"]
    output:
        papara_alignment = PLAGIO_PLACEMENT + "papara_alignment.default"
    params:
        placement_dir = PLAGIO_PLACEMENT
    threads: 17
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}

        cp {input.tree}       {params.placement_dir}
        cp {input.msa_phylip} {params.placement_dir}
        cp {input.msa_fasta}  {params.placement_dir}
        cp {input.clades}     {params.placement_dir}

        cd {params.placement_dir}

        module load papara_nt/2.5
        papara \
            -t $(basename {input.tree}) \
            -s $(basename {input.msa_phylip}) \
            -q $(basename {input.query_seqs}) \
            -j {threads} \
            -r

        """


rule epa_split_plagio:
    input:
        papara_alignment = PLAGIO_PLACEMENT + "papara_alignment.default",
        msa_fasta        = config["MSA_FASTA_PLAGIO"]
    output:
        query_fasta = PLAGIO_PLACEMENT + "query.fasta",
        ref_fasta   = PLAGIO_PLACEMENT + "reference.fasta"
    params:
        placement_dir = PLAGIO_PLACEMENT
    threads: 17
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}
        cp {input.msa_fasta}        {params.placement_dir}

        cd {params.placement_dir}

        epa-ng \
            --split $(basename {input.msa_fasta}) \
            $(basename {input.papara_alignment})

        """


rule raxml_evaluate_plagio:
    input:
        ref_fasta = PLAGIO_PLACEMENT + "reference.fasta",
        tree      = config["PLAGIO_TREE"]
    output:
        best_model = PLAGIO_PLACEMENT + "reference.fasta.raxml.bestModel",
        best_tree  = PLAGIO_PLACEMENT + "reference.fasta.raxml.bestTree",
        log        = PLAGIO_PLACEMENT + "reference.fasta.raxml.log",
        rba        = PLAGIO_PLACEMENT + "reference.fasta.raxml.rba",
        start_tree = PLAGIO_PLACEMENT + "reference.fasta.raxml.startTree"
    params:
        placement_dir = PLAGIO_PLACEMENT
    threads: 17
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}
        cp {input.tree}      {params.placement_dir}

        cd {params.placement_dir}

        raxml-ng \
            --evaluate \
            --msa $(basename {input.ref_fasta}) \
            --tree $(basename {input.tree}) \
            --model GTR+G \
            --threads {threads} --force perf_threads

        """


rule epa_placement_plagio:
    input:
        tree        = config["PLAGIO_TREE"],
        ref_fasta   = PLAGIO_PLACEMENT + "reference.fasta",
        query_fasta = PLAGIO_PLACEMENT + "query.fasta",
        best_model  = PLAGIO_PLACEMENT + "reference.fasta.raxml.bestModel"
    output:
        jplace   = PLAGIO_PLACEMENT + "epa_result.jplace",
        info_log = PLAGIO_PLACEMENT + "epa_info.log"
    params:
        placement_dir = PLAGIO_PLACEMENT,
        filter_acc_lwr = config["epa_filter_acc_lwr"],
        filter_max     = config["epa_filter_max"]
    threads: 17
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}
        cp {input.tree}        {params.placement_dir}

        rm -f {params.placement_dir}epa_result.jplace {params.placement_dir}epa_info.log

        cd {params.placement_dir}

        epa-ng \
            --filter-acc-lwr {params.filter_acc_lwr} \
            --filter-max {params.filter_max} \
            -t $(basename {input.tree}) \
            -s $(basename {input.ref_fasta}) \
            -q $(basename {input.query_fasta}) \
            --model $(basename {input.best_model})

        """


rule gappa_heat_tree_plagio:
    input:
        jplace = PLAGIO_PLACEMENT + "epa_result.jplace"
    output:
        svg    = PLAGIO_PLACEMENT + "tree.svg",
        newick = PLAGIO_PLACEMENT + "tree.newick",
        nexus  = PLAGIO_PLACEMENT + "tree.nexus"
    params:
        placement_dir = PLAGIO_PLACEMENT,
        mass_norm   = config["gappa_mass_norm"]
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}

        cd {params.placement_dir}

        gappa examine heat-tree \
            --jplace-path $(basename {input.jplace}) \
            --mass-norm {params.mass_norm} \
            --allow-file-overwriting \
            --write-svg-tree \
            --write-newick-tree \
            --write-nexus-tree

        """


rule gappa_assign_plagio:
    input:
        jplace = PLAGIO_PLACEMENT + "epa_result.jplace",
        clades = config["CLADES_PLAGIO"]
    output:
        per_query     = PLAGIO_PLACEMENT + "per_query.tsv",
        profile       = PLAGIO_PLACEMENT + "profile.tsv",
        labelled_tree = PLAGIO_PLACEMENT + "labelled_tree.newick"
    params:
        placement_dir = PLAGIO_PLACEMENT,
        consensus_thresh = config["gappa_consensus_thresh"]
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}
        cp {input.clades} {params.placement_dir}

        cd {params.placement_dir}

        gappa examine assign \
            --jplace-path $(basename {input.jplace}) \
            --consensus-thresh {params.consensus_thresh} \
            --taxon-file $(basename {input.clades}) \
            --best-hit \
            --per-query-results \
            --allow-file-overwriting

        """


rule gappa_lwr_edpl_plagio:
    input:
        jplace = PLAGIO_PLACEMENT + "epa_result.jplace"
    output:
        lwr_histogram  = PLAGIO_PLACEMENT + "lwr-histogram.csv",
        lwr_list       = PLAGIO_PLACEMENT + "lwr-list.csv",
        edpl_histogram = PLAGIO_PLACEMENT + "edpl_histogram.csv",
        edpl_list      = PLAGIO_PLACEMENT + "edpl_list.csv"
    params:
        placement_dir = PLAGIO_PLACEMENT,
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}

        cd {params.placement_dir}

        gappa examine lwr-histogram \
            --jplace-path $(basename {input.jplace}) \
            --allow-file-overwriting

        gappa examine lwr-list \
            --jplace-path $(basename {input.jplace}) \
            --allow-file-overwriting

        gappa examine edpl \
            --jplace-path $(basename {input.jplace}) \
            --allow-file-overwriting

        """


# ---- PLAGIOPYLEA PLACEMENT FILTERING (replicates filter_euk_placements / extract_unassigned_ciliates
#      on the Plagiopylea tree's own placement results) ----

rule filter_plagio_final:
    input:
        per_query = PLAGIO_PLACEMENT + "per_query.tsv",
        lwr_list  = PLAGIO_PLACEMENT + "lwr-list.csv",
        edpl_list = PLAGIO_PLACEMENT + "edpl_list.csv"
    output:
        tsv = PLAGIO_PLACEMENT + "filtered_Plagio_LWR_EDPL.tsv"
    params:
        edpl_threshold = config["edpl_threshold"]
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/filter_placements.py"


rule extract_plagio_final_seqs:
    input:
        tsv   = PLAGIO_PLACEMENT + "filtered_Plagio_LWR_EDPL.tsv",
        fasta = PLAGIO_PLACEMENT + "plagio_seqs.fasta"
    output:
        fasta = PLAGIO_PLACEMENT + "plagio_placed_seqs.fasta"
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/extract_fasta.py"


# ---- SCUTI PLACEMENT ----

rule extract_scuti_seqs:
    input:
        tsv   = CIL_PLACEMENT + "filtered_Scuti_LWR_EDPL.tsv",
        fasta = CIL_PLACEMENT + "ciliate_seqs.fasta"
    output:
        fasta = SCUTI_PLACEMENT + "scuti_seqs.fasta"
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/extract_fasta.py"


rule papara_scuti:
    input:
        query_seqs = SCUTI_PLACEMENT + "scuti_seqs.fasta",
        tree       = config["SCUTI_TREE"],
        msa_phylip = config["MSA_PHYLIP_SCUTI"],
        msa_fasta  = config["MSA_FASTA_SCUTI"],
        clades     = config["CLADES_SCUTI"]
    output:
        papara_alignment = SCUTI_PLACEMENT + "papara_alignment.default"
    params:
        placement_dir = SCUTI_PLACEMENT
    threads: 17
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}

        cp {input.tree}       {params.placement_dir}
        cp {input.msa_phylip} {params.placement_dir}
        cp {input.msa_fasta}  {params.placement_dir}
        cp {input.clades}     {params.placement_dir}

        cd {params.placement_dir}

        module load papara_nt/2.5
        papara \
            -t $(basename {input.tree}) \
            -s $(basename {input.msa_phylip}) \
            -q $(basename {input.query_seqs}) \
            -j {threads} \
            -r

        """


rule epa_split_scuti:
    input:
        papara_alignment = SCUTI_PLACEMENT + "papara_alignment.default",
        msa_fasta        = config["MSA_FASTA_SCUTI"]
    output:
        query_fasta = SCUTI_PLACEMENT + "query.fasta",
        ref_fasta   = SCUTI_PLACEMENT + "reference.fasta"
    params:
        placement_dir = SCUTI_PLACEMENT
    threads: 17
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}
        cp {input.msa_fasta}        {params.placement_dir}

        cd {params.placement_dir}

        epa-ng \
            --split $(basename {input.msa_fasta}) \
            $(basename {input.papara_alignment})

        """


rule raxml_evaluate_scuti:
    input:
        ref_fasta = SCUTI_PLACEMENT + "reference.fasta",
        tree      = config["SCUTI_TREE"]
    output:
        best_model = SCUTI_PLACEMENT + "reference.fasta.raxml.bestModel",
        best_tree  = SCUTI_PLACEMENT + "reference.fasta.raxml.bestTree",
        log        = SCUTI_PLACEMENT + "reference.fasta.raxml.log",
        rba        = SCUTI_PLACEMENT + "reference.fasta.raxml.rba",
        start_tree = SCUTI_PLACEMENT + "reference.fasta.raxml.startTree"
    params:
        placement_dir = SCUTI_PLACEMENT
    threads: 17
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}
        cp {input.tree}      {params.placement_dir}

        cd {params.placement_dir}

        raxml-ng \
            --evaluate \
            --msa $(basename {input.ref_fasta}) \
            --tree $(basename {input.tree}) \
            --model GTR+G \
            --threads {threads} --force perf_threads

        """


rule epa_placement_scuti:
    input:
        tree        = config["SCUTI_TREE"],
        ref_fasta   = SCUTI_PLACEMENT + "reference.fasta",
        query_fasta = SCUTI_PLACEMENT + "query.fasta",
        best_model  = SCUTI_PLACEMENT + "reference.fasta.raxml.bestModel"
    output:
        jplace   = SCUTI_PLACEMENT + "epa_result.jplace",
        info_log = SCUTI_PLACEMENT + "epa_info.log"
    params:
        placement_dir = SCUTI_PLACEMENT,
        filter_acc_lwr = config["epa_filter_acc_lwr"],
        filter_max     = config["epa_filter_max"]
    threads: 17
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}
        cp {input.tree}        {params.placement_dir}

        rm -f {params.placement_dir}epa_result.jplace {params.placement_dir}epa_info.log

        cd {params.placement_dir}

        epa-ng \
            --filter-acc-lwr {params.filter_acc_lwr} \
            --filter-max {params.filter_max} \
            -t $(basename {input.tree}) \
            -s $(basename {input.ref_fasta}) \
            -q $(basename {input.query_fasta}) \
            --model $(basename {input.best_model})

        """


rule gappa_heat_tree_scuti:
    input:
        jplace = SCUTI_PLACEMENT + "epa_result.jplace"
    output:
        svg    = SCUTI_PLACEMENT + "tree.svg",
        newick = SCUTI_PLACEMENT + "tree.newick",
        nexus  = SCUTI_PLACEMENT + "tree.nexus"
    params:
        placement_dir = SCUTI_PLACEMENT,
        mass_norm   = config["gappa_mass_norm"]
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}

        cd {params.placement_dir}

        gappa examine heat-tree \
            --jplace-path $(basename {input.jplace}) \
            --mass-norm {params.mass_norm} \
            --allow-file-overwriting \
            --write-svg-tree \
            --write-newick-tree \
            --write-nexus-tree

        """


rule gappa_assign_scuti:
    input:
        jplace = SCUTI_PLACEMENT + "epa_result.jplace",
        clades = config["CLADES_SCUTI"]
    output:
        per_query     = SCUTI_PLACEMENT + "per_query.tsv",
        profile       = SCUTI_PLACEMENT + "profile.tsv",
        labelled_tree = SCUTI_PLACEMENT + "labelled_tree.newick"
    params:
        placement_dir = SCUTI_PLACEMENT,
        consensus_thresh = config["gappa_consensus_thresh"]
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}
        cp {input.clades} {params.placement_dir}

        cd {params.placement_dir}

        gappa examine assign \
            --jplace-path $(basename {input.jplace}) \
            --consensus-thresh {params.consensus_thresh} \
            --taxon-file $(basename {input.clades}) \
            --best-hit \
            --per-query-results \
            --allow-file-overwriting

        """


rule gappa_lwr_edpl_scuti:
    input:
        jplace = SCUTI_PLACEMENT + "epa_result.jplace"
    output:
        lwr_histogram  = SCUTI_PLACEMENT + "lwr-histogram.csv",
        lwr_list       = SCUTI_PLACEMENT + "lwr-list.csv",
        edpl_histogram = SCUTI_PLACEMENT + "edpl_histogram.csv",
        edpl_list      = SCUTI_PLACEMENT + "edpl_list.csv"
    params:
        placement_dir = SCUTI_PLACEMENT,
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}

        cd {params.placement_dir}

        gappa examine lwr-histogram \
            --jplace-path $(basename {input.jplace}) \
            --allow-file-overwriting

        gappa examine lwr-list \
            --jplace-path $(basename {input.jplace}) \
            --allow-file-overwriting

        gappa examine edpl \
            --jplace-path $(basename {input.jplace}) \
            --allow-file-overwriting

        """


# ---- SCUTI PLACEMENT FILTERING (replicates filter_euk_placements / extract_unassigned_ciliates
#      on the Scuti tree's own placement results) ----

rule filter_scuti_final:
    input:
        per_query = SCUTI_PLACEMENT + "per_query.tsv",
        lwr_list  = SCUTI_PLACEMENT + "lwr-list.csv",
        edpl_list = SCUTI_PLACEMENT + "edpl_list.csv"
    output:
        tsv = SCUTI_PLACEMENT + "filtered_Scuti_LWR_EDPL.tsv"
    params:
        edpl_threshold = config["edpl_threshold"]
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/filter_placements.py"


rule extract_scuti_final_seqs:
    input:
        tsv   = SCUTI_PLACEMENT + "filtered_Scuti_LWR_EDPL.tsv",
        fasta = SCUTI_PLACEMENT + "scuti_seqs.fasta"
    output:
        fasta = SCUTI_PLACEMENT + "scuti_placed_seqs.fasta"
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/extract_fasta.py"


# ---- R FIGURE GENERATION ----
# Like every other rule, R and its packages (phyloseq, sf, rnaturalearth,
# etc.) are provisioned automatically via conda from envs/r_figures.yaml.

rule prepare_phyloseq_objects:
    input:
        metadata = "documents/merged/merged_metadata.csv",
        envo_map = "reference/ENVO_IDs.csv",
        apm = APM_PLACEMENT + "filtered_APM_LWR_EDPL.tsv",
        plagio = PLAGIO_PLACEMENT + "filtered_Plagio_LWR_EDPL.tsv",
        scuti = SCUTI_PLACEMENT + "filtered_Scuti_LWR_EDPL.tsv",
        cil_counts = MERGED + "export/table/merged-ciliophora-table.tsv",
        unassigned_counts = MERGED + "export/table/merged-unassigned-table.tsv"
    output:
        taxonomy = FIGURES + "pp_taxonomy.tsv",
        counts = FIGURES + "combined_count_table.tsv",
        meta_clean = FIGURES + "metadata_cleaned.csv",
        meta_rds = FIGURES + "metadata_sample_data.rds",
        tax_rds = FIGURES + "tax_table.rds",
        otu_rds = FIGURES + "otu_table.rds",
        ps_rds = FIGURES + "ps.rds",
        missing_report = FIGURES + "metadata_gaps.txt"
    conda:
        "envs/r_figures.yaml"
    shell:
        """
        mkdir -p {FIGURES}
        Rscript scripts/prepare_phyloseq_objects.R \
            {input.metadata} \
            {input.apm} {input.plagio} {input.scuti} \
            {input.cil_counts} {input.unassigned_counts} \
            {output.taxonomy} {output.counts} {output.meta_clean} \
            {output.meta_rds} {output.tax_rds} {output.otu_rds} {output.ps_rds} \
            {input.envo_map} {output.missing_report}
        """


rule generate_figures:
    input:
        ps_rds = FIGURES + "ps.rds",
        meta_clean = FIGURES + "metadata_cleaned.csv"
    output:
        bar_pdf = FIGURES + "taxa_barplot.pdf",
        bar_png = FIGURES + "taxa_barplot.png",
        map_pdf = FIGURES + "sample_map.pdf",
        map_png = FIGURES + "sample_map.png"
    conda:
        "envs/r_figures.yaml"
    shell:
        """
        Rscript scripts/generate_figures.R \
            {input.ps_rds} {input.meta_clean} \
            {output.bar_pdf} {output.bar_png} {output.map_pdf} {output.map_png}
        """