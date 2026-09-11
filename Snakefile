# Snakemake file - batch download and process multiple SRA studies

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
SRA_CACHE = config["sra_cache"]
FASTERQ_TEMP = config["fasterq_temp"]
QIIME_ENV = config["qiime_environment"]

# Merged paths
MERGED = OUTPUT + "_merged/"

# Placements
PLACEMENT        = OUTPUT + "_placement/"
TAXON_PLACEMENT    = OUTPUT + "_Taxon_placement/"
CLADE1_PLACEMENT    = OUTPUT + "_Clade1_placement/"
CLADE2_PLACEMENT = OUTPUT + "_Clade2_placement/"
CLADE3_PLACEMENT  = OUTPUT + "_Clade3_placement/"

# R figure generation (phyloseq objects, taxa barplots, sample map)
FIGURES = OUTPUT + "_figures/"

# ---- PAPARA CONFIGURATION ----
PAPARA_MODULE = config.get("papara_module")
PAPARA_SETUP = (
    f"module load {PAPARA_MODULE}"
    if PAPARA_MODULE
    else ":"
)
PAPARA_EXECUTABLE = config.get("papara_executable", "papara")
# ---- END PAPARA CONFIGURATION ----

# ---- CLADE TREE CONFIGURATION ----
MAX_CLADE = int(config.get("max_clade", 3))
if MAX_CLADE not in (1, 2, 3):
    raise ValueError("max_clade must be 1, 2, or 3")

CLADE_DIRS = [
    CLADE1_PLACEMENT,
    CLADE2_PLACEMENT,
    CLADE3_PLACEMENT,
]

SELECTED_CLADE_TARGETS = [
    target
    for number, clade_dir in enumerate(CLADE_DIRS[:MAX_CLADE], start=1)
    for target in (
        clade_dir + "tree.svg",
        clade_dir + f"Clade{number}_placed_seqs.fasta",
    )
]

# The current R script requires results from all three clades.
FIGURE_TARGETS = (
    [
        FIGURES + "taxa_barplot.pdf",
        FIGURES + "sample_map.pdf",
    ]
    if MAX_CLADE == 3
    else []
)
# ---- END CLADE TREE CONFIGURATION ----

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


# Targets for everything up to (but not including) phylogenetic placement.
# Shared by `rule all` and `rule preplacement` so the list only lives in one place.
PREPLACEMENT_TARGETS = dict(
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

    # Taxon merged artifacts and exports
    merged_Taxon_table = MERGED + "merged-Taxon-table.qza",
    merged_Taxon_seqs = MERGED + "merged-Taxon-seqs.qza",
    table_Taxon_biom = MERGED + "export/table/merged-Taxon-table.biom",
    table_Taxon_tsv = MERGED + "export/table/merged-Taxon-table.tsv",
    rep_seqs_Taxon_fasta = MERGED + "export/merged-Taxon-seqs.fasta",

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
)


rule all:
    input:
        **PREPLACEMENT_TARGETS,

        # Phylogenetic placement outputs
        # changed these to _unassigned
        unassigned_epa_jplace = PLACEMENT + "epa_result.jplace",
        unassigned_heat_tree_svg = PLACEMENT + "tree.svg",
        unassigned_per_query = PLACEMENT + "per_query.tsv",
        unassigned_lwr_histogram = PLACEMENT + "lwr-histogram.csv",
        unassigned_edpl_histogram = PLACEMENT + "edpl_histogram.csv",

        # EUK placement filtering
        per_query_Taxon = PLACEMENT + "per_query_Taxon.tsv",
        filtered_euk_placements = PLACEMENT + "filtered_Taxon_LWR_EDPL.tsv",
        unassigned_Taxon_fasta = PLACEMENT + "unassigned_Taxon.fasta",

        # Taxa placement outputs
        taxa_jplace = TAXON_PLACEMENT + "epa_result.jplace",
        taxa_heat_tree_svg = TAXON_PLACEMENT + "tree.svg",
        taxa_per_query = TAXON_PLACEMENT + "per_query.tsv",
        taxa_lwr_histogram = TAXON_PLACEMENT + "lwr-histogram.csv",
        taxa_edpl_histogram = TAXON_PLACEMENT + "edpl_histogram.csv",

        # Clade specification and figure gen
        selected_clades = SELECTED_CLADE_TARGETS,
        figures = FIGURE_TARGETS


rule preplacement:
    # Run everything up to (but not including) phylogenetic placement:
    #   snakemake preplacement --cores N --use-conda
    input:
        **PREPLACEMENT_TARGETS,
rule download_runinfo:
    output:
        runinfo = RAW_DATA + "{project}/runinfo.csv"
    conda:
        "envs/sra_download.yaml"
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
        cut -d ',' -f 1 {input.runinfo} | tail -n +2 > {output.SRRnumbers}
        """


rule fasterq_dump:
    input:
        SRRnumbers = RAW_DATA + "{project}/SRR.numbers"
    output:
        rawreads = directory(RAW_DATA + "{project}/reads")
    params:
        sra_dir = SRA_CACHE + "{project}",
        temp_dir = FASTERQ_TEMP + "{project}"
    log:
        RAW_DATA + "{project}/logs/fasterq_dump.log"
    threads: 10
    conda:
        "envs/sra_download.yaml"
    shell:
        """
        mkdir -p \
            "{output.rawreads}" \
            "{params.sra_dir}" \
            "{params.temp_dir}" \
            "$(dirname "{log}")"

        (
            cr=$(printf '\\r')

            tr -d "$cr" < "{input.SRRnumbers}" |
            while IFS= read -r srr; do
                [ -n "$srr" ] || continue

                echo "Downloading $srr"

                prefetch "$srr" \
                    --output-directory "{params.sra_dir}"

                fasterq-dump \
                    "{params.sra_dir}/$srr/$srr.sra" \
                    --threads {threads} \
                    --split-files \
                    --outdir "{output.rawreads}" \
                    --temp "{params.temp_dir}"

                for fastq in "{output.rawreads}/$srr"*.fastq; do
                    [ -e "$fastq" ] || continue
                    pigz -p {threads} "$fastq"
                done
            done
        ) 2>&1 | tee "{log}"
        """


rule rename_files_import:
    input:
        rawreads = RAW_DATA + "{project}/reads"
    output:
        q2_import = OUTPUT + "{project}/{project}-PE-demux.qza"
    log:
        OUTPUT + "{project}/logs/{project}_q2.log"
    conda:
        QIIME_ENV
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
        QIIME_ENV
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
        QIIME_ENV
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
        QIIME_ENV
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
        QIIME_ENV
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
        QIIME_ENV
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
        QIIME_ENV
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
        QIIME_ENV
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
        QIIME_ENV
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
        QIIME_ENV
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
        QIIME_ENV
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
        tables_Taxon = expand(OUTPUT + "{project}/{project}-table-taxa.qza", project=PROJECTS),
        seqs_Taxon = expand(OUTPUT + "{project}/{project}-rep-seqs-taxa.qza", project=PROJECTS),
        tables_unassigned = expand(OUTPUT + "{project}/{project}-table-unassigned.qza", project=PROJECTS),
        seqs_unassigned = expand(OUTPUT + "{project}/{project}-rep-seqs-unassigned.qza", project=PROJECTS)
    output:
        merged_Taxon_table = MERGED + "merged-Taxon-table.qza",
        merged_Taxon_seqs = MERGED + "merged-Taxon-seqs.qza",
        merged_unassigned_table = MERGED + "merged-unassigned-table.qza",
        merged_unassigned_seqs = MERGED + "merged-unassigned-seqs.qza"
    log:
        MERGED + "logs/merge_filtered_q2.log"
    conda:
        QIIME_ENV
    shell:
        """
        mkdir -p $(dirname {output.merged_Taxon_table})
        mkdir -p $(dirname {log})
        qiime feature-table merge \
            --i-tables {input.tables_Taxon} \
            --o-merged-table {output.merged_Taxon_table}
        qiime feature-table merge-seqs \
            --i-data {input.seqs_Taxon} \
            --o-merged-data {output.merged_Taxon_seqs}
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
        QIIME_ENV
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


rule export_Taxon:
    input:
        table_Taxon_merged = MERGED + "merged-Taxon-table.qza",
        rep_Taxon_merged = MERGED + "merged-Taxon-seqs.qza"
    output:
        table_Taxon_biom = MERGED + "export/table/merged-Taxon-table.biom",
        table_Taxon_tsv = MERGED + "export/table/merged-Taxon-table.tsv",
        rep_seqs_Taxon_fasta = MERGED + "export/merged-Taxon-seqs.fasta"
    log:
        MERGED + "logs/export_Taxon_q2.log"
    conda:
        QIIME_ENV
    params:
        table_outdir = MERGED + "export/table/",
        seqs_outdir = MERGED + "export/",
        # private temp dirs so a concurrently-running export_unassigned job can't
        # clobber our intermediate feature-table.biom / dna-sequences.fasta
        table_tmpdir = MERGED + "export/table/.tmp_Taxon/",
        seqs_tmpdir = MERGED + "export/.tmp_Taxon_seqs/"
    shell:
        """
        mkdir -p {params.table_outdir}
        mkdir -p {params.seqs_outdir}
        mkdir -p {params.table_tmpdir}
        mkdir -p {params.seqs_tmpdir}
        mkdir -p $(dirname {log})
        qiime tools export --input-path {input.table_Taxon_merged} --output-path {params.table_tmpdir}
        mv {params.table_tmpdir}feature-table.biom {output.table_Taxon_biom}
        rmdir {params.table_tmpdir}
        biom convert -i {output.table_Taxon_biom} -o {output.table_Taxon_tsv} --to-tsv
        qiime tools export --input-path {input.rep_Taxon_merged} --output-path {params.seqs_tmpdir}
        mv {params.seqs_tmpdir}dna-sequences.fasta {output.rep_seqs_Taxon_fasta}
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
        QIIME_ENV
    params:
        table_outdir = MERGED + "export/table/",
        seqs_outdir = MERGED + "export/",
        # private temp dirs so a concurrently-running export_Taxon job can't
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
        QIIME_ENV
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
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/merge_metadata.py"

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
        placement_dir = PLACEMENT,
        papara_setup = PAPARA_SETUP,
        papara_executable = PAPARA_EXECUTABLE
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

        {params.papara_setup}
        {params.papara_executable} \
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

rule filter_Taxon:
    input:
        tsv = PLACEMENT + "per_query.tsv"
    output:
        tsv = PLACEMENT + "per_query_Taxon.tsv"
    params:
        search_term = config["search_term_taxa"]
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/filter_by_taxopath.py"


rule filter_euk_placements:
    input:
        per_query = PLACEMENT + "per_query_Taxon.tsv",
        lwr_list  = PLACEMENT + "lwr-list.csv",
        edpl_list = PLACEMENT + "edpl_list.csv"
    output:
        tsv = PLACEMENT + "filtered_Taxon_LWR_EDPL.tsv"
    params:
        edpl_threshold = config["edpl_threshold"]
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/filter_placements.py"


rule extract_unassigned_Taxon:
    input:
        tsv   = PLACEMENT + "filtered_Taxon_LWR_EDPL.tsv",
        fasta = MERGED + "export/merged-unassigned-seqs.fasta"
    output:
        fasta = PLACEMENT + "unassigned_Taxon.fasta"
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/extract_fasta.py"


# ---- TAXON PLACEMENT ----

rule combine_Taxon_seqs:
    input:
        unassigned_Taxon = PLACEMENT + "unassigned_Taxon.fasta",
        asvs_Taxon       = MERGED + "export/merged-Taxon-seqs.fasta"
    output:
        combined = TAXON_PLACEMENT + "Taxon_seqs.fasta"
    shell:
        """
        mkdir -p {TAXON_PLACEMENT}
        cat {input.unassigned_Taxon} {input.asvs_Taxon} > {output.combined}
        """


rule papara_Taxon:
    input:
        Taxon_seqs = TAXON_PLACEMENT + "Taxon_seqs.fasta",
        Taxon_tree     = config["TAXON_TREE"],
        msa_phylip   = config["MSA_PHYLIP_TAXON"],
        msa_fasta    = config["MSA_FASTA_TAXON"],
        clades       = config["CLADES_TAXON"]
    output:
        papara_alignment = TAXON_PLACEMENT + "papara_alignment.default"
    params:
        placement_dir = TAXON_PLACEMENT,
        papara_setup = PAPARA_SETUP,
        papara_executable = PAPARA_EXECUTABLE
    threads: 17
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}

        cp {input.Taxon_tree}   {params.placement_dir}
        cp {input.msa_phylip} {params.placement_dir}
        cp {input.msa_fasta}  {params.placement_dir}
        cp {input.clades}     {params.placement_dir}

        cd {params.placement_dir}

        {params.papara_setup}
        {params.papara_executable} \
            -t $(basename {input.Taxon_tree}) \
            -s $(basename {input.msa_phylip}) \
            -q $(basename {input.Taxon_seqs}) \
            -j {threads} \
            -r

        """


rule epa_split_Taxon:
    input:
        papara_alignment = TAXON_PLACEMENT + "papara_alignment.default",
        msa_fasta        = config["MSA_FASTA_TAXON"]
    output:
        query_fasta = TAXON_PLACEMENT + "query.fasta",
        ref_fasta   = TAXON_PLACEMENT + "reference.fasta"
    params:
        placement_dir = TAXON_PLACEMENT
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


rule raxml_evaluate_Taxon:
    input:
        ref_fasta = TAXON_PLACEMENT + "reference.fasta",
        Taxon_tree  = config["TAXON_TREE"]
    output:
        best_model = TAXON_PLACEMENT + "reference.fasta.raxml.bestModel",
        best_tree  = TAXON_PLACEMENT + "reference.fasta.raxml.bestTree",
        log        = TAXON_PLACEMENT + "reference.fasta.raxml.log",
        rba        = TAXON_PLACEMENT + "reference.fasta.raxml.rba",
        start_tree = TAXON_PLACEMENT + "reference.fasta.raxml.startTree"
    params:
        placement_dir = TAXON_PLACEMENT
    threads: 17
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}
        cp {input.Taxon_tree}  {params.placement_dir}

        cd {params.placement_dir}

        raxml-ng \
            --evaluate \
            --msa $(basename {input.ref_fasta}) \
            --tree $(basename {input.Taxon_tree}) \
            --model GTR+G \
            --threads {threads} --force perf_threads

        """


rule epa_placement_Taxon:
    input:
        Taxon_tree   = config["TAXON_TREE"],
        ref_fasta  = TAXON_PLACEMENT + "reference.fasta",
        query_fasta = TAXON_PLACEMENT + "query.fasta",
        best_model = TAXON_PLACEMENT + "reference.fasta.raxml.bestModel"
    output:
        jplace   = TAXON_PLACEMENT + "epa_result.jplace",
        info_log = TAXON_PLACEMENT + "epa_info.log"
    params:
        placement_dir   = TAXON_PLACEMENT,
        filter_acc_lwr  = config["epa_filter_acc_lwr"],
        filter_max      = config["epa_filter_max"]
    threads: 17
    conda:
        "envs/phylo_placement.yaml"
    shell:
        """
        mkdir -p {params.placement_dir}
        cp {input.Taxon_tree}    {params.placement_dir}

        rm -f {params.placement_dir}epa_result.jplace \
              {params.placement_dir}epa_info.log

        cd {params.placement_dir}

        epa-ng \
            --filter-acc-lwr {params.filter_acc_lwr} \
            --filter-max {params.filter_max} \
            -t $(basename {input.Taxon_tree}) \
            -s $(basename {input.ref_fasta}) \
            -q $(basename {input.query_fasta}) \
            --model $(basename {input.best_model})

        """


rule gappa_heat_tree_Taxon:
    input:
        jplace = TAXON_PLACEMENT + "epa_result.jplace"
    output:
        svg    = TAXON_PLACEMENT + "tree.svg",
        newick = TAXON_PLACEMENT + "tree.newick",
        nexus  = TAXON_PLACEMENT + "tree.nexus"
    params:
        placement_dir = TAXON_PLACEMENT,
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


rule gappa_assign_Taxon:
    input:
        jplace = TAXON_PLACEMENT + "epa_result.jplace",
        clades = config["CLADES_TAXON"]
    output:
        per_query    = TAXON_PLACEMENT + "per_query.tsv",
        profile      = TAXON_PLACEMENT + "profile.tsv",
        labelled_tree = TAXON_PLACEMENT + "labelled_tree.newick"
    params:
        placement_dir = TAXON_PLACEMENT,
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


rule gappa_lwr_edpl_Taxon:
    input:
        jplace = TAXON_PLACEMENT + "epa_result.jplace"
    output:
        lwr_histogram  = TAXON_PLACEMENT + "lwr-histogram.csv",
        lwr_list       = TAXON_PLACEMENT + "lwr-list.csv",
        edpl_histogram = TAXON_PLACEMENT + "edpl_histogram.csv",
        edpl_list      = TAXON_PLACEMENT + "edpl_list.csv"
    params:
        placement_dir = TAXON_PLACEMENT,
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
        tsv = TAXON_PLACEMENT + "per_query.tsv"
    output:
        tsv = TAXON_PLACEMENT + "per_query_Clade1.tsv"
    params:
        search_term = config["search_term_Clade1"]
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/filter_by_taxopath.py"


rule filter_Clade2:
    input:
        tsv = TAXON_PLACEMENT + "per_query.tsv"
    output:
        tsv = TAXON_PLACEMENT + "per_query_Clade2.tsv"
    params:
        search_term = config["search_term_Clade2"]
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/filter_by_taxopath.py"


rule filter_Clade3:
    input:
        tsv = TAXON_PLACEMENT + "per_query.tsv"
    output:
        tsv = TAXON_PLACEMENT + "per_query_Clade3.tsv"
    params:
        search_term = config["search_term_Clade3"]
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/filter_by_taxopath.py"


rule filter_CLADE1_placements:
    input:
        per_query = TAXON_PLACEMENT + "per_query_Clade1.tsv",
        lwr_list  = TAXON_PLACEMENT + "lwr-list.csv",
        edpl_list = TAXON_PLACEMENT + "edpl_list.csv"
    output:
        tsv = TAXON_PLACEMENT + "filtered_Clade1_LWR_EDPL.tsv"
    params:
        edpl_threshold = config["edpl_threshold"]
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/filter_placements.py"


rule filter_Clade2_placements:
    input:
        per_query = TAXON_PLACEMENT + "per_query_Clade2.tsv",
        lwr_list  = TAXON_PLACEMENT + "lwr-list.csv",
        edpl_list = TAXON_PLACEMENT + "edpl_list.csv"
    output:
        tsv = TAXON_PLACEMENT + "filtered_Clade2_LWR_EDPL.tsv"
    params:
        edpl_threshold = config["edpl_threshold"]
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/filter_placements.py"


rule filter_Clade3_placements:
    input:
        per_query = TAXON_PLACEMENT + "per_query_Clade3.tsv",
        lwr_list  = TAXON_PLACEMENT + "lwr-list.csv",
        edpl_list = TAXON_PLACEMENT + "edpl_list.csv"
    output:
        tsv = TAXON_PLACEMENT + "filtered_Clade3_LWR_EDPL.tsv"
    params:
        edpl_threshold = config["edpl_threshold"]
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/filter_placements.py"


# ---- CLADE1 PLACEMENT ----

rule extract_Clade1_seqs:
    input:
        tsv   = TAXON_PLACEMENT + "filtered_Clade1_LWR_EDPL.tsv",
        fasta = TAXON_PLACEMENT + "Taxon_seqs.fasta"
    output:
        fasta = CLADE1_PLACEMENT + "Clade1_seqs.fasta"
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/extract_fasta.py"


rule papara_Clade1:
    input:
        query_seqs = CLADE1_PLACEMENT + "Clade1_seqs.fasta",
        tree       = config["CLADE1_TREE"],
        msa_phylip = config["MSA_PHYLIP_CLADE1"],
        msa_fasta  = config["MSA_FASTA_CLADE1"],
        clades     = config["CLADES_CLADE1"]
    output:
        papara_alignment = CLADE1_PLACEMENT + "papara_alignment.default"
    params:
        placement_dir = CLADE1_PLACEMENT,
        papara_setup = PAPARA_SETUP,
        papara_executable = PAPARA_EXECUTABLE
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

        {params.papara_setup}
        {params.papara_executable} \
            -t $(basename {input.tree}) \
            -s $(basename {input.msa_phylip}) \
            -q $(basename {input.query_seqs}) \
            -j {threads} \
            -r

        """


rule epa_split_Clade1:
    input:
        papara_alignment = CLADE1_PLACEMENT + "papara_alignment.default",
        msa_fasta        = config["MSA_FASTA_CLADE1"]
    output:
        query_fasta = CLADE1_PLACEMENT + "query.fasta",
        ref_fasta   = CLADE1_PLACEMENT + "reference.fasta"
    params:
        placement_dir = CLADE1_PLACEMENT
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


rule raxml_evaluate_Clade1:
    input:
        ref_fasta = CLADE1_PLACEMENT + "reference.fasta",
        tree      = config["CLADE1_TREE"]
    output:
        best_model = CLADE1_PLACEMENT + "reference.fasta.raxml.bestModel",
        best_tree  = CLADE1_PLACEMENT + "reference.fasta.raxml.bestTree",
        log        = CLADE1_PLACEMENT + "reference.fasta.raxml.log",
        rba        = CLADE1_PLACEMENT + "reference.fasta.raxml.rba",
        start_tree = CLADE1_PLACEMENT + "reference.fasta.raxml.startTree"
    params:
        placement_dir = CLADE1_PLACEMENT
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


rule epa_placement_Clade1:
    input:
        tree        = config["CLADE1_TREE"],
        ref_fasta   = CLADE1_PLACEMENT + "reference.fasta",
        query_fasta = CLADE1_PLACEMENT + "query.fasta",
        best_model  = CLADE1_PLACEMENT + "reference.fasta.raxml.bestModel"
    output:
        jplace   = CLADE1_PLACEMENT + "epa_result.jplace",
        info_log = CLADE1_PLACEMENT + "epa_info.log"
    params:
        placement_dir = CLADE1_PLACEMENT,
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


rule gappa_heat_tree_Clade1:
    input:
        jplace = CLADE1_PLACEMENT + "epa_result.jplace"
    output:
        svg    = CLADE1_PLACEMENT + "tree.svg",
        newick = CLADE1_PLACEMENT + "tree.newick",
        nexus  = CLADE1_PLACEMENT + "tree.nexus"
    params:
        placement_dir = CLADE1_PLACEMENT,
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


rule gappa_assign_Clade1:
    input:
        jplace = CLADE1_PLACEMENT + "epa_result.jplace",
        clades = config["CLADES_CLADE1"]
    output:
        per_query     = CLADE1_PLACEMENT + "per_query.tsv",
        profile       = CLADE1_PLACEMENT + "profile.tsv",
        labelled_tree = CLADE1_PLACEMENT + "labelled_tree.newick"
    params:
        placement_dir = CLADE1_PLACEMENT,
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


rule gappa_lwr_edpl_Clade1:
    input:
        jplace = CLADE1_PLACEMENT + "epa_result.jplace"
    output:
        lwr_histogram  = CLADE1_PLACEMENT + "lwr-histogram.csv",
        lwr_list       = CLADE1_PLACEMENT + "lwr-list.csv",
        edpl_histogram = CLADE1_PLACEMENT + "edpl_histogram.csv",
        edpl_list      = CLADE1_PLACEMENT + "edpl_list.csv"
    params:
        placement_dir = CLADE1_PLACEMENT,
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


# ---- CLADE1 PLACEMENT FILTERING

rule filter_Clade1_final:
    input:
        per_query = CLADE1_PLACEMENT + "per_query.tsv",
        lwr_list  = CLADE1_PLACEMENT + "lwr-list.csv",
        edpl_list = CLADE1_PLACEMENT + "edpl_list.csv"
    output:
        tsv = CLADE1_PLACEMENT + "filtered_Clade1_LWR_EDPL.tsv"
    params:
        edpl_threshold = config["edpl_threshold"]
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/filter_placements.py"


rule extract_Clade1_final_seqs:
    input:
        tsv   = CLADE1_PLACEMENT + "filtered_Clade1_LWR_EDPL.tsv",
        fasta = CLADE1_PLACEMENT + "Clade1_seqs.fasta"
    output:
        fasta = CLADE1_PLACEMENT + "Clade1_placed_seqs.fasta"
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/extract_fasta.py"


# ---- CLADE 2 PLACEMENT ----

rule extract_Clade2_seqs:
    input:
        tsv   = TAXON_PLACEMENT + "filtered_Clade2_LWR_EDPL.tsv",
        fasta = TAXON_PLACEMENT + "Taxon_seqs.fasta"
    output:
        fasta = CLADE2_PLACEMENT + "Clade2_seqs.fasta"
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/extract_fasta.py"


rule papara_Clade2:
    input:
        query_seqs = CLADE2_PLACEMENT + "Clade2_seqs.fasta",
        tree       = config["CLADE2_TREE"],
        msa_phylip = config["MSA_PHYLIP_CLADE2"],
        msa_fasta  = config["MSA_FASTA_CLADE2"],
        clades     = config["CLADES_CLADE2"]
    output:
        papara_alignment = CLADE2_PLACEMENT + "papara_alignment.default"
    params:
        placement_dir = CLADE2_PLACEMENT,
        papara_setup = PAPARA_SETUP,
        papara_executable = PAPARA_EXECUTABLE
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

        {params.papara_setup}
        {params.papara_executable} \
            -t $(basename {input.tree}) \
            -s $(basename {input.msa_phylip}) \
            -q $(basename {input.query_seqs}) \
            -j {threads} \
            -r

        """


rule epa_split_Clade2:
    input:
        papara_alignment = CLADE2_PLACEMENT + "papara_alignment.default",
        msa_fasta        = config["MSA_FASTA_CLADE2"]
    output:
        query_fasta = CLADE2_PLACEMENT + "query.fasta",
        ref_fasta   = CLADE2_PLACEMENT + "reference.fasta"
    params:
        placement_dir = CLADE2_PLACEMENT
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


rule raxml_evaluate_Clade2:
    input:
        ref_fasta = CLADE2_PLACEMENT + "reference.fasta",
        tree      = config["CLADE2_TREE"]
    output:
        best_model = CLADE2_PLACEMENT + "reference.fasta.raxml.bestModel",
        best_tree  = CLADE2_PLACEMENT + "reference.fasta.raxml.bestTree",
        log        = CLADE2_PLACEMENT + "reference.fasta.raxml.log",
        rba        = CLADE2_PLACEMENT + "reference.fasta.raxml.rba",
        start_tree = CLADE2_PLACEMENT + "reference.fasta.raxml.startTree"
    params:
        placement_dir = CLADE2_PLACEMENT
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


rule epa_placement_Clade2:
    input:
        tree        = config["CLADE2_TREE"],
        ref_fasta   = CLADE2_PLACEMENT + "reference.fasta",
        query_fasta = CLADE2_PLACEMENT + "query.fasta",
        best_model  = CLADE2_PLACEMENT + "reference.fasta.raxml.bestModel"
    output:
        jplace   = CLADE2_PLACEMENT + "epa_result.jplace",
        info_log = CLADE2_PLACEMENT + "epa_info.log"
    params:
        placement_dir = CLADE2_PLACEMENT,
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


rule gappa_heat_tree_Clade2:
    input:
        jplace = CLADE2_PLACEMENT + "epa_result.jplace"
    output:
        svg    = CLADE2_PLACEMENT + "tree.svg",
        newick = CLADE2_PLACEMENT + "tree.newick",
        nexus  = CLADE2_PLACEMENT + "tree.nexus"
    params:
        placement_dir = CLADE2_PLACEMENT,
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


rule gappa_assign_Clade2:
    input:
        jplace = CLADE2_PLACEMENT + "epa_result.jplace",
        clades = config["CLADES_CLADE2"]
    output:
        per_query     = CLADE2_PLACEMENT + "per_query.tsv",
        profile       = CLADE2_PLACEMENT + "profile.tsv",
        labelled_tree = CLADE2_PLACEMENT + "labelled_tree.newick"
    params:
        placement_dir = CLADE2_PLACEMENT,
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


rule gappa_lwr_edpl_Clade2:
    input:
        jplace = CLADE2_PLACEMENT + "epa_result.jplace"
    output:
        lwr_histogram  = CLADE2_PLACEMENT + "lwr-histogram.csv",
        lwr_list       = CLADE2_PLACEMENT + "lwr-list.csv",
        edpl_histogram = CLADE2_PLACEMENT + "edpl_histogram.csv",
        edpl_list      = CLADE2_PLACEMENT + "edpl_list.csv"
    params:
        placement_dir = CLADE2_PLACEMENT,
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


# ---- CLADE 2 PLACEMENT FILTERING

rule filter_Clade2_final:
    input:
        per_query = CLADE2_PLACEMENT + "per_query.tsv",
        lwr_list  = CLADE2_PLACEMENT + "lwr-list.csv",
        edpl_list = CLADE2_PLACEMENT + "edpl_list.csv"
    output:
        tsv = CLADE2_PLACEMENT + "filtered_Clade2_LWR_EDPL.tsv"
    params:
        edpl_threshold = config["edpl_threshold"]
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/filter_placements.py"


rule extract_Clade2_final_seqs:
    input:
        tsv   = CLADE2_PLACEMENT + "filtered_Clade2_LWR_EDPL.tsv",
        fasta = CLADE2_PLACEMENT + "Clade2_seqs.fasta"
    output:
        fasta = CLADE2_PLACEMENT + "Clade2_placed_seqs.fasta"
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/extract_fasta.py"


# ---- CLADE 3 PLACEMENT ----

rule extract_Clade3_seqs:
    input:
        tsv   = TAXON_PLACEMENT + "filtered_Clade3_LWR_EDPL.tsv",
        fasta = TAXON_PLACEMENT + "Taxon_seqs.fasta"
    output:
        fasta = CLADE3_PLACEMENT + "Clade3_seqs.fasta"
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/extract_fasta.py"


rule papara_Clade3:
    input:
        query_seqs = CLADE3_PLACEMENT + "Clade3_seqs.fasta",
        tree       = config["CLADE3_TREE"],
        msa_phylip = config["MSA_PHYLIP_CLADE3"],
        msa_fasta  = config["MSA_FASTA_CLADE3"],
        clades     = config["CLADES_CLADE3"]
    output:
        papara_alignment = CLADE3_PLACEMENT + "papara_alignment.default"
    params:
        placement_dir = CLADE3_PLACEMENT,
        papara_setup = PAPARA_SETUP,
        papara_executable = PAPARA_EXECUTABLE
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

        {params.papara_setup}
        {params.papara_executable} \
            -t $(basename {input.tree}) \
            -s $(basename {input.msa_phylip}) \
            -q $(basename {input.query_seqs}) \
            -j {threads} \
            -r

        """


rule epa_split_Clade3:
    input:
        papara_alignment = CLADE3_PLACEMENT + "papara_alignment.default",
        msa_fasta        = config["MSA_FASTA_CLADE3"]
    output:
        query_fasta = CLADE3_PLACEMENT + "query.fasta",
        ref_fasta   = CLADE3_PLACEMENT + "reference.fasta"
    params:
        placement_dir = CLADE3_PLACEMENT
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


rule raxml_evaluate_Clade3:
    input:
        ref_fasta = CLADE3_PLACEMENT + "reference.fasta",
        tree      = config["CLADE3_TREE"]
    output:
        best_model = CLADE3_PLACEMENT + "reference.fasta.raxml.bestModel",
        best_tree  = CLADE3_PLACEMENT + "reference.fasta.raxml.bestTree",
        log        = CLADE3_PLACEMENT + "reference.fasta.raxml.log",
        rba        = CLADE3_PLACEMENT + "reference.fasta.raxml.rba",
        start_tree = CLADE3_PLACEMENT + "reference.fasta.raxml.startTree"
    params:
        placement_dir = CLADE3_PLACEMENT
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


rule epa_placement_Clade3:
    input:
        tree        = config["CLADE3_TREE"],
        ref_fasta   = CLADE3_PLACEMENT + "reference.fasta",
        query_fasta = CLADE3_PLACEMENT + "query.fasta",
        best_model  = CLADE3_PLACEMENT + "reference.fasta.raxml.bestModel"
    output:
        jplace   = CLADE3_PLACEMENT + "epa_result.jplace",
        info_log = CLADE3_PLACEMENT + "epa_info.log"
    params:
        placement_dir = CLADE3_PLACEMENT,
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


rule gappa_heat_tree_Clade3:
    input:
        jplace = CLADE3_PLACEMENT + "epa_result.jplace"
    output:
        svg    = CLADE3_PLACEMENT + "tree.svg",
        newick = CLADE3_PLACEMENT + "tree.newick",
        nexus  = CLADE3_PLACEMENT + "tree.nexus"
    params:
        placement_dir = CLADE3_PLACEMENT,
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


rule gappa_assign_Clade3:
    input:
        jplace = CLADE3_PLACEMENT + "epa_result.jplace",
        clades = config["CLADES_CLADE3"]
    output:
        per_query     = CLADE3_PLACEMENT + "per_query.tsv",
        profile       = CLADE3_PLACEMENT + "profile.tsv",
        labelled_tree = CLADE3_PLACEMENT + "labelled_tree.newick"
    params:
        placement_dir = CLADE3_PLACEMENT,
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


rule gappa_lwr_edpl_Clade3:
    input:
        jplace = CLADE3_PLACEMENT + "epa_result.jplace"
    output:
        lwr_histogram  = CLADE3_PLACEMENT + "lwr-histogram.csv",
        lwr_list       = CLADE3_PLACEMENT + "lwr-list.csv",
        edpl_histogram = CLADE3_PLACEMENT + "edpl_histogram.csv",
        edpl_list      = CLADE3_PLACEMENT + "edpl_list.csv"
    params:
        placement_dir = CLADE3_PLACEMENT,
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


# ---- CLADE 3 PLACEMENT FILTERING

rule filter_Clade3_final:
    input:
        per_query = CLADE3_PLACEMENT + "per_query.tsv",
        lwr_list  = CLADE3_PLACEMENT + "lwr-list.csv",
        edpl_list = CLADE3_PLACEMENT + "edpl_list.csv"
    output:
        tsv = CLADE3_PLACEMENT + "filtered_Clade3_LWR_EDPL.tsv"
    params:
        edpl_threshold = config["edpl_threshold"]
    conda:
        "envs/pysradb.yaml"
    script:
        "scripts/filter_placements.py"


rule extract_Clade3_final_seqs:
    input:
        tsv   = CLADE3_PLACEMENT + "filtered_Clade3_LWR_EDPL.tsv",
        fasta = CLADE3_PLACEMENT + "Clade3_seqs.fasta"
    output:
        fasta = CLADE3_PLACEMENT + "Clade3_placed_seqs.fasta"
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
        Clade1 = CLADE1_PLACEMENT + "filtered_Clade1_LWR_EDPL.tsv",
        Clade2 = CLADE2_PLACEMENT + "filtered_Clade2_LWR_EDPL.tsv",
        Clade3 = CLADE3_PLACEMENT + "filtered_Clade3_LWR_EDPL.tsv",
        Taxon_counts = MERGED + "export/table/merged-Taxon-table.tsv",
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
            {input.Clade1} {input.Clade2} {input.Clade3} \
            {input.Taxon_counts} {input.unassigned_counts} \
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
