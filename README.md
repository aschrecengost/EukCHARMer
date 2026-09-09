# EukCHARMer

### **Euk**aryote-**C**entric **H**arvesting and **A**nalysis of **R**ibosomal **M**ark**er** genes

Written by Anna Schrecengost and Jaliyah Harrison, with help from this QIIME2 Snakemake tutorial from Sarah Hu (1) and this paper from Isabelle Ewers et al. (2). 

Contact us via email with questions: aschrecengost@uri.edu; jaliyahdharrison@gmail.com

## Table of Contents
- [Description](#Description)
- [Set-up](#Setup)
- [Tutorial](#Getting) 


<img width="2978" height="2284" alt="Pipeline" src="https://github.com/user-attachments/assets/454eedb8-7681-453a-a0f1-6cc68ed41568" />

**Figure 1.** Pipeline overview. Major snakemake steps are summarized in green and detailed in the solid boxes. Arrows represent outputs which become inputs for the next steps in the pipeline. White dotted-line boxes indicate user-supplied parameters or files, and green dotted-line boxes indicate end-point outputs. The major steps include:
1. Download metadata from NCBI SRA, format it, and merge across studies
2. Download raw .fastq files from SRA and import into QIIME2
3. Merge and denoise reads with QIIME2, vsearch, and deblur, filter and merge results across studies
4. Phylogenetically place ASVs onto reference trees (using a stepwise approach described below)
5. Import data into R and phyloseq and generate a simple taxonomic barplot grouped by habitat type and a map showing samples which recovered sequences from TOI

## Description

This Snakemake pipeline is designed as a tool to interrogate publicly available 18S rRNA metabarcoding datasets for sequences from any given eukaryotic of interest (TOI). It takes a list of SRA accession IDs as input, as well as several user-defined parameters and input files detailed below, and outputs taxonomically-assigned sequences and count tables from your TOI, along with sample metadata and a couple of basic figures generated in R which illustrate the geographic locations and habitat types of samples containing TOI and the relative abundances of TOI across habitat types. As it imports all of the relevant files into R/phyloseq, it is meant as a starting point for your analysis. Example analyses conducted using this pipeline can be seen in this preprint (3), which conducted a meta-analysis of global marine oxygen-depleted sequencing datasets to explore the diversity and distribution of marine anaerobic ciliates. 


This pipeline is designed to process Illumina short reads which were amplified from any region of the 18S rRNA gene. This works because we use a phylogenetic placement method, which places short reads onto a given full-length reference phylogenetic tree. In this way, ASVs from different primer sets or even different regions of the 18S rRNA gene can be analyzed together. As well as being a good tool for large scale meta-analyses, this method is particularly useful for groups that are not well-represented in reference databases and can retain reads which were originally designated as "unassigned" or "unclassified" with traditional pairwise taxononomic assignment methods. To learn more check out the following references: (2,4,5).


## Setup

It can be run either locally on your computer or on a high-performance computing cluster (HPC). If possible, we strongly recommend the latter, as some of the scripts are memory-intensive. If running on an HPC, you will need to submit it as a batch job (see `submit.sh`) and fill in the cluster configuration file (`cluster.yaml`), which details the computational resources requested for each step in the pipeline. The only requirements are snakemake, conda, and PaPaRa installations on the machine where it is running from.

If you are running locally on your computer, you will need to have [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html), [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html), and [PaPaRa](https://cme.h-its.org/exelixis/web/software/papara/index.html) installed. You need to add the PaPaRa executable to your PATH:

```
mkdir -p "$HOME/bin"
cp "/full/path/to/papara_nt-2.5/papara" "$HOME/bin/papara"
chmod +x "$HOME/bin/papara"
echo 'export PATH="$HOME/bin:$PATH"' >> "$HOME/.zshrc"
source "$HOME/.zshrc"
```
Then run the Snakefile with:

```
snakemake -s Snakefile --cores all --use-conda --rerun-incomplete
```

Every rule that needs a specific software environment declares its own `conda:` environment file (see `envs/`), so `--use-conda` builds each one automatically the first time it's needed - you do not need to create these environments by hand.

### Folder set-up:

```
Snakemake/
│   ├── submit.sh
│   ├── config.yaml
│   └── cluster_config.yaml
├── reference/
│   ├── reference_database.fasta*
│   ├── reference_database.qza*
│   ├── classifier.qza*
│   ├── ENVO_IDs.csv
│   ├── ENVO.csv
│   ├── euk_tree.tree**
│   ├── eukaryotic_reference_tax.txt**
│   ├── eukaryotic_reference_tree.fasta**
│   ├── eukaryotic_reference_tree.phy**
│   └── ***
├── scripts/
│   ├── clean_metadata.py
│   ├── extract_fasta.py
│   ├── filter_by_taxopath.py
│   ├── filter_placements.py
│   ├── generate_figures.R
│   ├── make_empty_qza.py
│   └── prepare_phyloseq_objects.R
└── envs/
│   ├── parallelfastqdump.yaml
│   ├── phylo_placement.yaml
│   ├── pysradb.yaml
│   ├── qiime2-amplicon-2026.1.yaml
│   └── r_figures.yaml
├── data/
│   └──raw/
├── results/
└── visualization/

```
**data/, results/, and visualization/ folders will populate while the pipeline is running and their names are specified in the config.yaml file (raw_data, output, and visualization).**

\* For `classifier.qza`, we used a taxonomic classifier trained on PR2 [PR2](https://pr2-database.org/) v5.1.1; [instructions to train your own classifier are here](#generating-the-large-reference-files-yourself). The QIIME2 version that you are running must be the same as the one you use to train the classifier. Similarly, the `reference_database.fasta` and `reference_database.qza` are from PR2 v5.1.1 and obtained with `qiime rescript get-pr2-data`. These three files are excluded from the git repo (see `.gitignore`) because they're too large to distribute via git.

\*\* Eukaryotic reference trees and files were obtained from [(4)](https://www.zotero.org/google-docs/?ilQoQ1).

\*\*\* You must provide phylogenetic reference trees for your TOI. If you are surveying a taxonomic group within Ciliophora, we provide the relevant files: `ciliate_reference_tax.txt`, `ciliate_reference_tree.fasta`, `ciliate_reference_tree.phy`, which were obtained and prepared from [(5)](https://www.zotero.org/google-docs/?ngbr9x). See the section "How do I find or generate appropriate reference trees for phylogenetic placement?"(#How do I find or generate appropriate reference trees for phylogenetic placement?)


## Getting started: a tutorial

This section walks through running the pipeline beginning to end. 

### 1. What gets installed automatically, and what doesn't

Everything that runs *inside* a Snakemake rule - QIIME2, `epa-ng`, `gappa`, `raxml-ng`, R and all of its packages (`phyloseq`, `sf`, `rnaturalearth`, etc.) - is provisioned automatically by `snakemake --use-conda` from the `.yaml` files in `envs/`. You never install these yourself; the first time a rule needs one, Snakemake builds it from the matching env file and reuses it on every subsequent run.

The things that are *not* handled this way, and that you are responsible for having on the machine that launches the pipeline:
- **Snakemake** itself and **conda** (see [Setup](#setup) above for install links).
- **PaPaRa** must be installed separately or made available by an HPC administrator: see [Setup](#setup)
- **`reference/reference_database.fasta`, `reference/reference_database.qza`, and `reference/classifier.qza`** - these three files are excluded from the git repo (see `.gitignore`) because they're too large to distribute via git. See the next section for how to generate them yourself; it only needs to be done once.

If you're on an HPC that uses Environment Modules/Lmod, you'll also typically need to load `conda` itself via a module before the `snakemake`/`conda` commands are available - check with your cluster's documentation. `submit.sh` has commented-out `module load` lines showing what we use on our HPC as a starting point.

#### Generating the large reference files yourself

These are built from [PR2](https://pr2-database.org/) via QIIME2's `rescript` plugin, using the exact same `envs/qiime2-amplicon-2026.1.yaml` environment the rest of the pipeline uses (important - see the note on version-matching below). Build and activate that environment once:

```
conda env create -f envs/qiime2-amplicon-2026.1.yaml -n qiime2-amplicon-2026.1
conda activate qiime2-amplicon-2026.1
```

Then, from the repository root:

```
# 1. Download and format the PR2 reference sequences + taxonomy
qiime rescript get-pr2-data \
    --p-version 5.1.0 \
    --o-pr2-sequences reference/reference_database.qza \
    --o-pr2-taxonomy reference/pr2_taxonomy.qza

# 2. Export a plain-FASTA copy alongside the .qza (some tooling expects the raw file)
qiime tools export \
    --input-path reference/reference_database.qza \
    --output-path reference/_export_tmp
mv reference/_export_tmp/dna-sequences.fasta reference/reference_database.fasta
rmdir reference/_export_tmp

# 3. Train the naive Bayes taxonomic classifier used by `assign_tax`
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads reference/reference_database.qza \
    --i-reference-taxonomy reference/pr2_taxonomy.qza \
    --o-classifier reference/classifier.qza
```

`get-pr2-data` requires a stable internet connection (it downloads directly from PR2) and step 3 is the most memory/time-intensive of the three - expect it to take a while and to need a reasonable amount of RAM, so run it on a compute node/allocation if you're on an HPC.

Code for how to train a taxonomic classifier on the SILVA database can be found [here](https://github.com/tripitakit/qiime2class). Some databases provide taxonomic classifiers for QIIME2 also on their websites. 

### 2. Configuring `config.yaml`

`config.yaml` is where you describe your input datasets and every tunable parameter. Under the "Projects" heading, list the names of the datasets you will be using. Then, fill in all the details for each project:

- **`SRAid`** - BioProject accession ID
- **`primerF`** and **`primerR`**: Forward and reverse primer sequence used to amplify DNA, 5' --> 3'
- **`max_diffs_merge`** - The maximum amount of bp pairwise differences allowed in the overlap region during vsearch's merging step. A good value for this is ~20% of the overlap length. Here the default is set to 40. If most of your reads are not merging, you may want to play around with this parameter and minovlen_merge; a discussion of these parameters can be found [here](https://forum.qiime2.org/t/question-regarding-parameters-used-in-qiime-vsearch-join-pairs/12289).  
- **`minovlen_merge`** - Minimum overlap length for vsearch to merge paired-end reads. Most of the problems with vsearch are reported to be with really short overlaps (like 10bp), so setting this value to 50 is essentially the same as setting it to 200. The default in our pipeline is 50; if you expect a shorter overlap then you should lower this value. Expected overlap is calculated with (2*read length) - amplicon length. 
- **`p_trim_length`** - This parameter trims all merged reads to this length before denoising with deblur; all reads which are shorter than this are dropped. It is important to include this if your read lengths are hetereogenous as deblur requires all reads are the same length. Also, for example, if you are comparing studies which used the same primer set, it is crucial to set this parameter to the same value for all studies so that you can compare the exact same region of the 18S rRNA gene. To determine what value to set this to, look at the -merged.qzv visualization: see section [5. Checking on a run / what to expect](#5. Checking on a run / what to expect). **If you do not want to trim, then set this to -1**.

Then there is a shared settings section, which is where you will detail the locations of your output folders and input reference files. Here are the values that you will likely need to change:

- **`p_include_unassigned`, `p_include_taxa`** - In `qiime taxa filter-table`/`filter-seqs` `--p-include` values are used to filter all of the denoised sequences from each project for downstream processing. We include unassigned sequences because phylogenetic placement is often able to assign reads where pairwise methods fail. We set `p_include_unassigned` to `Unassigned,Eukaryota,Eukaryota;TSAR,Eukaryota;TSAR;Alveolata` and `p_inlcude_taxa` to `Ciliophora` because we wanted to obtain ciliate sequences; **you need to change these values based on your TOI**.
- **`search_term_taxa`, `search_term_clade1`, `search_term_clade2`, `search_term_clade3`** - the taxopath substrings used to extract placements from each placement step. `search_term_taxa` is used to extract placements from the Eukaryote tree. Sequences assigned to this group will then be placed onto the Taxon tree, and then placements from that tree are extracted with `search_term_clade1`, and optionally `search_term_clade2` and `search_term_clade3`, and then placed on their corresponding Clade trees. In our example, `search_term_taxa` = `Ciliophora`, `search_term_clade1` = `Armophorea`, `search_term_clade2` = `Plagiopylea`, `search_term_clade3` = `Anaerocycliididae`.
- **`EUK_TREE`, `MSA_FASTA_EUK`, `MSA_PHYLIP_EUK`, `CLADES_EUKS`** and the equivalent `TAXON_*`, `CLADE1_*`, `CLADE2_*`, `CLADE3_*` blocks - the tree/alignment/taxonomy files for each phylogenetic placement step (see "How do I find or generate appropriate reference trees for phylogenetic placement?" below for how to build your own set for a different TOI).

Here are the remaining values: 
- **`taxonomy_fasta`, `trained_ref_database`, `REFFASTA`** - paths to the reference FASTA/classifier described in the [Setup](#setup).
- **`epa_filter_acc_lwr`, `epa_filter_max`** - EPA-NG's placement-filtering thresholds. EPA-ng sorts the placements by their likelihood weight ratio (LWR) in descending order and  adds branches to the output file until their combined sum meets or exceeds your specified threshold (epa_filter_acc_lwr), and then in the output .jplace file retains the top X number of placements as determined by epa_filter_max. We have set those to 0.99 and 100, respectively, in order to give us the best possible idea of the placement distribution, for downstream taxonomic assignment. For more reading about how the placement steps work, refer to (4).
- **`gappa_mass_norm`** - a parameter of `gappa examine heat-tree`, which generates a figure of your reference phylogenetic tree with branches colored by the total mass of placements at that branch. This parameter controls how the total placement masses are normalized before mapping them to the tree. We set this to "absolute"; in this case, only the top placement is shown, and so the heatmap corresponds to the total # of ASVs which were assigned to each branch. For more reading on how gappa tools work, refer to (6).
- **`gappa_consensus_thresh`** - a parameter of `gappa examine assign`, which uses the phylogenetic placement results and a given taxonomic reference file to taxonomically assign query ASVs. This parameter controls the minimum proportion of descendant nodes required to agree on a taxonomic label when resolving inner nodes on the reference tree. 
- **`edpl_threshold`** - The maximum Expected Distance between Placement Locations (EDPL) allowed for a placement to be counted as high-confidence; this is what separates the raw placement output from the `filtered_*_LWR_EDPL.tsv` files that feed into the next placement (and, for the clade-level trees, into the final R figures). We set this value to 0.05.
- **`p_min_length`, `filter_minquality`, `primer_err`** - default QC parameters (minimum read length, minimum quality score, and primer-matching error tolerance) applied across all projects unless overridden per-project.

### 3. Configuring `cluster.yaml` (HPC only)

If you're running on an HPC via `submit.sh`, `cluster.yaml` tells Snakemake what computational resources to request from your scheduler (we use Slurm) for each rule. There are inputs for `partition`, `time`, `mem`, `ntasks`, and `nodes` for the most resource-intensive rules and a `__default__` block used for any rule without its own entry. **You need to edit this file according to your HPC.**

### 4. Running it

**On an HPC:** review/edit `submit.sh` (partition names, time limits, and any `module load` lines specific to your cluster), then submit it as a batch job:
```
sbatch submit.sh
```

**Locally:** as shown in [Setup](#setup):
```
snakemake -s Snakefile --cores all --use-conda --rerun-incomplete
```

Either way, the first run will take a while, since every conda environment needs to be built from scratch and (for the SRA-derived projects) raw reads need to be downloaded. Subsequent runs reuse the built environments and any already-completed outputs.

### 5. Outputs

- **`results/`** fills in results per-project: `results/<Project name>/...`,
    - for all of the projects merged together: `results/_merged/`,
    - the phylogenetic placement results: `results/_placement/` (the Eukaryote tree),
    - `results/_Taxon_placement/` (the Taxon tree),
    - and `results/_clade1_placement/`, `results/_clade2_placement/`, `results/_clade3_placement/` (the clade-specific trees)
- **`results/visualizations/`** stores `.qzv` files you can upload to [view.qiime2.org](https://view.qiime2.org) to visualize to see what the reads look like after each step and diagnose any issues that might have occurred.
    - Check `-PE-demux-noprimer-merged-filtered.qzv` for read length distributions prior to deblur to choose `p_trim_length`
    - Check `-deblur-stats.qzv` to see the statistics for the deblur run and determine if it ran sucessfully. See [here](https://forum.qiime2.org/t/deblur-stats-qzv-file-meaning-and-interpretation/3748/7) for a discussion on what these stats mean.
- **`results/_figures/`** holds the final R outputs: the phyloseq object (`ps.rds`, which can be imported into R for further analyses), the taxa barplot and sample map (as both PDF and PNG), the cleaned/combined metadata and taxonomy tables, and **`metadata_gaps.txt`** - a report of exactly which metadata columns are incomplete for which projects, generated automatically every run (see [Known limitations](#known-limitations) for why we report gaps rather than trying to auto-fill them).
- **On an HPC**, per-rule logs land in `slurm-logs/<rule>.<jobid>.log`; the top-level `slurm-<jobid>.out` file is the master Snakemake driver log showing overall progress through the DAG.
- If a run stops partway through, both `sbatch submit.sh` and the local `snakemake` command above are safe to re-run - Snakemake picks up from whatever's already been produced rather than starting over.
    - If you cancel the run yourself, you will need to run `snakemake --unlock` from the directory where the snakefile is located prior to rerunning. 


### What kind of sequencing data can I include as input?

The first step is to find SRA BioProjects containing paired-end 18S rDNA reads which were sequenced with Illumina and which contain samples that you're interested in. They could be from a habitat type you're interested in or suspect contains TOI, for example. Since authors are mostly required to deposit their raw sequencing data into SRA or some other repository and report accession IDs, a literature review is a good place to start to collect SRA accession IDs and associated project information (we recommend you record the forward and reverse primer sequences used to amplify, the expected amplicon length, and the read length; and optionally, any other information you want to have associated with the project, like habitat type, location, whether they amplified DNA or cDNA, whether they used one primer pair or a nested primer strategy, etc.). The website [sra-explorer.info](http://sraexplorer.com) is also helpful to search NCBI SRA for samples from given habitat types. If you know of another tool that makes searching SRA easier, email us and let us know and we can add it here!

> If you want to include studies that are not on SRA, this is a feature that we will support in a future update. 

### How do I find or generate appropriate reference trees for phylogenetic placement? 

This Snakemake pipeline uses a multi-level phylogenetic placement scheme, as described in (5) and used in e.g. (3,6), and summarized in Figure 2: 

<img width="794" height="1123" alt="multilevel_placement" src="https://github.com/user-attachments/assets/89fd3415-d59a-4a88-9656-e54812fc86d3" />
**Figure 2** 3-tier multilevel placement scheme, adapted from (5). 1: Unassigned sequences (A) are placed onto the Eukaryote tree, and sequences from your Taxon of choice (purple) are extracted. 2: Sequences which are assigned to you Taxon of choice, either via placement onto the euk tree (A) or taxonomic assignment with QIIME2 (B, C) are placed onto the Taxon tree. Branches that are associated with a clade tree are colored accordingly (orange and green). 3: Clade trees. The backbone tree and clade trees overlap each other such that each clade tree is represented by branches in the backbone tree. Three sequences A, B, and C are placed 


## References:

1. Hu SK. shu251/tagseq-qiime2-snakemake. 2026. Available from: github.com/shu251/tagseq-qiime2-snakemake
2. Ewers I, Rajter L, Czech L, Mahé F, Stamatakis A, Dunthorn M. Interpreting phylogenetic placements for taxonomic assignment of environmental DNA. J Eukaryot Microbiol. 2023;70(5):e12990.
3. Schrecengost A, Frates E, Al-Haj A, Fulweiler RW, Beinart R. A meta-analysis of environmental sequencing data reveals the global distribution and hidden diversity of marine anaerobic ciliates. bioRxiv. 2025;2025–12.
4. Czech L, Stamatakis A, Dunthorn M, Barbera P. Metagenomic Analysis Using Phylogenetic Placement—A Review of the First Decade. Front Bioinforma. 2022
5. Czech L, Barbera P, Stamatakis A. Methods for automatic reference trees and multilevel phylogenetic placement. Bioinformatics. 2019 Apr 1;35(7):1151–8.
6. Mahé F, de Vargas C, Bass D, Czech L, Stamatakis A, Lara E, et al. Parasites dominate hyperdiverse soil protist communities in Neotropical rainforests. Nat Ecol Evol. 2017 Apr;1(4):0091.
7. Rajter Ľ, Dunthorn M. Ciliate SSU-rDNA reference alignments and trees for phylogenetic placements of metabarcoding data. Metabarcoding Metagenomics. 2021 Aug 30;5:e69602.
8. Czech L, Barbera P, Mahe F. lzech/gappa. 2025. Available from: github.com/lczech/gappa
