# Metacontam v0.0.2
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/metacontam/README.html)

**Metacontam** is a contaminant detection tool for shotgun metagenomic sequencing data.  
Metacontam identifies contaminant species using a blacklist-guided Louvain algorithm and conANI comparison. It is best suited for low-biomass metagenomic samples.

<img width="434" height="144" alt="Metacontam" src="Metacontam_image.png" />

---

## Overview

Metacontam detects microbial contamination through a 12-step pipeline:

```
Kraken2 → Bracken → Count Matrix → Prevalence Threshold
→ Network Analysis (NetCoMi + Louvain) → Genome Retrieval
→ Bowtie2 Alignment → MASH Distance → inStrain → Final Prediction
```

---

## The blacklist

Metacontam ships with a **blacklist** of 1,782 species-level NCBI taxids
(`Metacontam/species_under_blacklist.pkl`), compiled from published surveys of
reagent and laboratory contaminants. It is a fixed resource built independently
of any dataset analysed here, so it is not tuned to the benchmark.

**The blacklist is a prior, not a filter.** Metacontam never removes a species
just because it is on the list, and never keeps one just because it is absent.
A blacklisted species that behaves like a genuine community member is retained,
and a species absent from the list is called a contaminant whenever the network
and ANI evidence say so. That is the point of the method: the blacklist only
tells the algorithm *where to start looking*.

It enters the pipeline in three places:

| Stage | How the blacklist is used |
|---|---|
| **Prevalence threshold** | The adaptive cutoff starts from `median_black`, the median prevalence of the (up to) 100 most prevalent blacklisted taxa present in the data, then shifts by `0.7 × median_corr × 10 / sqrt(n)`. Contaminants are prevalent across samples, so the blacklist supplies a data-driven anchor instead of an arbitrary constant. |
| **Community seeding** | Blacklisted taxa that survive the filter become **seeds** for the seeded Louvain step. Among the detected communities, Metacontam selects the one containing the seeds; every other member of that community is a contaminant *candidate*, whether or not it is blacklisted. This is how species missing from the list are recovered. |
| **ANI cutoff** | The conANI percentile is adjusted as `effective_percentile = 0.6 − 0.6 × bl_ratio`, where `bl_ratio` is the fraction of ANI-evaluable candidates that are blacklisted. A candidate community dominated by known contaminants relaxes the cutoff; one with few known contaminants keeps it strict. |

If no blacklisted taxon survives the prevalence filter, the adaptive threshold
cannot be defined and the run stops with an explicit message rather than looping.
This usually means the classifier database uses a different taxonomy; use
`--filtered-matrix` to supply a taxon set directly and skip the adaptive filter.

---

## Kraken2 database requirement

Metacontam requires a Kraken2 database that includes raw `library.fna` files (needed for genome retrieval) and has `bracken-build` applied before use.

### Option A — Download pre-built database

Building the standard Kraken2 database can take many hours and is prone to download errors. We provide a pre-built version (August 2023, includes bacteria, archaea, viral, human, plasmid, and UniVec_Core; Bracken indices for 50/100/125/150 bp read lengths included):

```bash
wget https://pub-3e5157a62ab94e1eb91f43f8fdc5c5f4.r2.dev/kraken2_standard_db_Aug2023.tar
tar -xf kraken2_standard_db_Aug2023.tar
```

> **Note**: The database is ~284 GB. Make sure you have sufficient disk space before downloading.

---

### Option B — Build the latest database

```bash
# Build standard Kraken2 DB (bacteria + archaea + viral + human + plasmid + UniVec_Core)
kraken2-build --standard --db /path/to/kraken2_db --threads 16

# Run bracken-build — set -l to your actual read length (e.g. 100, 150, 250)
bracken-build -d /path/to/kraken2_db -t 16 -k 35 -l <READ_LENGTH>
```

> **Note**: If download fails due to rsync issues, add `--use-ftp`:
> ```bash
> kraken2-build --standard --db /path/to/kraken2_db --threads 16 --use-ftp
> ```

---

## Installation

### Option 1 — Install with conda

**Step 1. Install Metacontam**
```bash
conda create -n metacontam python=3.10
conda activate metacontam
conda install -c bioconda -c conda-forge metacontam
pip install instrain
```

> **Note**: inStrain is installed separately via pip to avoid conda dependency conflicts.

**Step 2. Verify installation**
```bash
metacontam --help
```

<br>

---


### Option 2 — Install with environment file (mamba)

> All external tools and Python dependencies are installed automatically via a single environment file.

**Step 1. Install Miniforge (skip if mamba is already installed)**
```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
# Restart your terminal after installation
```

**Step 2. Clone and install Metacontam**
```bash
git clone https://github.com/Lifemining-lab/Metacontam.git
cd Metacontam
mamba env create -f environment.yml
conda activate metacontam
pip install .
```

**Step 3. Verify installation**
```bash
metacontam --help
```

<br>

---

### Option 3 — Install from source (manual)

> Most flexible option — install each external tool yourself. Useful if you already have some tools installed or need fine-grained control, but requires more manual steps.

**Step 1. Install external tools**

Ensure the following tools are installed and available in your `$PATH`:

| Tool | Tested version |
|------|----------------|
| [Kraken2](https://github.com/DerrickWood/kraken2) | 2.1.3 |
| [Bracken](https://github.com/jenniferlu717/Bracken) | 2.8 |
| [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2) | 2.5 |
| [Samtools](http://www.htslib.org/) | 1.10 |
| [MASH](https://github.com/marbl/Mash) | 2.1 |
| [inStrain](https://github.com/MrOlm/inStrain) | 1.9.1 |
| [NetCoMi](https://github.com/stefpeschel/NetCoMi) (R package) | 1.1.0 |

**Step 2. Clone and install Metacontam**
```bash
git clone https://github.com/Lifemining-lab/Metacontam.git
cd Metacontam
mamba env create -f environment_min.yml
conda activate metacontam
pip install instrain
pip install .
```

**Step 3. Verify installation**
```bash
metacontam --help
```

---

## Quick Start

```bash
metacontam \
  --metadata    metadata.tsv \
  --output      ./output \
  --DB          /path/to/kraken2_db \
  --read-length <READ_LENGTH>   # set to your actual read length (e.g. 150)
```

### Metadata format (`metadata.tsv`)

This is the only input file you write yourself. It lists the samples to analyse
and where their reads live. **Tab-separated, four columns, no header row:**

| # | Column | Description |
|---|--------|-------------|
| 1 | Sample name | Identifier used for every output file of that sample. Must be unique. Avoid whitespace; `_` is used internally to split the sample ID from the rest of a filename, so a name without `_` is safest. |
| 2 | Sample type | Group label (e.g. body site, batch). **Reserved for future batch-aware decontamination and unused by the current pipeline** — any placeholder such as `sample` is accepted. |
| 3 | R1 path | Forward reads, FASTQ (`.fastq` or `.fastq.gz`). Absolute paths are safest. |
| 4 | R2 path | Reverse reads. Paired-end data is required. |

```
SRR1000001	skin	/data/SRR1000001_R1.fastq.gz	/data/SRR1000001_R2.fastq.gz
SRR1000002	skin	/data/SRR1000002_R1.fastq.gz	/data/SRR1000002_R2.fastq.gz
SRR1000003	nasal	/data/SRR1000003_R1.fastq.gz	/data/SRR1000003_R2.fastq.gz
```

> Separate the columns with **tab characters, not spaces**. A spreadsheet export
> saved as "Tab delimited text" works; a CSV does not.

**How many samples do I need?** Metacontam infers contamination from
co-occurrence across samples, so it needs enough samples for the correlation
network to be meaningful. The datasets in the paper range from 20 to 344
samples. Below roughly 10 samples the network stage becomes unreliable; the
BH-FDR on the edges (`--fdr-alpha`, on by default) is there to keep the edge set
honest at small *n*.

> **Note**: No negative controls are required. Metacontam is designed for
> datasets where blanks were never sequenced.

---

## Arguments

### Required
| Argument | Description |
|----------|-------------|
| `--metadata` | Metadata TSV (sample name, sample type, R1 path, R2 path) |
| `--output` | Output directory |
| `--DB` | Kraken2 database directory |
| `--read-length` | Read length for Bracken re-estimation |

### Optional
| Argument | Default | Description |
|----------|---------|-------------|
| `--numcore` | 1 | Number of threads |
| `--min-reads` | 2 | Minimum reads to call a species |
| `--min-cor` | 0.45 | Pearson correlation threshold for network edges |
| `--kraken-report` |  | Pre-computed Kraken2 report directory (skips Kraken2) |
| `--bracken-report` |  | Pre-computed Bracken report directory (skips Bracken) |
| `--dist-matrix` |  | Pre-computed MASH distance matrix (skips MASH) |
| `--candidate-genome` |  | Pre-fetched candidate FASTA (skips genome retrieval) |
| `--bam-dir` |  | Pre-made BAM directory (skips alignment) |
| `--filtered-matrix` |  | Reuse an existing `kraken_filtered_matrix.txt` (skips the adaptive prevalence filter) |

### Network (new in 0.0.2)
| Argument | Default | Description |
|----------|---------|-------------|
| `--fdr-alpha` | 0.1 | Benjamini-Hochberg FDR level applied to the edges that pass `--min-cor`. `0` disables it, reproducing v0.0.1. A fixed correlation threshold is statistically permissive at small sample sizes; the FDR corrects for that. It is **not** applied inside the adaptive prevalence filter, so the taxon set stays independent of the edge rule. |
| `--zero-method` | `multRepl` | Zero handling before the CLR. `multRepl` replaces zeros using a detection limit, `bayesMult` is Bayesian-multiplicative and uses none, `pseudo` adds a pseudocount at the count level. |
| `--zero-dl` | NetCoMi's `1e-3` | Detection limit for `multRepl`, as a proportion of the composition. Zeros are filled with `0.65 * dl`. |
| `--pseudocount` | 0.5 | Pseudocount for `--zero-method pseudo`. |

> On all ten real datasets in the paper (20-344 samples) `--fdr-alpha 0.1` removes **zero** edges, so the published results are unchanged. It only bites at small *n*.

---

## Output

| File / Directory | Description |
|-----------------|-------------|
| `Kraken_dir/` | Kraken2 classification reports |
| `Bracken_dir/` | Bracken species-level reports |
| `Abundance.txt` | Mean abundance per species across samples |
| `Prevalence.txt` | Prevalence per species per sample type |
| `kraken_filtered_matrix.txt` | Filtered count matrix used for network analysis |
| `Network_Output/` | Edge list, Louvain partition, community figure |
| `Genome_dir/Candidate.fasta` | Candidate contaminant genome sequences |
| `Bamfiles/` | Sorted BAM files |
| `mash_sketches/` | MASH sketch files |
| `dist_matrix.txt` | Pairwise MASH distance matrix |
| `pair_output.tsv` | Stratified sample pairs for inStrain |
| `ISfiles/` | inStrain profile outputs |
| `IScompare/` | inStrain compare outputs |
| `merged_IS_compare_Table.tsv` | Merged pairwise ANI comparison table |
| `Final_prediction.txt` | Final contaminant/non-contaminant classification |

---

## Changelog

### 0.0.2

- **`--fdr-alpha` (default 0.1)** — Benjamini-Hochberg FDR on the correlation
  edges, on top of `--min-cor`. Pass `0` for v0.0.1 behaviour.
- **`--zero-method` / `--zero-dl` / `--pseudocount`** — the zero-replacement
  step before the CLR is now selectable. The default is unchanged from v0.0.1.
- **`--filtered-matrix`** — reuse an existing filtered matrix to hold the taxon
  set fixed across runs.
- **No more infinite loop when the blacklist is empty.** The adaptive prevalence
  search used to raise `min_reads` and retry forever if no blacklist taxon
  survived the filter. Raising `min_reads` can only remove taxa, so it now
  aborts with an explicit message.
- **Bracken failures are no longer silently treated as "no species detected".**
  Each sample is recorded in `Bracken_dir/bracken_status.tsv` as `ok`, `empty`
  (genuine non-detection), `failed_no_report`, or `failed_bracken`, with a
  summary printed at the end.
- **Explicit "not evaluable" handling.** A run with no valid ANI comparison, or
  with no comparison file at all, used to die with an `IndexError` or
  `ValueError`. Both are normal outcomes for low-biomass data, so the run now
  says what is missing and writes an empty `Final_prediction.txt`.
- **The `min_reads` retry cap was removed.**

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Citation


If you use Metacontam, please cite our preprint:

> Jo J, Lee H, Baek JW, Lee S, Singh V, Shoaie S, Mardinoglu A, Choi J, Lee S. 
> **Metacontam: A Negative Control-Free Decontamination Method for Metagenomic Analysis.** 
> *bioRxiv* 2026.04.26.720876; doi: https://doi.org/10.64898/2026.04.26.720876

> **Note**: This work is currently a preprint and has not yet been peer-reviewed. 
> Citation will be updated once the peer-reviewed version is published.
