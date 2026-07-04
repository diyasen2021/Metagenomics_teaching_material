
# From Raw Reads to Microbial Genomes

**Based on:** [GTN Tutorial — Building and Annotating MAGs from Short Metagenomics Paired Reads](https://training.galaxyproject.org/training-material/topics/microbiome/tutorials/mags-building/tutorial.html)  
**Environment:** Linux / HPC / WSL2 (command-line)  
**Binner shown:** MetaBAT2 (one of four in the full workflow; sufficient to learn the concept)  
**Estimated time:** ~6–10 hours for a first run  
**Level:** Intermediate

---

## What is a MAG and why do we build them?

A **Metagenome-Assembled Genome (MAG)** is a near-complete or complete microbial genome reconstructed computationally from a mixed environmental sample — without ever culturing the organism. Because most environmental microbes resist cultivation, MAGs are often our only window into their genomes.

**The pipeline goal:** raw shotgun reads → trimmed reads → host-depleted reads → assembled contigs → binned genomes → refined, de-replicated MAGs → annotated microbial genomes.

---

## Dataset

This tutorial uses two honey bee gut metagenome samples from Sbaghdi et al. 2024:

| Sample | Treatment | SRA Accession |
|--------|-----------|---------------|
| IP2 | Infected + probiotic | SRR24759598 |
| I1 | Infected only | SRR24759616 |

The honey bee gut is dominated by ~5 bacterial lineages, making it manageable as a teaching dataset while still representative of real complexity.

---

## Pipeline Overview

```
Raw reads (FASTQ)
    │
    ├─── Step 1: Quality control & trimming     (fastp + MultiQC)
    ├─── Step 2: Host & contaminant removal     (Bowtie2)
    ├─── Step 3: Assembly                       (MEGAHIT)
    ├─── Step 4: Assembly QC                   (QUAST)
    ├─── Step 5: Binning                       (MetaBAT2)
    ├─── Step 6: Bin quality check             (CheckM2)
    ├─── Step 7: Bin refinement                (Binette)
    ├─── Step 8: De-replication                (dRep)
    ├─── Step 9: Final QC                      (CheckM2 + QUAST)
    ├─── Step 10: Abundance estimation         (CoverM)
    ├─── Step 11: Taxonomic assignment         (GTDB-Tk)
    └─── Step 12: Functional annotation        (Bakta)
```

---

## Environment Setup

All tools are installed via conda. If you don't have conda:

```bash
# Install Miniconda (Linux / WSL2)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict
```

### Create environments

We use separate environments to avoid dependency conflicts — a real-world best practice.

```bash
# QC and preprocessing
conda create -n qc -y fastp multiqc bowtie2 samtools sra-tools
conda activate qc

# Assembly
conda create -n assembly -y megahit quast
conda activate assembly

# Binning
conda create -n binning -y metabat2
conda activate binning

# Bin quality and refinement
conda create -n checkm2 -y checkm2
conda activate checkm2

conda create -n binette -y binette
conda activate binette

# De-replication
conda create -n drep -y drep
conda activate drep

# Abundance
conda create -n coverm -y coverm
conda activate coverm

# Taxonomy
conda create -n gtdbtk -y gtdbtk
conda activate gtdbtk

# Functional annotation
conda create -n bakta -y bakta
conda activate bakta
```

### Set up project structure

```bash
mkdir -p mag_tutorial/{00_raw,01_qc,02_hostfree,03_assembly,04_quast,\
05_binning,06_checkm2_raw,07_refinement,08_drep,09_final_qc,\
10_abundance,11_taxonomy,12_annotation}
cd mag_tutorial
```

---

## Step 0: Download the data

```bash
conda activate qc

# Download both SRA runs (paired-end FASTQ)
cd 00_raw
fasterq-dump --split-files --threads 8 SRR24759598
fasterq-dump --split-files --threads 8 SRR24759616
gzip *.fastq        # compress to save space
cd ..
```

> **Tip — on a slow connection:** download directly from ENA instead:
> ```bash
> wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR247/098/SRR24759598/SRR24759598_1.fastq.gz
> wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR247/098/SRR24759598/SRR24759598_2.fastq.gz
> wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR247/016/SRR24759616/SRR24759616_1.fastq.gz
> wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR247/016/SRR24759616/SRR24759616_2.fastq.gz
> mv *.fastq.gz 00_raw/
> ```

> **Working on a laptop?** Subsample to 500K reads first — full metagenomes need HPC memory:
> ```bash
> conda install -n qc -y seqtk
> conda activate qc
> seqtk sample -s 42 00_raw/SRR24759598_1.fastq.gz 500000 > 00_raw/SRR24759598_sub_1.fastq
> seqtk sample -s 42 00_raw/SRR24759598_2.fastq.gz 500000 > 00_raw/SRR24759598_sub_2.fastq
> # repeat for SRR24759616
> ```

---

## Step 1: Quality Control and Trimming

### Why this step matters

Raw Illumina reads contain:
- **Adapter sequences** — synthetic oligonucleotides from library prep, not biology. If left in, the assembler treats them as real sequence and creates chimeric contigs.
- **Low-quality bases** — mostly at the 3' end of reads. In a metagenome with low-coverage organisms, errors look like real variants.
- **Short reads** — reads trimmed below ~50 bp become nearly useless for assembly (can't span repeats, map ambiguously).

### Install (if not already done)

```bash
conda activate qc
# fastp and multiqc already installed in qc environment above
```

### Run fastp on each sample

```bash
for SAMPLE in SRR24759598 SRR24759616; do
    fastp \
        --in1  00_raw/${SAMPLE}_1.fastq.gz \
        --in2  00_raw/${SAMPLE}_2.fastq.gz \
        --out1 01_qc/${SAMPLE}_1.fastq.gz \
        --out2 01_qc/${SAMPLE}_2.fastq.gz \
        --detect_adapter_for_pe \          # auto-detect adapters for paired-end
        --qualified_quality_phred 15 \     # discard bases below Q15
        --length_required 15 \             # discard reads shorter than 15 bp after trimming
        --thread 8 \
        --html 01_qc/${SAMPLE}_fastp.html \
        --json 01_qc/${SAMPLE}_fastp.json
done
```

### Aggregate reports with MultiQC

```bash
multiqc 01_qc/ -o 01_qc/multiqc_report
```

Open `01_qc/multiqc_report/multiqc_report.html` in a browser.

### What to look for

| Metric | Acceptable range |
|---|---|
| % reads passing filter | > 80% |
| Per-base quality (after trimming) | Q ≥ 30 for most of the read |
| Adapter content | Near 0% after trimming |
| Mean read length | Close to 150 bp (not dramatically shorter) |

> **Expected for this dataset:** ~43.9M reads (SRR24759598), ~42.8M reads (SRR24759616). Both already have high average Q-scores (~34), so minimal trimming occurs.

---

## Step 2: Host and Contaminant Read Removal

### Why this step matters

Your sample contains honey bee gut DNA — meaning there is bee (host) DNA present. If you assemble it, you assemble the bee genome alongside the microbes, wasting compute and polluting your gene catalog. Human contamination from lab handling is also possible. We remove both sequentially.

### Install reference genomes

```bash
conda activate qc
mkdir -p refs

# Honey bee genome (apiMel3)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/195/\
GCF_000002195.4_Amel_4.5/GCF_000002195.4_Amel_4.5_genomic.fna.gz \
-O refs/bee_genome.fna.gz

# Human genome (hg38) — large file, ~3 GB compressed
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/\
GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz \
-O refs/human_genome.fna.gz

# Build Bowtie2 indices
bowtie2-build refs/bee_genome.fna.gz   refs/bee_index
bowtie2-build refs/human_genome.fna.gz refs/human_index
```

> **Note:** Index building is slow (20–60 min for human). Pre-built indices are sometimes available from your HPC admin — ask before building.

### Remove bee reads

```bash
for SAMPLE in SRR24759598 SRR24759616; do
    bowtie2 \
        -x refs/bee_index \
        -1 01_qc/${SAMPLE}_1.fastq.gz \
        -2 01_qc/${SAMPLE}_2.fastq.gz \
        --un-conc-gz 02_hostfree/${SAMPLE}_nobee_%.fastq.gz \
        --threads 8 \
        -S /dev/null   # we discard the SAM output — we only want the unmapped reads
done
```

`--un-conc-gz` writes the unmapped read pairs to the specified files. These are reads that did **not** map to the bee genome — the ones we want.

### Remove human reads

```bash
for SAMPLE in SRR24759598 SRR24759616; do
    bowtie2 \
        -x refs/human_index \
        -1 02_hostfree/${SAMPLE}_nobee_1.fastq.gz \
        -2 02_hostfree/${SAMPLE}_nobee_2.fastq.gz \
        --un-conc-gz 02_hostfree/${SAMPLE}_clean_%.fastq.gz \
        --threads 8 \
        -S /dev/null
done
```

Your final decontaminated reads are `02_hostfree/{SAMPLE}_clean_1.fastq.gz` and `_clean_2.fastq.gz`.

> **Expected:** ~2.8% of reads map to bee (SRR24759598), ~1.1% (SRR24759616). ~0.1% map to human for both. These are clean samples.

---

## Step 3: Assembly

### What is assembly?

Assembly stitches short reads back into longer contiguous sequences (**contigs**) by finding overlapping fragments. Metagenomics makes this harder because: reads come from hundreds of organisms simultaneously; rare organisms have low coverage (few overlapping reads); and closely related strains confuse assemblers with ambiguous overlaps.

### Tool: MEGAHIT

MEGAHIT is the recommended starting point. It uses low memory relative to its alternatives, is fast, and produces reliable contigs for most community types.

| | MEGAHIT | metaSPAdes |
|---|---|---|
| RAM required | Moderate (~32 GB for large samples) | Very high (>100 GB) |
| Speed | Fast | Slow |
| Contig quality | Good | Often better for complex samples |
| Use when | Default choice, laptop/HPC | HPC, you need longer contigs |

### Run MEGAHIT

```bash
conda activate assembly

for SAMPLE in SRR24759598 SRR24759616; do
    megahit \
        -1 02_hostfree/${SAMPLE}_clean_1.fastq.gz \
        -2 02_hostfree/${SAMPLE}_clean_2.fastq.gz \
        --min-contig-len 200 \       # discard contigs shorter than 200 bp
        --out-dir 03_assembly/${SAMPLE} \
        -t 16
    
    # Copy final contigs to a named file for clarity
    cp 03_assembly/${SAMPLE}/final.contigs.fa \
       03_assembly/${SAMPLE}_contigs.fa
done
```

> **Expected:** SRR24759598 → ~534K contigs; SRR24759616 → ~572K contigs. This is normal for metagenomes — most contigs will be short and fragmented from low-abundance organisms.

### Step 4: Assess assembly quality with QUAST

```bash
for SAMPLE in SRR24759598 SRR24759616; do
    metaquast.py \
        03_assembly/${SAMPLE}_contigs.fa \
        --output-dir 04_quast/${SAMPLE} \
        --threads 8 \
        --min-contig 500   # QUAST reports stats for contigs >= 500 bp by default
done
```

Open the HTML report at `04_quast/{SAMPLE}/report.html`.

### Key metrics to understand

| Metric | What it means | Expected for metagenomes |
|---|---|---|
| # contigs (≥500 bp) | Contigs long enough for binning | ~85–90K per sample here |
| Total length | Sum of all contig lengths | ~120–140 Mb |
| Largest contig | Longest single assembled sequence | ~200–430 Kb |
| **N50** | Length where cumulative assembly = 50% total. **Higher = more contiguous.** | ~1,800–2,200 bp (normal) |
| **L50** | Minimum number of contigs to cover 50% of assembly. **Lower = better.** | ~11–13K contigs |
| % reads mapped | Fraction of input reads that ended up assembled | Should be >85% |

> **Expected:** N50 ~ 2,213 bp (SRR24759598) and 1,751 bp (SRR24759616). These low N50 values are normal for metagenomes and do not indicate a failed assembly. They reflect the presence of many low-abundance organisms with insufficient coverage for long contig extension.

---

## Step 5: Binning with MetaBAT2

### What is binning?

After assembly, your FASTA file contains hundreds of thousands of contigs from many different organisms all mixed together. **Binning** groups contigs that likely came from the same organism using two signals:

1. **Tetranucleotide frequency (TNF)** — each organism has a characteristic DNA composition "fingerprint"
2. **Coverage depth** — contigs from the same organism should have similar read depth across samples

### Prepare coverage information

MetaBAT2 needs to know how deeply each contig is covered by reads. We generate this by mapping the cleaned reads back to the assembly:

```bash
conda activate binning

for SAMPLE in SRR24759598 SRR24759616; do
    # Index the assembly
    bowtie2-build 03_assembly/${SAMPLE}_contigs.fa \
                  03_assembly/${SAMPLE}_index

    # Map reads back to contigs
    bowtie2 \
        -x 03_assembly/${SAMPLE}_index \
        -1 02_hostfree/${SAMPLE}_clean_1.fastq.gz \
        -2 02_hostfree/${SAMPLE}_clean_2.fastq.gz \
        --threads 8 | \
    samtools sort -o 05_binning/${SAMPLE}.sorted.bam
    
    samtools index 05_binning/${SAMPLE}.sorted.bam
done
```

### Generate depth profile

```bash
# For each sample, calculate per-contig coverage depth
for SAMPLE in SRR24759598 SRR24759616; do
    jgi_summarize_bam_contig_depths \
        --outputDepth 05_binning/${SAMPLE}_depth.txt \
        05_binning/${SAMPLE}.sorted.bam
done
```

`jgi_summarize_bam_contig_depths` is bundled with MetaBAT2.

### Run MetaBAT2

```bash
for SAMPLE in SRR24759598 SRR24759616; do
    mkdir -p 05_binning/${SAMPLE}_bins
    metabat2 \
        -i 03_assembly/${SAMPLE}_contigs.fa \
        -a 05_binning/${SAMPLE}_depth.txt \
        -o 05_binning/${SAMPLE}_bins/bin \
        --minContig 1500 \   # minimum contig length to include in binning
        --threads 8 \
        -v                   # verbose output
done
```

Each output file is one bin: `05_binning/SRR24759598_bins/bin.1.fa`, `bin.2.fa`, etc.

> **Expected:** MetaBAT2 produces ~34 bins for SRR24759598 and ~25 for SRR24759616. Compare this to CONCOCT (80 and 59) or SemiBin (16 and 20) — MetaBAT2 is moderate and generally reliable. The full workflow uses all four binners and refines; here we continue with MetaBAT2 bins.

---

## Step 6: Assess Raw Bin Quality with CheckM2

### What is CheckM2?

CheckM2 estimates completeness and contamination of each bin by checking for the presence of universal single-copy marker genes. These genes appear exactly once in virtually all bacteria/archaea.

- **Completeness (%)** — what fraction of expected marker genes are present in this bin
- **Contamination (%)** — what fraction of marker genes appear more than once (signals from multiple organisms mixed in the bin)

### Install CheckM2 database

```bash
conda activate checkm2
checkm2 database --download --path ~/.checkm2/
```

### Run CheckM2 on raw bins

```bash
for SAMPLE in SRR24759598 SRR24759616; do
    checkm2 predict \
        --input 05_binning/${SAMPLE}_bins/ \
        --output-directory 06_checkm2_raw/${SAMPLE} \
        --extension .fa \
        --threads 8
done
```

Look at `06_checkm2_raw/{SAMPLE}/quality_report.tsv`.

### Interpret the output

| Column | Meaning |
|---|---|
| Name | Bin filename |
| Completeness | % of marker genes present (higher = better) |
| Contamination | % of duplicated marker genes (lower = better) |
| Completeness_Model_Used | General or specific model |

> **What to expect:** MetaBAT2 bins typically have lower contamination than CONCOCT bins. Some bins will be excellent (>90% completeness, <5% contamination). Others will be poor. That is why we refine next.

---

## Step 7: Bin Refinement with Binette

### Why refine?

Raw bins from a single tool often:
- Contain contigs from multiple organisms (**contamination**)
- Have the same organism's genome split across multiple bins (**fragmentation**)

Binette (successor to metaWRAP) addresses both. It takes bins from one or more binners, scores them using CheckM2, and produces a cleaner, more complete set.

**Scoring formula:** `score = completeness − (contamination × 2)`  
Contamination is penalised twice as heavily as incompleteness.

### Install Binette database

```bash
conda activate binette
# Binette uses CheckM2 under the hood — point it to the database
export CHECKM2DB=~/.checkm2/uniref100.KO.1.dmnd
```

### Run Binette

```bash
# Pool bin directories from both samples
# (in the full workflow, bins from all samples are pooled together)
for SAMPLE in SRR24759598 SRR24759616; do
    binette \
        --bin_dirs 05_binning/${SAMPLE}_bins/ \
        --contigs  03_assembly/${SAMPLE}_contigs.fa \
        --outdir   07_refinement/${SAMPLE} \
        --min_completeness 75 \    # discard bins below 75% completeness
        --contamination_weight 2 \ # penalise contamination twice as much as incompleteness
        --threads 8
done
```

Refined bins are written to `07_refinement/{SAMPLE}/final_bins/`.

> **Expected:** Bin count drops substantially — from ~34 MetaBAT2 raw bins to ~9 refined bins for SRR24759598. The discarded bins had completeness < 75% or were absorbed into better bins. Refined bins have contamination generally < 10%.

---

## Step 8: De-replication with dRep

### Why de-replicate?

When you assemble each sample individually, the same microbial species may produce a MAG in sample 1 and again in sample 2. Keeping both:
- Splits reads between two near-identical MAGs → artificially low abundance estimates for each
- Inflates the apparent diversity of your dataset

De-replication keeps only the single best representative of each species across all samples.

### Key concept: Average Nucleotide Identity (ANI)

ANI measures genomic similarity between two microbial genomes (percentage of identical nucleotides across all aligned regions):

| ANI | Interpretation |
|---|---|
| ~100% | Nearly identical strains / clones |
| ≥95–96% | Same species (the standard threshold) |
| 80–95% | Same genus, different species |
| <80% | Different genus or higher rank |

### Pool all refined bins

```bash
mkdir -p 08_drep/all_bins

# Copy all refined bins from both samples into one directory
for SAMPLE in SRR24759598 SRR24759616; do
    cp 07_refinement/${SAMPLE}/final_bins/*.fa \
       08_drep/all_bins/${SAMPLE}_$(basename {}
done

# Simpler: just copy them all with a glob
cp 07_refinement/SRR24759598/final_bins/*.fa 08_drep/all_bins/
cp 07_refinement/SRR24759616/final_bins/*.fa 08_drep/all_bins/
```

### Run CheckM2 on the pooled bins (required by dRep)

```bash
conda activate checkm2
checkm2 predict \
    --input 08_drep/all_bins/ \
    --output-directory 08_drep/checkm2_pooled \
    --extension .fa \
    --threads 8
```

### Run dRep

```bash
conda activate drep

dRep dereplicate 08_drep/drep_output \
    --genomes 08_drep/all_bins/*.fa \
    --checkM_method checkm2 \
    --genomeInfo 08_drep/checkm2_pooled/quality_report.tsv \
    --completeness 75 \       # minimum completeness to include in dereplication
    --contamination 25 \      # maximum contamination to include
    --length 50000 \          # minimum genome length (50 Kb)
    --P_ani 0.9 \             # primary clustering ANI threshold (Mash, fast)
    --S_ani 0.95 \            # secondary clustering ANI threshold (precise, species boundary)
    --processors 8
```

**Key output:** `08_drep/drep_output/dereplicated_genomes/` — your final non-redundant MAG set.

> **Expected:** The 13 pooled bins (9 + 4) are de-replicated down to a smaller set. Bins from both samples representing the same species are collapsed to one representative.

---

## Step 9: Final Quality Assessment

### CheckM2 on final MAGs

```bash
conda activate checkm2
checkm2 predict \
    --input 08_drep/drep_output/dereplicated_genomes/ \
    --output-directory 09_final_qc/checkm2 \
    --extension .fa \
    --threads 8
```

### QUAST on final MAGs

```bash
conda activate assembly
for MAG in 08_drep/drep_output/dereplicated_genomes/*.fa; do
    NAME=$(basename $MAG .fa)
    metaquast.py $MAG \
        --output-dir 09_final_qc/quast/${NAME} \
        --threads 4
done
```

### MIMAG quality tiers — the community standard

Use these tiers to describe your MAGs in a paper or report:

| Tier | Completeness | Contamination | When to use |
|---|---|---|---|
| **High quality** | >90% | <5% | Functional annotation, comparative genomics |
| **Medium quality** | ≥50% | <10% | Taxonomic profiling (with caveats) |
| **Low quality** | <50% | — | Taxonomy only, mention limitations |

> Only high-quality MAGs should be used for Bakta functional annotation or any gene-level analysis.

---

## Step 10: Abundance Estimation with CoverM

### Why estimate abundance?

Coverage depth is proportional to the relative abundance of each organism in the original community. This tells you not just *who* is there, but *how much* of each organism is present — and whether that changes between your conditions (probiotic-treated vs. untreated).

```bash
conda activate coverm

coverm genome \
    --bam-files 05_binning/SRR24759598.sorted.bam \
                05_binning/SRR24759616.sorted.bam \
    --genome-fasta-files 08_drep/drep_output/dereplicated_genomes/*.fa \
    --methods relative_abundance trimmed_mean \
    --output-file 10_abundance/coverage_table.tsv \
    --threads 8
```

**Output:** A TSV table with one row per MAG and one column per sample — the starting point for any differential abundance analysis (DESeq2, ANCOM-BC).

---

## Step 11: Taxonomic Assignment with GTDB-Tk

### What is GTDB-Tk?

GTDB-Tk places each MAG into the Genome Taxonomy Database (GTDB) — a phylogenetically consistent bacterial/archaeal classification system that corrects many inconsistencies in traditional NCBI taxonomy.

**How it works:**
1. Identifies conserved marker genes in the MAG
2. Builds a placement tree against thousands of reference genomes
3. Assigns taxonomy based on phylogenetic placement

### Download GTDB database (~66 GB — do this on HPC)

```bash
conda activate gtdbtk
download-db.sh   # automated download script bundled with GTDB-Tk
# Or manually: https://gtdb.ecogenomic.org/downloads
# Set the database path:
export GTDBTK_DATA_PATH=/path/to/gtdbtk_data/
```

### Run GTDB-Tk

```bash
gtdbtk classify_wf \
    --genome_dir 08_drep/drep_output/dereplicated_genomes/ \
    --extension fa \
    --out_dir 11_taxonomy \
    --cpus 16 \
    --skip_ani_screen    # use this flag if the database download check fails
```

**Key output file:** `11_taxonomy/gtdbtk.bac120.summary.tsv`

Each row is one MAG; the `classification` column gives the full lineage:
`d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__...;g__Snodgrassella;s__Snodgrassella alvi`

> **What to watch for:**
> - Classification only to phylum/class level = likely novel lineage, or MAG too incomplete for confident placement
> - Familiar genus names may differ from NCBI — GTDB reclassifies based on phylogenetics, not historical naming

---

## Step 12: Functional Annotation with Bakta

### What does Bakta do?

Bakta annotates all genomic features in each MAG:

| Feature type | What it is |
|---|---|
| CDS (coding sequences) | Protein-coding genes, assigned functions via UniRef90/KEGG/Pfam |
| rRNA genes | 16S, 23S, 5S ribosomal RNA |
| tRNA genes | Transfer RNA |
| ncRNA | Non-coding regulatory RNAs |
| CRISPR arrays | Bacterial adaptive immune elements |
| AMR genes | Antimicrobial resistance (via AMRFinderPlus) |

### Download Bakta database (~30 GB)

```bash
conda activate bakta
bakta_db download --output ~/.bakta_db/ --type full
```

### Run Bakta on each high-quality MAG

```bash
# Run only on high-quality MAGs (>90% complete, <5% contamination)
# Identify these from your CheckM2 output first, then loop

for MAG in 08_drep/drep_output/dereplicated_genomes/*.fa; do
    NAME=$(basename $MAG .fa)
    bakta \
        --db ~/.bakta_db/db \
        --output 12_annotation/${NAME} \
        --prefix ${NAME} \
        --threads 8 \
        --genus $(grep ${NAME} 11_taxonomy/gtdbtk.bac120.summary.tsv | cut -f18) \
        $MAG
done
```

### Key output files per MAG

| File | Contents |
|---|---|
| `{name}.gff3` | Genome annotation (standard format, load into genome browsers) |
| `{name}.tsv` | Tab-separated table of all annotated features |
| `{name}.faa` | Protein sequences of all predicted CDS |
| `{name}.hypotheticals.tsv` | Genes with no known function |

**Biological questions you can now ask:**
- Do any MAGs carry nitrogen fixation genes (*nif* operon)?
- Which MAGs encode antibiotic resistance genes?
- Do MAGs from probiotic-treated bees (IP2) show enriched genes for immune modulation?

---

## Summary: Tools at a Glance

| Step | Tool | conda env | Key input | Key output |
|---|---|---|---|---|
| QC & trim | fastp, MultiQC | `qc` | Raw FASTQ | Trimmed FASTQ |
| Host removal | Bowtie2 | `qc` | Trimmed FASTQ + genome | Clean FASTQ |
| Assembly | MEGAHIT | `assembly` | Clean FASTQ | contigs.fa |
| Assembly QC | QUAST | `assembly` | contigs.fa | HTML report |
| Coverage | Bowtie2, samtools | `qc` | Clean FASTQ + contigs | sorted.bam |
| Binning | MetaBAT2 | `binning` | contigs.fa + depth.txt | bin.N.fa |
| Bin QC | CheckM2 | `checkm2` | Bin FASTA files | quality_report.tsv |
| Refinement | Binette | `binette` | Bin dirs + contigs | Refined bins |
| De-replication | dRep | `drep` | All bins + CheckM2 report | Non-redundant MAGs |
| Abundance | CoverM | `coverm` | BAM + MAG FASTA | Coverage TSV |
| Taxonomy | GTDB-Tk | `gtdbtk` | MAG FASTA | Classification TSV |
| Annotation | Bakta | `bakta` | MAG FASTA | GFF3, FAA, TSV |

---

## Common Pitfalls

**Assembly:**
- Always test your pipeline on a subsampled dataset (100K–500K reads) before committing to a full run
- A low N50 (~2,000 bp) is *normal* for metagenomes — it does not mean something went wrong
- Never use metaSPAdes if you are memory-limited; MEGAHIT is the safer default

**Binning:**
- The `--minContig 1500` flag in MetaBAT2 is important — very short contigs have unreliable TNF signals and produce noisy bins
- More bins ≠ better quality; MetaBAT2's moderate bin count is usually more reliable than CONCOCT's high count

**De-replication:**
- Pool bins from *all* samples before running dRep — the whole point is cross-sample comparison
- The 95% ANI threshold is for species-level work. If you need strain-level distinctions, use 99% ANI and pair with InStrain post-hoc

**Taxonomy:**
- Familiar genus names may look different in GTDB — this is by design, not an error
- Classification only to phylum or class level is biologically interesting (novel lineage!) but requires care in downstream interpretation

**Annotation:**
- Bakta predictions are computational — always verify biologically important hits with BLAST or literature
- Skip Bakta on low-quality MAGs (<50% completeness); fragmented genomes produce fragmented, unreliable annotations

---

## Checkpoint Questions

After completing this workflow, you should be able to answer:

1. You start with 43 million reads and after host removal only 40 million remain. What percentage were removed, and is this concerning?
2. QUAST reports N50 = 2,200 bp for your assembly. A colleague says this means the assembly failed. How do you respond?
3. MetaBAT2 produces 34 bins. Binette refines this to 9. Where did the other 25 go?
4. A MAG has 88% completeness and 7% contamination. What quality tier is it in, and is it suitable for Bakta annotation?
5. Two MAGs from different samples have 96% ANI. What happens to them during dRep, and why?
6. Your GTDB-Tk output assigns one MAG only to phylum level. What are two possible biological explanations?

---

## Next Steps

Once comfortable with this pipeline:

- **Multiple binners → Binette ensemble:** Run MaxBin2, SemiBin, and CONCOCT in addition to MetaBAT2, then feed all four bin sets to Binette. The full GTN workflow does this automatically.
- **Differential abundance analysis:** Take the CoverM coverage table into R (DESeq2 or ANCOM-BC) to test whether MAG abundances differ significantly between treatments.
- **Strain-level resolution:** After de-replication, use InStrain to detect within-species variation across your samples.
- **Workflow automation:** Convert this bash script series into a Snakemake or Nextflow pipeline for reproducibility and HPC scaling.

---

*Tutorial summary written for graduate bioinformatics students. Based on the Galaxy Training Network (GTN) tutorial by Bérénice Batut (CC BY 4.0). Original: https://training.galaxyproject.org/training-material/topics/microbiome/tutorials/mags-building/tutorial.html*
