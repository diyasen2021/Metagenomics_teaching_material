# Nanopore Metagenome Assembly Pipeline — Older Chemistry (R9.4.1)

**Tool:** Flye (metaFlye mode)  
**Designed for:** High-error-rate reads (~90–95% accuracy, NOT Q20+)  
**Organism:** Drosophila melanogaster gut microbiome (or any insect gut)

---

## Why This Pipeline Differs from Illumina

Older Nanopore chemistry (R9.4.1, R9.5) produces reads with:

- **Mean accuracy:** ~90–95% (vs >99% for modern R10.4.1 / Q20+)
- **Error profile:** predominantly insertions and deletions (indels), not substitutions — very different from Illumina errors
- **Read length:** 1–100 kb (median ~8–15 kb depending on prep)
- **Coverage bias:** homopolymer runs are difficult to resolve

**Consequences for the pipeline:**
1. Use `--nano-raw` in Flye (NOT `--nano-hq` — that flag assumes Q20+)
2. Polish aggressively: multiple Racon rounds then Medaka
3. Medaka model **must** match the flowcell chemistry — wrong model = worse results
4. Quality filtering removes short/low-quality reads BEFORE assembly
5. No adapter trimming needed (Guppy/MinKNOW already strips adapters during basecalling)

---

## Medaka Model Quick Reference (Older Chemistry)

| Basecaller & Chemistry | Model |
|---|---|
| Guppy 3.x, R9.4.1 HAC | `r941_min_high_g303` |
| Guppy 4.x–5.x, R9.4.1 HAC | `r941_min_high_g360` ← most common |
| Guppy ≥5.0, R9.4.1 SUP | `r941_min_sup_g507` |
| Guppy ≥5.0, PromethION R9.4.1 SUP | `r941_prom_sup_g507` |

> Check your Guppy log or MinKNOW run report for the exact basecaller version. If unsure, `r941_min_high_g360` is the safest default for older R9.4.1 data.

---

## Environment Installation

```bash
conda create -n nanopore_qc  -y nanoplot filtlong
conda create -n flye         -y flye
conda create -n assembly     -y quast
conda create -n polishing    -y minimap2 samtools racon
conda create -n medaka       -y medaka
conda create -n binning      -y metabat2
conda create -n checkm2      -y checkm2
conda create -n prokka       -y prokka
```

> **Medaka GPU note:** Medaka runs on CPU but is very slow (~8–12 hours for large assembly). On GPU (CUDA): 30–60 minutes. If your HPC has GPUs, request one.

---

## Configuration

Edit these variables for your dataset before running:

```bash
SAMPLE="drosophila_gut"
RAW_FASTQ="raw_reads/${SAMPLE}.fastq.gz"
THREADS=16
MIN_READ_LEN=1000           # discard reads shorter than 1 kb
MIN_READ_QUAL=8             # discard reads below Q8 (lenient — older reads are noisy)
MEDAKA_MODEL="r941_min_high_g360"   # CHANGE to match your chemistry
RACON_ROUNDS=3              # number of Racon polishing iterations

mkdir -p 00_raw 01_qc 02_filtered 03_assembly 04_polishing/racon \
         04_polishing/medaka 05_quast 06_binning 07_annotation logs
```

---

## Step 0 — Basecalling (skip if you already have FASTQ)

If you have raw FAST5 files not yet basecalled, use Guppy:

```bash
guppy_basecaller \
    --input_path fast5/ \
    --save_path basecalled/ \
    --config dna_r9.4.1_450bps_hac.cfg \   # HAC = high-accuracy config
    --device cuda:0 \                        # GPU — remove if CPU only
    --recursive \
    --compress_fastq

# Merge all passing reads into one file
cat basecalled/pass/*.fastq.gz > raw_reads/${SAMPLE}.fastq.gz
```

> **Note:** Guppy is proprietary (Oxford Nanopore) and NOT available on conda. Download from the [Nanopore Community portal](https://community.nanoporetech.com/downloads). Modern open-source alternative: Dorado — but use older models if working with R9.4.1 data.

---

## Step 1 — Raw Read QC with NanoPlot

NanoPlot gives you a detailed HTML report on read length, quality, and yield — essential before filtering so you know what you are working with.

```bash
conda activate nanopore_qc

NanoPlot \
    --fastq  ${RAW_FASTQ} \
    --outdir 01_qc/raw_nanoplot \
    --threads ${THREADS} \
    --plots dot \
    --N50 \
    --title "${SAMPLE} — Raw reads"
```

**Key metrics to check in the HTML report:**

| Metric | Expected for older R9.4.1 |
|---|---|
| Mean read quality | Q8–Q12 |
| N50 read length | >5 kb for good metagenomic assemblies |
| Total yield | >50× coverage of expected community size |
| Read length histogram | Long tail — not truncated at short lengths |

---

## Step 2 — Read Filtering with Filtlong

Filtlong keeps the longest, highest-quality reads up to a target total bases. For older Nanopore reads this is more critical than for Illumina because short reads (<1 kb) have poor signal-to-noise, and very long reads with low quality drag down assembly accuracy.

```bash
filtlong \
    --min_length   ${MIN_READ_LEN} \
    --min_mean_q   ${MIN_READ_QUAL} \
    --target_bases 10000000000 \
    ${RAW_FASTQ} | gzip > 02_filtered/${SAMPLE}_filtered.fastq.gz
```

> **`--target_bases`:** Set to ~100× your expected total community genome size. For a simple gut community (~50 Mb total), 5 Gb is sufficient. For a complex soil community, use 20–50 Gb.

Check filtered output:

```bash
NanoPlot \
    --fastq  02_filtered/${SAMPLE}_filtered.fastq.gz \
    --outdir 01_qc/filtered_nanoplot \
    --threads ${THREADS} \
    --N50 \
    --title "${SAMPLE} — Filtered reads"

echo "Reads before filtering:"
zcat ${RAW_FASTQ} | awk 'NR%4==1' | wc -l

echo "Reads after filtering:"
zcat 02_filtered/${SAMPLE}_filtered.fastq.gz | awk 'NR%4==1' | wc -l
```

---

## Step 3 — Metagenome Assembly with Flye

Flye uses a repeat graph approach that handles the high error rate of older reads well.

### Critical flag: `--nano-raw` vs `--nano-hq`

| Flag | When to use |
|---|---|
| `--nano-raw` | R9.4.1 reads, ~90–95% accuracy ← **USE THIS** |
| `--nano-hq` | R10.4.1 / Q20+ reads, >97% accuracy |
| `--nano-corr` | Pre-corrected reads (e.g. corrected by Canu first) |

Using `--nano-hq` on R9.4.1 reads will make Flye expect higher quality than exists — resulting in a fragmented or failed assembly.

```bash
conda activate flye

flye \
    --nano-raw 02_filtered/${SAMPLE}_filtered.fastq.gz \
    --out-dir  03_assembly \
    --threads  ${THREADS} \
    --meta \
    --min-overlap 2000 \
    --keep-haplotypes

cp 03_assembly/assembly.fasta 03_assembly/${SAMPLE}_assembly.fasta
```

**`--meta`** enables metagenomic mode: adjusts minimum overlap parameters for uneven coverage and allows assembly of low-coverage organisms that would otherwise fail.

**`--min-overlap 2000`** is reduced from the default 3000 because high error rates mean overlaps look less clean — the default can miss valid overlaps in noisy older reads.

**`--keep-haplotypes`** is important for metagenomes — do not collapse strain variants, as different strains are biologically meaningful.

**Flye outputs:**

| File | Contents |
|---|---|
| `assembly.fasta` | Final assembled contigs |
| `assembly_info.txt` | Per-contig stats (coverage, length, circular?) |
| `assembly_graph.gfa` | Assembly graph — visualise in Bandage |
| `flye.log` | Detailed log |

Check for circular contigs — these are likely complete genomes:

```bash
echo "Total contigs assembled:"
grep -c ">" 03_assembly/assembly.fasta

echo "Circular contigs (likely complete genomes):"
grep "Y" 03_assembly/assembly_info.txt
```

---

## Step 4 — Assembly QC with QUAST (Pre-polishing)

Run QUAST before polishing so you can compare before/after metrics.

```bash
conda activate assembly

metaquast.py \
    03_assembly/assembly.fasta \
    --output-dir 05_quast/pre_polish \
    --threads    ${THREADS} \
    --min-contig 1000
```

**What to expect from a Nanopore assembly (very different from Illumina):**

| Metric | Illumina metagenome | Older Nanopore metagenome |
|---|---|---|
| N50 | ~2,000 bp | ~200,000–800,000 bp |
| Total contigs | 50,000–500,000 | 200–5,000 |
| Largest contig | ~100–500 kb | Can be several Mb |
| Circular contigs | Rare | Common for dominant organisms |

---

## Step 5 — Polishing: Racon (3 rounds)

Racon uses the raw reads to correct errors in the assembly. For older R9.4.1 reads, 3 rounds of Racon followed by 1 round of Medaka is the established best practice.

**Why multiple Racon rounds?** Each round corrects errors from the previous assembly but with diminishing returns. Three rounds captures >95% of correctable errors before Medaka handles the rest.

**Why not just use Medaka directly?** Medaka works best on assemblies already partially corrected by Racon. Raw Flye → Medaka leaves more errors than Racon ×3 → Medaka.

```bash
conda activate polishing

CURRENT_ASSEMBLY="03_assembly/assembly.fasta"

for ROUND in 1 2 3; do
    echo "--- Racon round ${ROUND} ---"

    # Map reads to current assembly
    minimap2 \
        -ax map-ont \
        -t ${THREADS} \
        ${CURRENT_ASSEMBLY} \
        02_filtered/${SAMPLE}_filtered.fastq.gz \
        > 04_polishing/racon/round${ROUND}_mapping.sam

    # Run Racon
    racon \
        --threads ${THREADS} \
        --quality-threshold 5 \
        02_filtered/${SAMPLE}_filtered.fastq.gz \
        04_polishing/racon/round${ROUND}_mapping.sam \
        ${CURRENT_ASSEMBLY} \
        > 04_polishing/racon/${SAMPLE}_racon_r${ROUND}.fasta

    CURRENT_ASSEMBLY="04_polishing/racon/${SAMPLE}_racon_r${ROUND}.fasta"

    # Clean up large SAM file
    rm 04_polishing/racon/round${ROUND}_mapping.sam

    echo "Round ${ROUND} complete — assembly: ${CURRENT_ASSEMBLY}"
done
```

> **`--quality-threshold 5`** is lower than the default to be more tolerant of older noisy reads — at the default threshold Racon may discard too many alignments.

---

## Step 6 — Polishing: Medaka

Medaka is Oxford Nanopore's neural network polisher trained on actual Nanopore error profiles. It corrects systematic errors that Racon cannot — particularly homopolymer runs.

### Model selection — the most important parameter

The model **must** match the flowcell chemistry AND the basecaller version used. Using the wrong model introduces more errors than it corrects.

| Guppy version & chemistry | Medaka model |
|---|---|
| Guppy 3.x, R9.4.1 HAC | `r941_min_high_g303` |
| Guppy 4.x–5.x, R9.4.1 HAC | `r941_min_high_g360` |
| Guppy ≥5.0, R9.4.1 SUP | `r941_min_sup_g507` |
| Guppy ≥5.0, PromethION R9.4.1 SUP | `r941_prom_sup_g507` |

```bash
conda activate medaka

medaka_haploid_variant \
    -i 02_filtered/${SAMPLE}_filtered.fastq.gz \
    -r ${CURRENT_ASSEMBLY} \
    -o 04_polishing/medaka \
    -t ${THREADS} \
    -m ${MEDAKA_MODEL}

POLISHED="04_polishing/medaka/consensus.fasta"
cp ${POLISHED} ${SAMPLE}_final_assembly.fasta

echo "Final polished assembly: ${SAMPLE}_final_assembly.fasta"
echo "Contig count:"
grep -c ">" ${SAMPLE}_final_assembly.fasta
```

---

## Step 7 — Assembly QC Post-Polishing

```bash
conda activate assembly

metaquast.py \
    ${SAMPLE}_final_assembly.fasta \
    --output-dir 05_quast/post_polish \
    --threads    ${THREADS} \
    --min-contig 1000
```

Polishing should improve gene completeness and genome fraction. Contig count and N50 should remain similar — polishing fixes bases, not structure.

---

## Step 8 — Binning with MetaBAT2

The binning approach is similar to Illumina but with key differences: coverage is calculated from long reads using minimap2 (not bowtie2), and some bins may be complete circular genomes from the assembly.

```bash
conda activate polishing

# Map reads back to polished assembly
minimap2 \
    -ax map-ont \
    -t ${THREADS} \
    ${SAMPLE}_final_assembly.fasta \
    02_filtered/${SAMPLE}_filtered.fastq.gz | \
samtools sort \
    -@ ${THREADS} \
    -o 06_binning/${SAMPLE}.sorted.bam

samtools index 06_binning/${SAMPLE}.sorted.bam

conda activate binning

# Generate depth file
jgi_summarize_bam_contig_depths \
    --outputDepth 06_binning/${SAMPLE}_depth.txt \
    06_binning/${SAMPLE}.sorted.bam

# Run MetaBAT2
mkdir -p 06_binning/bins
metabat2 \
    -i ${SAMPLE}_final_assembly.fasta \
    -a 06_binning/${SAMPLE}_depth.txt \
    -o 06_binning/bins/bin \
    --minContig 1000 \
    --maxEdges  500 \
    --seed 42 \
    --numThreads ${THREADS} \
    -v
```

> **`--minContig 1000`** is lower than the Illumina default of 1500 because Nanopore produces far fewer total contigs — every contig is valuable and worth including in binning.

---

## Step 9 — Bin Quality with CheckM2

```bash
conda activate checkm2

checkm2 predict \
    --input     06_binning/bins/ \
    --output-directory 06_binning/checkm2 \
    --extension .fa \
    --threads   ${THREADS}

cat 06_binning/checkm2/quality_report.tsv | column -t
```

> For Nanopore assemblies, expect **higher completeness** and **lower contamination** than equivalent Illumina bins — because long reads produce less fragmentation and fewer chimeric contigs during binning. High-quality Nanopore MAGs (>95% complete, <2% contamination) are common for dominant community members.

---

## Step 10 — Functional Annotation with Prokka

For older Nanopore assemblies, Prokka is preferred over Bakta because the polished assembly may still contain occasional indel errors causing frameshifts. Prokka's `--metagenome` mode is more tolerant of partial ORFs.

```bash
conda activate prokka

for BIN in 06_binning/bins/*.fa; do
    BINNAME=$(basename ${BIN} .fa)

    COMPLETENESS=$(grep "${BINNAME}" 06_binning/checkm2/quality_report.tsv | awk '{print $2}')
    echo "Bin: ${BINNAME}  Completeness: ${COMPLETENESS}%"

    if (( $(echo "${COMPLETENESS} >= 75" | bc -l) )); then
        prokka \
            --outdir    07_annotation/${BINNAME} \
            --prefix    ${BINNAME} \
            --metagenome \
            --rfam \
            --cpus      ${THREADS} \
            --kingdom   Bacteria \
            ${BIN}
        echo "Annotated: ${BINNAME}"
    else
        echo "Skipped (low completeness): ${BINNAME}"
    fi
done
```

---

## Key Differences from Illumina Pipeline — Summary

| Step | Illumina pipeline | Older Nanopore (this pipeline) |
|---|---|---|
| QC tool | fastp | NanoPlot |
| Trim/filter | fastp (adapters) | filtlong (length + quality) |
| Adapter removal | Required | Done by Guppy already |
| Assembler | MEGAHIT | Flye `--meta --nano-raw` |
| Polishing | None needed | Racon ×3 + Medaka (essential) |
| Medaka model | n/a | Must match chemistry/basecaller |
| Read mapper | bowtie2 | minimap2 `-ax map-ont` |
| Expected N50 | ~2,000 bp | ~200,000–800,000 bp |
| Expected contigs | 50,000–500,000 | 200–5,000 |
| Circular contigs | Rare | Common for dominant organisms |
| Binning min length | 1,500 bp | 1,000 bp |
| Annotation tool | Bakta | Prokka `--metagenome` |

---

## Common Pitfalls for Older Nanopore Data

**1. Using `--nano-hq` instead of `--nano-raw` in Flye**
Assembly will be fragmented or fail entirely. Always match the flag to your actual read accuracy.

**2. Wrong Medaka model**
Using the wrong model introduces new errors. Check your Guppy log. When uncertain, `r941_min_high_g360` is the safest fallback for R9.4.1.

**3. Skipping Racon before Medaka**
Medaka alone on a raw Flye assembly leaves more uncorrected errors than the Racon ×3 → Medaka combination.

**4. Using bowtie2 instead of minimap2 for read mapping**
bowtie2 is a short-read aligner — it cannot align reads >~1 kb. Always use `minimap2 -ax map-ont` for Nanopore reads.

**5. Setting `--target_bases` too low in filtlong**
If you filter down to too few total bases, low-abundance organisms lose all coverage and disappear from the assembly. Aim for >50× coverage of your expected total community size.

**6. Expecting Illumina-like contig counts**
Nanopore assemblies have far fewer, far longer contigs. 300 contigs covering 5 Mb is excellent. Do not compare N50 or contig count directly between Illumina and Nanopore assemblies of the same sample.

---

*Pipeline designed for older Nanopore chemistry (R9.4.1). For modern R10.4.1 / Q20+ reads, use `--nano-hq` in Flye, reduce Racon to 1 round, and select the appropriate Medaka or Dorado polishing model.*
