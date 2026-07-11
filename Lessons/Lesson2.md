# Hybrid Metagenome Assembly and MAG Refinement Tutorial
## Reconstructing High-Quality MAGs Using Illumina and Oxford Nanopore Sequencing

This tutorial walks through the hybrid metagenomic workflow used in the 2019 Nature paper describing *Candidatus Argoarchaeum ethanivorans*. The goal is to recover high-quality metagenome-assembled genomes (MAGs) from an enrichment culture using Illumina short reads and Oxford Nanopore long reads.

---

# Workflow Overview

```text
DNA Extraction
      │
      ▼
Illumina Sequencing
      │
      ▼
Quality Control
      │
      ▼
Metagenome Assembly (SPAdes)
      │
      ▼
Assembly Assessment (metaQUAST)
      │
      ▼
Genome Binning (MaxBin2)
      │
      ▼
Identify MAG of Interest
      │
      ▼
Read Mapping (BBMap)
      │
      ▼
Extract Mapped Reads
      │
      ▼
Reassemble
      │
      ▼
Re-bin
      │
      ▼
Repeat Mapping (Higher Identity)
      │
      ▼
Nanopore Scaffolding (npScarf)
      │
      ▼
Assembly Polishing (Pilon)
      │
      ▼
Genome Quality Assessment (CheckM)
      │
      ▼
Genome Annotation
```

---

# Software Required

- Trimmomatic
- SPAdes
- metaQUAST
- MaxBin2
- BBMap
- Oxford Nanopore Basecaller (Guppy or Dorado; the paper used Albacore)
- Porechop
- npScarf
- Pilon
- CheckM
- Prodigal
- RAST (optional)
- EggNOG
- KEGG
- Pfam

---

# Directory Structure

```text
project/

├── raw_data/
│      illumina/
│      nanopore/
│
├── trimmed/
├── assembly/
├── bins/
├── refinement/
├── nanopore/
├── polished/
└── annotation/
```

---

# Step 1 — Illumina Quality Control

Remove adapters and low-quality bases.

```bash
trimmomatic PE \
sample_R1.fastq.gz \
sample_R2.fastq.gz \
R1_paired.fastq.gz \
R1_unpaired.fastq.gz \
R2_paired.fastq.gz \
R2_unpaired.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
LEADING:20 \
TRAILING:20 \
SLIDINGWINDOW:4:20 \
MINLEN:36
```

Output

```
R1_paired.fastq.gz
R2_paired.fastq.gz
```

These reads are used for assembly.

---

# Step 2 — Assemble the Metagenome

Run SPAdes in metagenomic mode.

```bash
spades.py \
--meta \
-1 R1_paired.fastq.gz \
-2 R2_paired.fastq.gz \
-o assembly/
```

The paper used multiple k-mers.

Equivalent command:

```bash
spades.py \
--meta \
-k 21,33,55,77,99,127 \
-1 R1_paired.fastq.gz \
-2 R2_paired.fastq.gz \
-o assembly/
```

Output

```
assembly/scaffolds.fasta
```

---

# Step 3 — Evaluate Assembly Quality

Use metaQUAST.

```bash
metaquast.py \
assembly/scaffolds.fasta \
-o metaquast_output
```

Important statistics

- Total assembly length
- Number of contigs
- N50
- Largest contig
- GC%

Although the Nature paper does not report the N50, they used metaQUAST to inspect the assembly.

---

# Step 4 — Genome Binning

Cluster contigs into draft genomes.

```bash
run_MaxBin.pl \
-contig assembly/scaffolds.fasta \
-abund abundance.txt \
-out bins/maxbin
```

Output

```
bin.001.fasta
bin.002.fasta
bin.003.fasta
```

Each bin represents a draft MAG.

---

# Step 5 — Identify the MAG of Interest

Extract the 16S rRNA gene.

```bash
rnammer \
-S bac \
-m tsu,ssu \
< bin001.fasta
```

Compare against SILVA.

Alternatively

```bash
blastn \
-query 16S.fasta \
-db SILVA_database \
-out results.txt
```

Suppose

```
Bin001 = archaeon

Bin002 = sulfate reducer
```

---

# Step 6 — First Read Recruitment

Map ALL Illumina reads back to the archaeal MAG.

The paper used

```
Minimum identity = 90%
```

BBMap

```bash
bbmap.sh \
ref=bin001.fasta \
in1=R1_paired.fastq.gz \
in2=R2_paired.fastq.gz \
outm=mapped.sam \
minid=0.90
```

This recruits reads belonging to the archaeon.

---

# Step 7 — Extract Mapped Reads

Convert SAM to BAM.

```bash
samtools view -Sb mapped.sam > mapped.bam
```

Extract mapped read pairs.

```bash
reformat.sh \
in=mapped.bam \
out1=archaea_R1.fastq.gz \
out2=archaea_R2.fastq.gz
```

Now you have reads almost exclusively from the archaeon.

---

# Step 8 — Reassemble the Archaeal Reads

Assemble only the recruited reads.

```bash
spades.py \
--meta \
-1 archaea_R1.fastq.gz \
-2 archaea_R2.fastq.gz \
-o refinement_round1
```

Advantages

- fewer contaminating reads
- longer contigs
- improved coverage
- fewer chimeras

---

# Step 9 — Re-bin

Run MaxBin again.

```bash
run_MaxBin.pl \
-contig refinement_round1/scaffolds.fasta \
-abund abundance.txt \
-out refinement_bins
```

The MAG should now be cleaner.

---

# Step 10 — Second Read Recruitment

Repeat the process.

This time increase mapping stringency.

The paper used

```
97% identity
```

```bash
bbmap.sh \
ref=refinement_bins/bin001.fasta \
in1=R1_paired.fastq.gz \
in2=R2_paired.fastq.gz \
outm=round2.sam \
minid=0.97
```

Why?

The first refinement produced a much cleaner genome.

Now only reads almost identical to the MAG are retained.

---

# Step 11 — Extract Reads Again

```bash
samtools view -Sb round2.sam > round2.bam

reformat.sh \
in=round2.bam \
out1=round2_R1.fastq.gz \
out2=round2_R2.fastq.gz
```

---

# Step 12 — Final Reassembly

```bash
spades.py \
--meta \
-1 round2_R1.fastq.gz \
-2 round2_R2.fastq.gz \
-o final_archaea
```

At this point you have a refined MAG assembled solely from reads recruited to that genome.

---

# Step 13 — Hybrid Scaffolding

Use Nanopore reads to connect Illumina contigs.

```bash
npScarf \
-contig final_archaea/scaffolds.fasta \
-long nanopore_trimmed.fastq \
-prefix scaffolded
```

Output

```
scaffolded.fasta
```

Long reads bridge repetitive regions that Illumina cannot resolve.

---

# Step 14 — Polish the Genome

Map Illumina reads.

```bash
bwa index scaffolded.fasta

bwa mem \
scaffolded.fasta \
R1_paired.fastq.gz \
R2_paired.fastq.gz \
> polish.sam
```

Convert

```bash
samtools view -Sb polish.sam > polish.bam

samtools sort polish.bam > polish.sorted.bam

samtools index polish.sorted.bam
```

Run Pilon.

```bash
java -jar pilon.jar \
--genome scaffolded.fasta \
--frags polish.sorted.bam \
--output polished
```

Pilon corrects

- SNPs
- insertions
- deletions
- local assembly errors

---

# Step 15 — Assess Genome Quality

```bash
checkm lineage_wf \
polished/ \
checkm_output/
```

Typical output

```
Completeness: 93%

Contamination: 2%

Strain heterogeneity: 0%
```

These values indicate a high-quality MAG.

---

# Step 16 — Gene Prediction

```bash
prodigal \
-i polished.fasta \
-a proteins.faa \
-d genes.fna \
-p single
```

---

# Step 17 — Functional Annotation

Example with EggNOG.

```bash
emapper.py \
-i proteins.faa \
-o annotation
```

Other annotation databases

- KEGG
- Pfam
- RAST

---

# Final Workflow

```text
Illumina Reads
      │
      ▼
Quality Control
      │
      ▼
SPAdes Assembly
      │
      ▼
MaxBin
      │
      ▼
Initial MAG
      │
      ▼
Map Reads (90%)
      │
      ▼
Extract Reads
      │
      ▼
Reassemble
      │
      ▼
Re-bin
      │
      ▼
Map Reads (97%)
      │
      ▼
Extract Reads
      │
      ▼
Final Reassembly
      │
      ▼
Nanopore Sequence Data
      │
      ▼
Hybrid Scaffolding (npScarf)
      │
      ▼
Pilon Polishing
      │
      ▼
CheckM
      │
      ▼
Final High-Quality MAG
      │
      ▼
Functional Annotation
```

# Key takeaways

- Initial metagenome assemblies often contain contamination and fragmented genomes.
- Iterative read recruitment (90% then 97% identity) enriches reads belonging to the target organism and improves MAG quality.
- Reassembling recruited reads reduces chimeric contigs and increases completeness.
- Nanopore long reads are used to scaffold Illumina contigs, increasing genome continuity.
- Illumina reads are then used to polish the assembly, correcting errors introduced by long-read sequencing.
- The final product is a high-quality metagenome-assembled genome suitable for downstream metabolic reconstruction and comparative genomics.

# Improvemets with long reads
Nanopore reads
      │
      ▼
metaFlye (--meta)
      │
      ▼
Medaka
      │
      ▼
Map Illumina reads
      │
      ▼
Polypolish (or Pilon)
      │
      ▼
MetaBAT2
      │
      ▼
CheckM2
      │
      ▼
Annotation
