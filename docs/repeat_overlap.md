# Gene-Repeat Overlap Analyzer

`mandragora repeat-overlap` analyzes overlap between gene coordinates and repeat annotations.

This module is designed for exploring repeat-associated annotation complexity in large, repeat-rich, and unusual plant genomes.

---

## Purpose

The gene-repeat overlap analyzer helps answer questions such as:

| Question | Output |
|---|---|
| How many genes overlap repeats? | `gene_repeat_overlap_summary.tsv` |
| Which genes overlap repeats? | `genes_with_repeat_overlap.tsv` |
| How much of each gene overlaps repeats? | `genes_with_repeat_overlap.tsv` |
| Which repeat classes overlap genes? | `repeat_class_summary.tsv` |
| Are many genes repeat-entangled? | All summary outputs |

This is useful for genomes where repeat content may complicate gene prediction and interpretation.

Examples include:

| Genome type | Why this matters |
|---|---|
| Repeat-rich plant genomes | Genes may be embedded in TE-rich regions |
| Parasitic plant genomes | Gene structure may be unusual or degraded |
| Draft assemblies | Repeats may inflate or confuse gene models |
| Large genomes | Repeat-gene overlap may be widespread |
| Reannotated genomes | Different annotation pipelines may produce different repeat-overlap patterns |

---

## Basic Usage

```bash
python -m mandragora.cli repeat-overlap \
  --genes examples/toy_annotation.gff3 \
  --repeats examples/toy_repeats.bed \
  --outdir results/repeat_overlap
```

Arguments:

| Argument | Description |
|---|---|
| `--genes` / `-g` | Gene annotation file in GFF3/GTF-like format or BED |
| `--repeats` / `-r` | Repeat annotation file in BED format |
| `--outdir` / `-o` | Output directory |

---

## Input

### Gene Annotation

The gene input can be:

```text
GFF3/GTF-like annotation
BED file
```

For GFF3 input, MANDRAGORA extracts features of type:

```text
gene
```

Example GFF3 gene:

```gff3
scaffold_1	toy	gene	101	500	.	+	.	ID=gene1;Name=Toy_gene_1
```

Internally, this becomes:

```bed
scaffold_1	100	500	gene1	0	+
```

### Repeat Annotation

In v0.1, repeat input should be BED.

Recommended six-column BED:

```bed
scaffold_1	180	260	LTR/Gypsy	0	+
scaffold_1	430	470	DNA/TIR	0	+
scaffold_1	720	760	LINE/L1	0	-
scaffold_2	120	160	LTR/Copia	0	+
```

The fourth column is treated as the repeat class/name.

Future versions may support direct RepeatMasker output.

---

## Output Files

Running the repeat-overlap analyzer produces:

```text
results/repeat_overlap/
├── genes_with_repeat_overlap.tsv
├── gene_repeat_overlap_summary.tsv
└── repeat_class_summary.tsv
```

---

## `genes_with_repeat_overlap.tsv`

This file reports repeat overlap per gene.

Example:

```tsv
gene_id	chrom	start	end	strand	gene_length	repeat_overlap_bp	repeat_fraction	repeat_count	repeat_classes
gene1	scaffold_1	100	500	+	400	120	0.300000	2	DNA/TIR,LTR/Gypsy
gene2	scaffold_1	600	900	-	300	40	0.133333	1	LINE/L1
gene3	scaffold_2	100	250	+	150	40	0.266667	1	LTR/Copia
```

Columns:

| Column | Meaning |
|---|---|
| `gene_id` | Gene ID or BED interval name |
| `chrom` | Scaffold/chromosome |
| `start` | Gene start in BED-style coordinates |
| `end` | Gene end in BED-style coordinates |
| `strand` | Gene strand |
| `gene_length` | Gene length in bp |
| `repeat_overlap_bp` | Total bp of the gene overlapping repeats |
| `repeat_fraction` | Fraction of the gene overlapping repeats |
| `repeat_count` | Number of repeat intervals overlapping the gene |
| `repeat_classes` | Repeat classes overlapping the gene |

---

## `gene_repeat_overlap_summary.tsv`

This file summarizes gene-repeat overlap globally.

Example:

```tsv
metric	value
total_genes	3
total_repeats	5
genes_overlapping_repeats	3
percent_genes_overlapping_repeats	100.000000
total_gene_bp	850
total_gene_repeat_overlap_bp	200
```

Metrics:

| Metric | Meaning |
|---|---|
| `total_genes` | Number of genes analyzed |
| `total_repeats` | Number of repeat intervals analyzed |
| `genes_overlapping_repeats` | Number of genes with at least one repeat overlap |
| `percent_genes_overlapping_repeats` | Percentage of genes overlapping repeats |
| `total_gene_bp` | Total bp covered by genes |
| `total_gene_repeat_overlap_bp` | Total bp of gene regions overlapping repeats |

---

## `repeat_class_summary.tsv`

This file summarizes repeat classes that overlap genes.

Example:

```tsv
repeat_class	repeats_overlapping_genes	genes_overlapped	total_overlap_bp
DNA/TIR	1	1	40
LINE/L1	1	1	40
LTR/Copia	1	1	40
LTR/Gypsy	1	1	80
```

Columns:

| Column | Meaning |
|---|---|
| `repeat_class` | Repeat class/name from BED column 4 |
| `repeats_overlapping_genes` | Number of repeat intervals of this class overlapping genes |
| `genes_overlapped` | Number of genes overlapped by this repeat class |
| `total_overlap_bp` | Total bp overlap between this repeat class and genes |

---

## Overlap Logic

MANDRAGORA calculates interval overlap using BED-style coordinates.

Example:

Gene:

```bed
scaffold_1	100	500	gene1	0	+
```

Repeat:

```bed
scaffold_1	180	260	LTR/Gypsy	0	+
```

Overlap:

```text
260 - 180 = 80 bp
```

If multiple repeats overlap the same gene, MANDRAGORA merges overlapping repeat segments before calculating total repeat-covered bp per gene.

This avoids double-counting overlapping repeat intervals for the same gene.

---

## Interpretation for Weird Plant Genomes

High gene-repeat overlap may indicate:

| Observation | Possible interpretation |
|---|---|
| Many genes overlap repeats | Repeat-rich genome or annotation inflation |
| CDS/gene regions heavily overlap repeats | Possible TE-derived or false gene models |
| Specific repeat class dominates gene overlap | TE family may influence gene architecture |
| Some genes are mostly repeat-covered | Candidate suspicious gene models |
| Repeat overlap differs between assemblies | Annotation or assembly pipeline effect |

MANDRAGORA does not automatically classify a gene as false or real. It produces diagnostic tables for interpretation.

---

## Current Limitations

MANDRAGORA v0.1:

| Feature | Status |
|---|---|
| Gene overlap from GFF3 `gene` features | Supported |
| Repeat BED input | Supported |
| RepeatMasker `.out` parsing | Planned |
| RepeatMasker `.gff` parsing | Planned |
| Exon/CDS/intron-specific repeat overlap | Planned |
| Promoter-repeat overlap | Planned |
| Statistical enrichment tests | Planned |
| Visualization | Planned |

---

## Recommended Follow-up Analyses

After running `repeat-overlap`, useful next steps include:

| Follow-up | Purpose |
|---|---|
| Compare repeat-overlapping vs non-overlapping genes | Identify annotation patterns |
| Check protein/domain support | Distinguish real genes from repeat-derived artifacts |
| Compare repeat overlap between assemblies | Detect assembly/annotation effects |
| Run intron analysis | Test whether long introns are repeat-enriched |
| Inspect high-overlap genes manually | Validate suspicious models |

---

## Example Test Run

```bash
python -m mandragora.cli repeat-overlap \
  --genes examples/toy_annotation.gff3 \
  --repeats examples/toy_repeats.bed \
  --outdir results/repeat_overlap
```

Inspect output:

```bash
cat results/repeat_overlap/genes_with_repeat_overlap.tsv
cat results/repeat_overlap/gene_repeat_overlap_summary.tsv
cat results/repeat_overlap/repeat_class_summary.tsv
```
