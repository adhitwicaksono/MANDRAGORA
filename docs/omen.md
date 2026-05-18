# Gene Omen Scoring

`mandragora omen` scores genes for suspicious architecture and repeat-associated weirdness.

This module is designed as a diagnostic prioritization tool for strange, repeat-rich, parasitic, or poorly resolved plant genomes. It does **not** classify genes as real or false. Instead, it helps identify gene models that deserve closer inspection.

---

## Purpose

The omen module helps answer questions such as:

| Question | Output |
|---|---|
| Which genes look structurally suspicious? | `gene_omen_scores.tsv` |
| Which genes have high repeat overlap? | `gene_omen_scores.tsv` |
| Which short or single-exon genes overlap repeats? | `gene_omen_scores.tsv` |
| Which genes should be manually inspected first? | `gene_omen_scores.tsv` |
| How many genes fall into each warning category? | `gene_omen_summary.tsv` |

This is useful for genomes where gene prediction may be complicated by:

| Factor | Why it matters |
|---|---|
| Large genome size | Gene spans may be inflated |
| Repeat-rich regions | Repeats may overlap or obscure gene models |
| Parasitic evolution | Gene loss, pseudogenization, and unusual gene architecture may occur |
| Fragmented assemblies | Gene models may be incomplete or incorrectly joined |
| Draft annotations | Predicted genes may need triage before biological interpretation |

---

## Basic Usage

```bash
python -m mandragora.cli omen \
  --annotation examples/toy_annotation.gff3 \
  --repeats examples/toy_repeats.bed \
  --outdir results/omen
```

Arguments:

| Argument | Description |
|---|---|
| `--annotation` / `-a` | Input genome annotation in GFF3/GTF-like format |
| `--repeats` / `-r` | Input repeat annotation in BED format |
| `--outdir` / `-o` | Output directory |

---

## Optional Thresholds

The omen module uses configurable thresholds.

```bash
python -m mandragora.cli omen \
  --annotation examples/toy_annotation.gff3 \
  --repeats examples/toy_repeats.bed \
  --outdir results/omen \
  --short-gene 300 \
  --long-gene 100000 \
  --very-long-gene 500000 \
  --long-intron 100000 \
  --very-long-intron 500000 \
  --repeat-warning 0.25 \
  --repeat-high 0.50 \
  --repeat-severe 0.75
```

Thresholds:

| Option | Default | Meaning |
|---|---:|---|
| `--short-gene` | `300` | Gene length below this value receives a short-gene warning |
| `--long-gene` | `100000` | Gene length at or above this value receives a long-gene warning |
| `--very-long-gene` | `500000` | Gene length at or above this value receives a stronger warning |
| `--long-intron` | `100000` | Maximum intron length at or above this value receives a long-intron warning |
| `--very-long-intron` | `500000` | Maximum intron length at or above this value receives a stronger warning |
| `--repeat-warning` | `0.25` | Repeat fraction at or above this value receives a warning |
| `--repeat-high` | `0.50` | Repeat fraction at or above this value receives a stronger warning |
| `--repeat-severe` | `0.75` | Repeat fraction at or above this value receives the strongest repeat warning |

---

## Input

The omen module currently requires:

```text
annotation.gff3
repeats.bed
```

The annotation file is used to extract:

| Feature | Use |
|---|---|
| `gene` | Gene coordinates and gene length |
| `mRNA` / `transcript` | Transcript structure |
| `exon` | Intron inference |

The repeat BED file is used to calculate:

| Feature | Use |
|---|---|
| repeat overlap bp | How many bp of each gene overlap repeats |
| repeat fraction | Fraction of gene length overlapping repeats |
| repeat classes | Which repeat classes overlap each gene |

---

## Output Files

Running the omen module produces:

```text
results/omen/
├── gene_omen_scores.tsv
└── gene_omen_summary.tsv
```

---

## `gene_omen_scores.tsv`

This file reports gene-level omen scores.

Example:

```tsv
gene_id	chrom	start	end	strand	gene_length	exon_count	intron_count	max_intron_length	repeat_overlap_bp	repeat_fraction	repeat_count	repeat_classes	omen_score	omen_level	flags
gene1	scaffold_1	100	500	+	400	3	2	100	120	0.300000	2	DNA/TIR,LTR/Gypsy	2	MODERATE	repeat_fraction_ge_0.25
gene2	scaffold_1	600	900	-	300	2	1	150	40	0.133333	1	LINE/L1	0	NONE	.
gene3	scaffold_2	100	250	+	150	1	0	0	40	0.266667	1	LTR/Copia	4	HIGH	short_gene,repeat_fraction_ge_0.25,single_exon_repeat_overlap
```

Columns:

| Column | Meaning |
|---|---|
| `gene_id` | Gene ID |
| `chrom` | Scaffold/chromosome |
| `start` | Gene start in BED-style coordinates |
| `end` | Gene end in BED-style coordinates |
| `strand` | Gene strand |
| `gene_length` | Gene length in bp |
| `exon_count` | Maximum exon count observed for the gene |
| `intron_count` | Maximum intron count observed for the gene |
| `max_intron_length` | Longest inferred intron for the gene |
| `repeat_overlap_bp` | Total bp of gene region overlapping repeats |
| `repeat_fraction` | Fraction of gene length overlapping repeats |
| `repeat_count` | Number of repeat intervals overlapping the gene |
| `repeat_classes` | Repeat classes overlapping the gene |
| `omen_score` | Numeric warning score |
| `omen_level` | Qualitative warning level |
| `flags` | Reasons contributing to the omen score |

---

## `gene_omen_summary.tsv`

This file summarizes how many genes fall into each warning category.

Example:

```tsv
metric	value
total_genes	3
none_omen_genes	1
low_omen_genes	0
moderate_omen_genes	1
high_omen_genes	1
severe_omen_genes	0
```

Metrics:

| Metric | Meaning |
|---|---|
| `total_genes` | Total number of genes analyzed |
| `none_omen_genes` | Genes with no warning |
| `low_omen_genes` | Genes with low warning level |
| `moderate_omen_genes` | Genes with moderate warning level |
| `high_omen_genes` | Genes with high warning level |
| `severe_omen_genes` | Genes with severe warning level |

---

## Omen Score Logic

The omen score is additive.

Current warning signals include:

| Signal | Score effect |
|---|---:|
| Short gene | `+1` |
| Long gene | `+1` |
| Very long gene | `+2` |
| Long intron | `+1` |
| Very long intron | `+2` |
| Repeat fraction >= 0.25 | `+2` |
| Repeat fraction >= 0.50 | `+3` |
| Repeat fraction >= 0.75 | `+4` |
| Single-exon gene with repeat overlap | `+1` |

The final numeric score is converted into a warning level:

| Omen score | Omen level |
|---:|---|
| `0` | `NONE` |
| `1` | `LOW` |
| `2-3` | `MODERATE` |
| `4-5` | `HIGH` |
| `>=6` | `SEVERE` |

---

## Interpretation

The omen score is a **triage score**, not a biological truth label.

A high omen score may indicate:

| Observation | Possible interpretation |
|---|---|
| Very short gene | Fragment, pseudogene, annotation artifact, or real small gene |
| Very long gene span | Repeat expansion, annotation inflation, or genuine large gene |
| Very long intron | Repeat-rich intron, poor gene prediction, or unusual biology |
| High repeat fraction | TE-associated gene, false gene model, or repeat-adjacent gene |
| Single-exon repeat-overlapping gene | Possible TE-derived prediction or suspicious annotation |

MANDRAGORA does not automatically remove or reject genes. It helps prioritize gene models for manual inspection.

---

## Example Interpretation

In the toy dataset:

| Gene | Pattern | Omen level |
|---|---|---|
| `gene1` | Moderate repeat overlap | `MODERATE` |
| `gene2` | Low repeat overlap, otherwise ordinary | `NONE` |
| `gene3` | Short, single-exon, repeat-overlapping | `HIGH` |

This does not mean `gene3` is false. It means `gene3` should be inspected before being used for biological interpretation.

---

## Current Limitations

MANDRAGORA v0.2 omen scoring is an early prototype.

Current limitations include:

- It does not check protein/domain support.
- It does not check transcript or RNA-seq evidence.
- It does not distinguish repeat overlap in exon, CDS, intron, UTR, or promoter regions yet.
- It does not yet use BUSCO/ortholog support.
- It does not yet include functional annotation support.
- Thresholds are heuristic and should be adjusted for each genome.

---

## Recommended Follow-up Analyses

After running `mandragora omen`, useful next steps include:

| Follow-up | Purpose |
|---|---|
| Inspect high/severe omen genes manually | Validate suspicious gene models |
| Compare against protein/domain annotation | Check biological support |
| Compare against RNA-seq evidence | Confirm expression support |
| Run repeat-overlap analysis | Inspect repeat classes and overlap fractions |
| Run intron analysis | Check whether long introns drive the warning |
| Compare between assemblies | Detect assembly/annotation pipeline effects |

---

## Example Test Run

```bash
python -m mandragora.cli omen \
  --annotation examples/toy_annotation.gff3 \
  --repeats examples/toy_repeats.bed \
  --outdir results/omen
```

Inspect output:

```bash
cat results/omen/gene_omen_scores.tsv
cat results/omen/gene_omen_summary.tsv
```
