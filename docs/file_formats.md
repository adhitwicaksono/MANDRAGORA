# File Formats

Project **MANDRAGORA** currently works with simple, standard bioinformatics file formats commonly used in genome annotation, intron analysis, and repeat-overlap analysis.

This document describes the expected input formats for MANDRAGORA v0.1.

---

## Coordinate Conventions

MANDRAGORA internally uses **BED-style coordinates**:

| Format | Coordinate system | Meaning |
|---|---|---|
| GFF3/GTF | 1-based, closed | Both start and end positions are included |
| BED | 0-based, half-open | Start is included, end is excluded |
| MANDRAGORA internal representation | 0-based, half-open | Same as BED |

Example:

A GFF3 feature:

```text
scaffold_1    source    gene    101    500    .    +    .    ID=gene1
```

is internally represented as:

```text
scaffold_1    100    500    gene1    0    +
```

This is important because a feature from GFF3 position `101–500` has length:

```text
500 - 100 = 400 bp
```

---

## Genome FASTA

Genome FASTA is optional for some modules and required for future sequence extraction features.

Example:

```fasta
>scaffold_1
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>scaffold_2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
```

Requirements:

| Requirement | Description |
|---|---|
| FASTA headers | Should match scaffold/chromosome names in annotation files |
| Sequence names | Must be consistent across FASTA, GFF3/GTF, and BED |
| Line wrapping | Any standard FASTA line wrapping is accepted |

Example of matching names:

```text
FASTA: >scaffold_1
GFF3:  scaffold_1
BED:   scaffold_1
```

Example of problematic mismatch:

```text
FASTA: >chr1
GFF3:  scaffold_1
BED:   1
```

MANDRAGORA does not yet automatically repair scaffold-name mismatches.

---

## GFF3 Annotation

GFF3 is the preferred annotation input format for intron analysis.

Minimum useful features:

```text
gene
mRNA or transcript
exon
```

Example:

```gff3
##gff-version 3
scaffold_1	toy	gene	101	500	.	+	.	ID=gene1;Name=Toy_gene_1
scaffold_1	toy	mRNA	101	500	.	+	.	ID=gene1.t1;Parent=gene1
scaffold_1	toy	exon	101	150	.	+	.	ID=gene1.exon1;Parent=gene1.t1
scaffold_1	toy	exon	251	300	.	+	.	ID=gene1.exon2;Parent=gene1.t1
scaffold_1	toy	exon	401	500	.	+	.	ID=gene1.exon3;Parent=gene1.t1
```

MANDRAGORA uses exon coordinates to infer introns.

For the example above:

| Exon | GFF3 coordinates |
|---|---|
| Exon 1 | 101–150 |
| Exon 2 | 251–300 |
| Exon 3 | 401–500 |

The inferred introns are:

| Intron | BED-style coordinates |
|---|---|
| Intron 1 | 150–250 |
| Intron 2 | 300–400 |

---

## GTF Annotation

GTF support is planned but currently limited.

MANDRAGORA v0.1 can read some GTF-like lines, but GFF3 is recommended for now.

Recommended:

```text
Use GFF3 when possible.
```

---

## BED Files

MANDRAGORA reads BED-like interval files with at least three columns:

```text
chrom    start    end
```

Recommended six-column BED:

```text
chrom    start    end    name    score    strand
```

Example:

```bed
scaffold_1	180	260	LTR/Gypsy	0	+
scaffold_1	430	470	DNA/TIR	0	+
scaffold_1	720	760	LINE/L1	0	-
```

Column meanings:

| Column | Meaning |
|---|---|
| `chrom` | Scaffold or chromosome name |
| `start` | 0-based start coordinate |
| `end` | Half-open end coordinate |
| `name` | Feature name, repeat class, or interval ID |
| `score` | Optional score; use `0` if not needed |
| `strand` | `+`, `-`, or `.` |

---

## Repeat BED

For `mandragora repeat-overlap`, repeats should currently be provided as BED.

Example:

```bed
scaffold_1	180	260	LTR/Gypsy	0	+
scaffold_1	430	470	DNA/TIR	0	+
scaffold_2	120	160	LTR/Copia	0	+
```

MANDRAGORA currently treats the fourth BED column as the repeat class/name.

For example:

```text
LTR/Gypsy
DNA/TIR
LTR/Copia
LINE/L1
Simple_repeat
```

Future versions may support direct parsing of:

```text
RepeatMasker .out
RepeatMasker .gff
RepeatMasker .tbl
```

---

## Toy Example Files

The repository includes synthetic toy files for testing:

```text
examples/
├── toy_genome.fa
├── toy_annotation.gff3
├── toy_repeats.bed
└── expected_outputs/
```

These files are not biological data. They are controlled test files designed to verify MANDRAGORA behavior.

---

## Common Format Problems

| Problem | Symptom | Suggested fix |
|---|---|---|
| Scaffold names do not match | No overlaps detected | Make FASTA/GFF3/BED scaffold names identical |
| BED file uses 1-based coordinates | Off-by-one output | Convert BED to true 0-based half-open coordinates |
| GFF3 lacks exon features | No introns inferred | Provide annotation with exon rows |
| GFF3 lacks `ID` or `Parent` attributes | Gene/transcript parsing may fail | Clean annotation attributes first |
| Repeat file has spaces instead of tabs | Parsing may fail or behave unexpectedly | Use tab-delimited BED |
| RepeatMasker `.out` used directly | Not yet supported in v0.1 | Convert RepeatMasker output to BED first |

---

## Recommended Workflow

For v0.1:

```bash
python -m mandragora.cli prepare \
  --annotation examples/toy_annotation.gff3 \
  --repeats examples/toy_repeats.bed \
  --genome examples/toy_genome.fa \
  --outdir results/prepare
```

Then:

```bash
python -m mandragora.cli intron \
  --annotation examples/toy_annotation.gff3 \
  --outdir results/intron
```

And:

```bash
python -m mandragora.cli repeat-overlap \
  --genes examples/toy_annotation.gff3 \
  --repeats examples/toy_repeats.bed \
  --outdir results/repeat_overlap
```
