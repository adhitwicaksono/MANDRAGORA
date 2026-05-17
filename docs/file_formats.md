cat > docs/file_formats.md <<'EOF'
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
