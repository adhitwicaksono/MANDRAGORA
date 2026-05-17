from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List

from mandragora.utils import (
    Interval,
    ensure_outdir,
    extract_gene_intervals_from_gff3,
    format_float,
    merge_intervals,
    overlap_bp,
    read_bed,
    write_tsv,
)


def load_gene_intervals(gene_path: str | Path) -> List[Interval]:
    """
    Load gene intervals from GFF3/GTF-like annotation or BED.

    For v0.1:
    - GFF3/GTF input: extracts feature type 'gene'
    - BED input: uses intervals directly
    """
    gene_path = Path(gene_path)
    suffixes = {suffix.lower() for suffix in gene_path.suffixes}

    if ".gff" in suffixes or ".gff3" in suffixes or ".gtf" in suffixes:
        return extract_gene_intervals_from_gff3(gene_path)

    return read_bed(gene_path)


def analyze_gene_repeat_overlap(
    gene_path: str | Path,
    repeat_bed_path: str | Path,
) -> tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]]]:
    """
    Analyze overlap between gene intervals and repeat intervals.

    Repeats are expected as BED-like intervals:
      chrom, start, end, repeat_class/name, score, strand
    """
    genes = load_gene_intervals(gene_path)
    repeats = read_bed(repeat_bed_path)

    gene_overlap_rows: List[Dict[str, Any]] = []

    repeat_class_counts: Dict[str, Dict[str, Any]] = defaultdict(
        lambda: {
            "repeat_indices": set(),
            "genes": set(),
            "total_overlap_bp": 0,
        }
    )

    genes_overlapping_repeats = 0
    total_gene_bp = sum(gene.length for gene in genes)
    total_gene_repeat_overlap_bp = 0

    for gene in genes:
        overlapping_repeat_segments: List[tuple[int, int]] = []
        overlapping_repeat_classes = set()
        overlapping_repeat_count = 0

        for repeat_index, repeat in enumerate(repeats):
            bp = overlap_bp(gene, repeat)

            if bp <= 0:
                continue

            overlapping_repeat_count += 1
            overlapping_repeat_classes.add(repeat.name)

            segment_start = max(gene.start, repeat.start)
            segment_end = min(gene.end, repeat.end)
            overlapping_repeat_segments.append((segment_start, segment_end))

            repeat_class_counts[repeat.name]["repeat_indices"].add(repeat_index)
            repeat_class_counts[repeat.name]["genes"].add(gene.name)
            repeat_class_counts[repeat.name]["total_overlap_bp"] += bp

        merged_segments = merge_intervals(overlapping_repeat_segments)
        repeat_overlap_bp = sum(end - start for start, end in merged_segments)

        if repeat_overlap_bp > 0:
            genes_overlapping_repeats += 1

        total_gene_repeat_overlap_bp += repeat_overlap_bp

        repeat_fraction = repeat_overlap_bp / gene.length if gene.length > 0 else 0

        gene_overlap_rows.append(
            {
                "gene_id": gene.name,
                "chrom": gene.chrom,
                "start": gene.start,
                "end": gene.end,
                "strand": gene.strand,
                "gene_length": gene.length,
                "repeat_overlap_bp": repeat_overlap_bp,
                "repeat_fraction": format_float(repeat_fraction),
                "repeat_count": overlapping_repeat_count,
                "repeat_classes": ",".join(sorted(overlapping_repeat_classes))
                if overlapping_repeat_classes
                else ".",
            }
        )

    percent_genes_overlapping_repeats = (
        genes_overlapping_repeats / len(genes) * 100 if genes else 0
    )

    summary_rows = [
        {"metric": "total_genes", "value": len(genes)},
        {"metric": "total_repeats", "value": len(repeats)},
        {"metric": "genes_overlapping_repeats", "value": genes_overlapping_repeats},
        {
            "metric": "percent_genes_overlapping_repeats",
            "value": format_float(percent_genes_overlapping_repeats),
        },
        {"metric": "total_gene_bp", "value": total_gene_bp},
        {
            "metric": "total_gene_repeat_overlap_bp",
            "value": total_gene_repeat_overlap_bp,
        },
    ]

    repeat_class_summary_rows: List[Dict[str, Any]] = []

    for repeat_class in sorted(repeat_class_counts):
        info = repeat_class_counts[repeat_class]

        repeat_class_summary_rows.append(
            {
                "repeat_class": repeat_class,
                "repeats_overlapping_genes": len(info["repeat_indices"]),
                "genes_overlapped": len(info["genes"]),
                "total_overlap_bp": info["total_overlap_bp"],
            }
        )

    return gene_overlap_rows, summary_rows, repeat_class_summary_rows


def run_repeat_overlap_analysis(
    gene_path: str | Path,
    repeat_bed_path: str | Path,
    outdir: str | Path,
) -> Dict[str, Path]:
    """
    Run gene-repeat overlap analysis and write output files.
    """
    outdir = ensure_outdir(outdir)

    gene_overlap_rows, summary_rows, repeat_class_summary_rows = analyze_gene_repeat_overlap(
        gene_path=gene_path,
        repeat_bed_path=repeat_bed_path,
    )

    gene_overlap_tsv = outdir / "genes_with_repeat_overlap.tsv"
    summary_tsv = outdir / "gene_repeat_overlap_summary.tsv"
    repeat_class_summary_tsv = outdir / "repeat_class_summary.tsv"

    write_tsv(
        gene_overlap_rows,
        gene_overlap_tsv,
        fieldnames=[
            "gene_id",
            "chrom",
            "start",
            "end",
            "strand",
            "gene_length",
            "repeat_overlap_bp",
            "repeat_fraction",
            "repeat_count",
            "repeat_classes",
        ],
    )

    write_tsv(
        summary_rows,
        summary_tsv,
        fieldnames=["metric", "value"],
    )

    write_tsv(
        repeat_class_summary_rows,
        repeat_class_summary_tsv,
        fieldnames=[
            "repeat_class",
            "repeats_overlapping_genes",
            "genes_overlapped",
            "total_overlap_bp",
        ],
    )

    return {
        "genes_with_repeat_overlap_tsv": gene_overlap_tsv,
        "gene_repeat_overlap_summary_tsv": summary_tsv,
        "repeat_class_summary_tsv": repeat_class_summary_tsv,
    }
