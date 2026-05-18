from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any, Dict, List

from mandragora.intron import infer_introns_from_annotation
from mandragora.repeat_overlap import analyze_gene_repeat_overlap, load_gene_intervals
from mandragora.utils import ensure_outdir, format_float, write_tsv


def classify_omen_level(score: int) -> str:
    """
    Convert numeric omen score into a qualitative warning level.
    """
    if score <= 0:
        return "NONE"
    if score == 1:
        return "LOW"
    if 2 <= score <= 3:
        return "MODERATE"
    if 4 <= score <= 5:
        return "HIGH"
    return "SEVERE"


def analyze_gene_omen(
    annotation_path: str | Path,
    repeat_bed_path: str | Path,
    long_gene_threshold: int = 100_000,
    very_long_gene_threshold: int = 500_000,
    long_intron_threshold: int = 100_000,
    very_long_intron_threshold: int = 500_000,
    repeat_warning_fraction: float = 0.25,
    repeat_high_fraction: float = 0.50,
    repeat_severe_fraction: float = 0.75,
    short_gene_threshold: int = 300,
) -> tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """
    Generate gene-level MANDRAGORA omen scores.

    The omen score is not a biological truth label.
    It is a diagnostic warning score to prioritize manual inspection.
    """
    genes = load_gene_intervals(annotation_path)

    _, gene_intron_summary, _ = infer_introns_from_annotation(annotation_path)

    repeat_rows, _, _ = analyze_gene_repeat_overlap(
        gene_path=annotation_path,
        repeat_bed_path=repeat_bed_path,
    )

    intron_by_gene: Dict[str, Dict[str, Any]] = {}

    for row in gene_intron_summary:
        gene_id = row["gene_id"]

        current = intron_by_gene.get(
            gene_id,
            {
                "exon_count": 0,
                "intron_count": 0,
                "total_intron_bp": 0,
                "max_intron_length": 0,
            },
        )

        # For genes with multiple transcript models, keep the maximum values.
        current["exon_count"] = max(current["exon_count"], int(row["exon_count"]))
        current["intron_count"] = max(current["intron_count"], int(row["intron_count"]))
        current["total_intron_bp"] = max(
            current["total_intron_bp"], int(row["total_intron_bp"])
        )
        current["max_intron_length"] = max(
            current["max_intron_length"], int(row["max_intron_length"])
        )

        intron_by_gene[gene_id] = current

    repeat_by_gene = {row["gene_id"]: row for row in repeat_rows}

    omen_rows: List[Dict[str, Any]] = []

    for gene in genes:
        gene_id = gene.name
        intron_info = intron_by_gene.get(
            gene_id,
            {
                "exon_count": 0,
                "intron_count": 0,
                "total_intron_bp": 0,
                "max_intron_length": 0,
            },
        )

        repeat_info = repeat_by_gene.get(
            gene_id,
            {
                "repeat_overlap_bp": 0,
                "repeat_fraction": "0.000000",
                "repeat_count": 0,
                "repeat_classes": ".",
            },
        )

        repeat_fraction = float(repeat_info["repeat_fraction"])
        max_intron_length = int(intron_info["max_intron_length"])
        exon_count = int(intron_info["exon_count"])
        intron_count = int(intron_info["intron_count"])
        repeat_overlap_bp = int(repeat_info["repeat_overlap_bp"])
        repeat_count = int(repeat_info["repeat_count"])

        score = 0
        flags: List[str] = []

        if gene.length < short_gene_threshold:
            score += 1
            flags.append("short_gene")

        if gene.length >= very_long_gene_threshold:
            score += 2
            flags.append("very_long_gene")
        elif gene.length >= long_gene_threshold:
            score += 1
            flags.append("long_gene")

        if max_intron_length >= very_long_intron_threshold:
            score += 2
            flags.append("very_long_intron")
        elif max_intron_length >= long_intron_threshold:
            score += 1
            flags.append("long_intron")

        if repeat_fraction >= repeat_severe_fraction:
            score += 4
            flags.append("repeat_fraction_ge_0.75")
        elif repeat_fraction >= repeat_high_fraction:
            score += 3
            flags.append("repeat_fraction_ge_0.50")
        elif repeat_fraction >= repeat_warning_fraction:
            score += 2
            flags.append("repeat_fraction_ge_0.25")

        if exon_count == 1 and repeat_overlap_bp > 0:
            score += 1
            flags.append("single_exon_repeat_overlap")

        omen_rows.append(
            {
                "gene_id": gene_id,
                "chrom": gene.chrom,
                "start": gene.start,
                "end": gene.end,
                "strand": gene.strand,
                "gene_length": gene.length,
                "exon_count": exon_count,
                "intron_count": intron_count,
                "max_intron_length": max_intron_length,
                "repeat_overlap_bp": repeat_overlap_bp,
                "repeat_fraction": format_float(repeat_fraction),
                "repeat_count": repeat_count,
                "repeat_classes": repeat_info["repeat_classes"],
                "omen_score": score,
                "omen_level": classify_omen_level(score),
                "flags": ",".join(flags) if flags else ".",
            }
        )

    level_counts = Counter(row["omen_level"] for row in omen_rows)

    summary_rows = [
        {"metric": "total_genes", "value": len(omen_rows)},
        {"metric": "none_omen_genes", "value": level_counts.get("NONE", 0)},
        {"metric": "low_omen_genes", "value": level_counts.get("LOW", 0)},
        {"metric": "moderate_omen_genes", "value": level_counts.get("MODERATE", 0)},
        {"metric": "high_omen_genes", "value": level_counts.get("HIGH", 0)},
        {"metric": "severe_omen_genes", "value": level_counts.get("SEVERE", 0)},
    ]

    return omen_rows, summary_rows


def run_omen_analysis(
    annotation_path: str | Path,
    repeat_bed_path: str | Path,
    outdir: str | Path,
    long_gene_threshold: int = 100_000,
    very_long_gene_threshold: int = 500_000,
    long_intron_threshold: int = 100_000,
    very_long_intron_threshold: int = 500_000,
    repeat_warning_fraction: float = 0.25,
    repeat_high_fraction: float = 0.50,
    repeat_severe_fraction: float = 0.75,
    short_gene_threshold: int = 300,
) -> Dict[str, Path]:
    """
    Run gene omen analysis and write output files.
    """
    outdir = ensure_outdir(outdir)

    omen_rows, summary_rows = analyze_gene_omen(
        annotation_path=annotation_path,
        repeat_bed_path=repeat_bed_path,
        long_gene_threshold=long_gene_threshold,
        very_long_gene_threshold=very_long_gene_threshold,
        long_intron_threshold=long_intron_threshold,
        very_long_intron_threshold=very_long_intron_threshold,
        repeat_warning_fraction=repeat_warning_fraction,
        repeat_high_fraction=repeat_high_fraction,
        repeat_severe_fraction=repeat_severe_fraction,
        short_gene_threshold=short_gene_threshold,
    )

    omen_scores_tsv = outdir / "gene_omen_scores.tsv"
    omen_summary_tsv = outdir / "gene_omen_summary.tsv"

    write_tsv(
        omen_rows,
        omen_scores_tsv,
        fieldnames=[
            "gene_id",
            "chrom",
            "start",
            "end",
            "strand",
            "gene_length",
            "exon_count",
            "intron_count",
            "max_intron_length",
            "repeat_overlap_bp",
            "repeat_fraction",
            "repeat_count",
            "repeat_classes",
            "omen_score",
            "omen_level",
            "flags",
        ],
    )

    write_tsv(
        summary_rows,
        omen_summary_tsv,
        fieldnames=["metric", "value"],
    )

    return {
        "gene_omen_scores_tsv": omen_scores_tsv,
        "gene_omen_summary_tsv": omen_summary_tsv,
    }
