from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional

from mandragora.repeat_overlap import load_gene_intervals
from mandragora.utils import (
    Interval,
    ensure_outdir,
    format_float,
    merge_intervals,
    overlap_bp,
    read_bed,
    write_bed,
    write_tsv,
)


def read_fasta_lengths(genome_path: str | Path) -> Dict[str, int]:
    """
    Read FASTA record lengths without loading full sequences into memory.
    """
    genome_path = Path(genome_path)
    lengths: Dict[str, int] = {}
    current_name: Optional[str] = None
    current_length = 0

    with genome_path.open() as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_name is not None:
                    lengths[current_name] = current_length

                current_name = line[1:].split()[0]
                current_length = 0
            else:
                current_length += len(line)

    if current_name is not None:
        lengths[current_name] = current_length

    return lengths


def reverse_complement(sequence: str) -> str:
    """
    Reverse-complement a DNA sequence.
    """
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return sequence.translate(table)[::-1]


def extract_fasta_sequence(
    genome_path: str | Path,
    chrom: str,
    start: int,
    end: int,
    strand: str = "+",
) -> str:
    """
    Extract sequence from FASTA using pyfaidx.

    Coordinates are 0-based, half-open.
    If strand is '-', the reverse-complement is returned.
    """
    try:
        from pyfaidx import Fasta
    except ImportError as exc:
        raise ImportError(
            "pyfaidx is required for FASTA extraction. "
            "Install it with: pip install pyfaidx"
        ) from exc

    fasta = Fasta(str(genome_path), rebuild=False)

    sequence = str(fasta[chrom][start:end])

    if strand == "-":
        sequence = reverse_complement(sequence)

    return sequence


def infer_promoters_from_annotation(
    annotation_path: str | Path,
    upstream: int = 2000,
    genome_path: Optional[str | Path] = None,
) -> List[Interval]:
    """
    Infer promoter/upstream intervals from gene coordinates.

    Coordinates are BED-style:
    0-based, half-open.

    Plus-strand gene:
      promoter = gene.start - upstream to gene.start

    Minus-strand gene:
      promoter = gene.end to gene.end + upstream

    If genome_path is provided, promoter coordinates are clamped to scaffold bounds.
    """
    genes = load_gene_intervals(annotation_path)

    scaffold_lengths: Dict[str, int] = {}
    if genome_path is not None:
        scaffold_lengths = read_fasta_lengths(genome_path)

    promoters: List[Interval] = []

    for gene in genes:
        scaffold_length = scaffold_lengths.get(gene.chrom)

        if gene.strand == "-":
            promoter_start = gene.end
            promoter_end = gene.end + upstream

            if scaffold_length is not None:
                promoter_end = min(promoter_end, scaffold_length)
        else:
            promoter_start = max(0, gene.start - upstream)
            promoter_end = gene.start

            if scaffold_length is not None:
                promoter_start = max(0, promoter_start)
                promoter_end = min(promoter_end, scaffold_length)

        if promoter_end < promoter_start:
            promoter_end = promoter_start

        promoters.append(
            Interval(
                chrom=gene.chrom,
                start=promoter_start,
                end=promoter_end,
                name=f"{gene.name}.promoter",
                score="0",
                strand=gene.strand,
            )
        )

    return promoters


def write_promoter_fasta(
    promoters: List[Interval],
    genome_path: str | Path,
    output_path: str | Path,
) -> None:
    """
    Write promoter sequences to FASTA.
    """
    output_path = Path(output_path)

    with output_path.open("w") as handle:
        for promoter in promoters:
            sequence = extract_fasta_sequence(
                genome_path=genome_path,
                chrom=promoter.chrom,
                start=promoter.start,
                end=promoter.end,
                strand=promoter.strand,
            )

            header = (
                f">{promoter.name} "
                f"{promoter.chrom}:{promoter.start}-{promoter.end}({promoter.strand})"
            )

            handle.write(header + "\n")

            for i in range(0, len(sequence), 80):
                handle.write(sequence[i : i + 80] + "\n")


def analyze_promoter_repeat_overlap(
    promoters: List[Interval],
    repeat_bed_path: str | Path,
) -> List[Dict[str, Any]]:
    """
    Analyze overlap between promoter intervals and repeat intervals.
    """
    repeats = read_bed(repeat_bed_path)
    rows: List[Dict[str, Any]] = []

    for promoter in promoters:
        overlapping_segments: List[tuple[int, int]] = []
        repeat_classes = set()
        repeat_count = 0

        for repeat in repeats:
            bp = overlap_bp(promoter, repeat)

            if bp <= 0:
                continue

            repeat_count += 1
            repeat_classes.add(repeat.name)

            overlapping_segments.append(
                (
                    max(promoter.start, repeat.start),
                    min(promoter.end, repeat.end),
                )
            )

        merged_segments = merge_intervals(overlapping_segments)
        repeat_overlap_bp = sum(end - start for start, end in merged_segments)

        promoter_length = promoter.length
        repeat_fraction = (
            repeat_overlap_bp / promoter_length if promoter_length > 0 else 0.0
        )

        rows.append(
            {
                "promoter_id": promoter.name,
                "chrom": promoter.chrom,
                "start": promoter.start,
                "end": promoter.end,
                "strand": promoter.strand,
                "promoter_length": promoter_length,
                "repeat_overlap_bp": repeat_overlap_bp,
                "repeat_fraction": format_float(repeat_fraction),
                "repeat_count": repeat_count,
                "repeat_classes": ",".join(sorted(repeat_classes))
                if repeat_classes
                else ".",
            }
        )

    return rows


def summarize_promoters(
    promoters: List[Interval],
    promoter_repeat_rows: Optional[List[Dict[str, Any]]] = None,
) -> List[Dict[str, Any]]:
    """
    Create promoter summary rows.
    """
    lengths = [promoter.length for promoter in promoters]

    rows = [
        {"metric": "total_promoters", "value": len(promoters)},
        {"metric": "total_promoter_bp", "value": sum(lengths)},
        {"metric": "min_promoter_length", "value": min(lengths) if lengths else 0},
        {"metric": "max_promoter_length", "value": max(lengths) if lengths else 0},
    ]

    if promoter_repeat_rows is not None:
        promoters_with_repeats = sum(
            1 for row in promoter_repeat_rows if int(row["repeat_overlap_bp"]) > 0
        )
        total_repeat_overlap_bp = sum(
            int(row["repeat_overlap_bp"]) for row in promoter_repeat_rows
        )

        rows.extend(
            [
                {
                    "metric": "promoters_overlapping_repeats",
                    "value": promoters_with_repeats,
                },
                {
                    "metric": "percent_promoters_overlapping_repeats",
                    "value": format_float(
                        promoters_with_repeats / len(promoters) * 100
                        if promoters
                        else 0.0
                    ),
                },
                {
                    "metric": "total_promoter_repeat_overlap_bp",
                    "value": total_repeat_overlap_bp,
                },
            ]
        )

    return rows


def run_promoter_analysis(
    annotation_path: str | Path,
    outdir: str | Path,
    upstream: int = 2000,
    genome_path: Optional[str | Path] = None,
    repeat_bed_path: Optional[str | Path] = None,
) -> Dict[str, Path]:
    """
    Run promoter extraction and optional promoter-repeat analysis.
    """
    outdir = ensure_outdir(outdir)

    promoters = infer_promoters_from_annotation(
        annotation_path=annotation_path,
        upstream=upstream,
        genome_path=genome_path,
    )

    outputs: Dict[str, Path] = {}

    promoters_bed = outdir / "promoters.bed"
    write_bed(promoters, promoters_bed)
    outputs["promoters_bed"] = promoters_bed

    promoter_repeat_rows = None

    if genome_path is not None:
        promoters_fa = outdir / "promoters.fa"
        write_promoter_fasta(
            promoters=promoters,
            genome_path=genome_path,
            output_path=promoters_fa,
        )
        outputs["promoters_fa"] = promoters_fa

    if repeat_bed_path is not None:
        promoter_repeat_rows = analyze_promoter_repeat_overlap(
            promoters=promoters,
            repeat_bed_path=repeat_bed_path,
        )

        promoter_repeat_tsv = outdir / "promoter_repeat_overlap.tsv"

        write_tsv(
            promoter_repeat_rows,
            promoter_repeat_tsv,
            fieldnames=[
                "promoter_id",
                "chrom",
                "start",
                "end",
                "strand",
                "promoter_length",
                "repeat_overlap_bp",
                "repeat_fraction",
                "repeat_count",
                "repeat_classes",
            ],
        )

        outputs["promoter_repeat_overlap_tsv"] = promoter_repeat_tsv

    promoter_summary_rows = summarize_promoters(
        promoters=promoters,
        promoter_repeat_rows=promoter_repeat_rows,
    )

    promoter_summary_tsv = outdir / "promoter_summary.tsv"

    write_tsv(
        promoter_summary_rows,
        promoter_summary_tsv,
        fieldnames=["metric", "value"],
    )

    outputs["promoter_summary_tsv"] = promoter_summary_tsv

    return outputs
