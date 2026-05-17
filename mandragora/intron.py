from __future__ import annotations

from pathlib import Path
from statistics import mean, median
from typing import Any, Dict, List

from mandragora.utils import (
    Interval,
    ensure_outdir,
    extract_transcript_exons_from_gff3,
    write_bed,
    write_tsv,
)


def infer_introns_from_annotation(annotation_path: str | Path) -> tuple[List[Interval], List[Dict[str, Any]], List[Dict[str, Any]]]:
    """
    Infer introns from exon coordinates in a GFF3 annotation.

    Coordinates are reported in BED convention:
    0-based, half-open.
    """
    transcript_models = extract_transcript_exons_from_gff3(annotation_path)

    all_introns: List[Interval] = []
    gene_intron_summary: List[Dict[str, Any]] = []

    for model in transcript_models:
        gene_id = model["gene_id"]
        transcript_id = model["transcript_id"]
        chrom = model["chrom"]
        strand = model["strand"]
        exons = model["exons"]

        introns_for_transcript: List[Interval] = []

        if len(exons) >= 2:
            exons_sorted = sorted(exons, key=lambda x: (x.chrom, x.start, x.end))

            for index in range(len(exons_sorted) - 1):
                left_exon = exons_sorted[index]
                right_exon = exons_sorted[index + 1]

                intron_start = left_exon.end
                intron_end = right_exon.start

                if intron_end <= intron_start:
                    continue

                intron = Interval(
                    chrom=chrom,
                    start=intron_start,
                    end=intron_end,
                    name=f"{gene_id}.intron{index + 1}",
                    score="0",
                    strand=strand,
                )

                introns_for_transcript.append(intron)
                all_introns.append(intron)

        intron_lengths = [intron.length for intron in introns_for_transcript]

        gene_intron_summary.append(
            {
                "gene_id": gene_id,
                "transcript_id": transcript_id,
                "chrom": chrom,
                "strand": strand,
                "exon_count": len(exons),
                "intron_count": len(introns_for_transcript),
                "total_intron_bp": sum(intron_lengths),
                "max_intron_length": max(intron_lengths) if intron_lengths else 0,
            }
        )

    intron_lengths_all = [intron.length for intron in all_introns]

    if intron_lengths_all:
        intron_stats = [
            {"metric": "total_transcripts", "value": len(transcript_models)},
            {"metric": "total_introns", "value": len(all_introns)},
            {"metric": "total_intron_bp", "value": sum(intron_lengths_all)},
            {"metric": "min_intron_length", "value": min(intron_lengths_all)},
            {"metric": "max_intron_length", "value": max(intron_lengths_all)},
            {"metric": "mean_intron_length", "value": f"{mean(intron_lengths_all):.6f}"},
            {"metric": "median_intron_length", "value": median(intron_lengths_all)},
        ]
    else:
        intron_stats = [
            {"metric": "total_transcripts", "value": len(transcript_models)},
            {"metric": "total_introns", "value": 0},
            {"metric": "total_intron_bp", "value": 0},
            {"metric": "min_intron_length", "value": 0},
            {"metric": "max_intron_length", "value": 0},
            {"metric": "mean_intron_length", "value": "0.000000"},
            {"metric": "median_intron_length", "value": 0},
        ]

    return all_introns, gene_intron_summary, intron_stats


def run_intron_analysis(annotation_path: str | Path, outdir: str | Path) -> Dict[str, Path]:
    """
    Run intron analysis and write output files.
    """
    outdir = ensure_outdir(outdir)

    introns, gene_intron_summary, intron_stats = infer_introns_from_annotation(annotation_path)

    introns_bed = outdir / "introns.bed"
    gene_summary_tsv = outdir / "gene_intron_summary.tsv"
    intron_stats_tsv = outdir / "intron_stats.tsv"

    write_bed(introns, introns_bed)

    write_tsv(
        gene_intron_summary,
        gene_summary_tsv,
        fieldnames=[
            "gene_id",
            "transcript_id",
            "chrom",
            "strand",
            "exon_count",
            "intron_count",
            "total_intron_bp",
            "max_intron_length",
        ],
    )

    write_tsv(
        intron_stats,
        intron_stats_tsv,
        fieldnames=["metric", "value"],
    )

    return {
        "introns_bed": introns_bed,
        "gene_intron_summary_tsv": gene_summary_tsv,
        "intron_stats_tsv": intron_stats_tsv,
    }
