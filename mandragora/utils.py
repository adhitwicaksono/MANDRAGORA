from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional


@dataclass(frozen=True)
class Interval:
    """
    Simple BED-like interval.

    Coordinates are 0-based, half-open:
    start is included, end is excluded.
    """

    chrom: str
    start: int
    end: int
    name: str = "."
    score: str = "0"
    strand: str = "."

    @property
    def length(self) -> int:
        return max(0, self.end - self.start)


def ensure_outdir(outdir: str | Path) -> Path:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    return outdir


def parse_gff_attributes(attribute_text: str) -> Dict[str, str]:
    """
    Parse GFF3 attributes, e.g.:

    ID=gene1;Name=ABC;Parent=gene0

    Returns a dictionary.
    """
    attributes: Dict[str, str] = {}

    if not attribute_text or attribute_text == ".":
        return attributes

    for item in attribute_text.strip().split(";"):
        if not item:
            continue

        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            # loose GTF-ish fallback
            key, value = item.split(" ", 1)
            value = value.strip('"')
        else:
            continue

        attributes[key.strip()] = value.strip()

    return attributes


def read_gff3_rows(gff_path: str | Path) -> List[Dict[str, Any]]:
    """
    Read GFF3-like rows and convert coordinates to BED convention.

    GFF3:
      1-based, inclusive

    Internal/MANDRAGORA:
      0-based, half-open
    """
    rows: List[Dict[str, Any]] = []
    gff_path = Path(gff_path)

    with gff_path.open() as handle:
        for line in handle:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")

            if len(parts) != 9:
                continue

            chrom, source, feature_type, start, end, score, strand, phase, attrs = parts

            start_i = int(start)
            end_i = int(end)

            rows.append(
                {
                    "chrom": chrom,
                    "source": source,
                    "type": feature_type,
                    "start": start_i - 1,
                    "end": end_i,
                    "score": score,
                    "strand": strand if strand in {"+", "-"} else ".",
                    "phase": phase,
                    "attributes": parse_gff_attributes(attrs),
                    "raw_attributes": attrs,
                }
            )

    return rows


def extract_gene_intervals_from_gff3(gff_path: str | Path) -> List[Interval]:
    """
    Extract gene intervals from a GFF3 file.
    """
    rows = read_gff3_rows(gff_path)
    genes: List[Interval] = []

    for row in rows:
        if row["type"] != "gene":
            continue

        attrs = row["attributes"]
        gene_id = attrs.get("ID") or attrs.get("gene_id") or attrs.get("Name") or "unknown_gene"

        genes.append(
            Interval(
                chrom=row["chrom"],
                start=row["start"],
                end=row["end"],
                name=gene_id,
                score="0",
                strand=row["strand"],
            )
        )

    return genes


def extract_transcript_exons_from_gff3(gff_path: str | Path) -> List[Dict[str, Any]]:
    """
    Extract transcript-exon structures from a GFF3 file.

    Returns a list of dictionaries:

    {
      "gene_id": "...",
      "transcript_id": "...",
      "chrom": "...",
      "strand": "+",
      "exons": [Interval(...), ...]
    }
    """
    rows = read_gff3_rows(gff_path)

    genes: Dict[str, Dict[str, Any]] = {}
    transcripts: Dict[str, Dict[str, Any]] = {}
    exons_by_transcript: Dict[str, List[Interval]] = defaultdict(list)

    transcript_feature_types = {"mRNA", "transcript"}

    for row in rows:
        attrs = row["attributes"]
        feature_type = row["type"]

        if feature_type == "gene":
            gene_id = attrs.get("ID") or attrs.get("gene_id") or attrs.get("Name")

            if gene_id:
                genes[gene_id] = {
                    "gene_id": gene_id,
                    "chrom": row["chrom"],
                    "start": row["start"],
                    "end": row["end"],
                    "strand": row["strand"],
                }

        elif feature_type in transcript_feature_types:
            transcript_id = attrs.get("ID") or attrs.get("transcript_id")
            parent_gene = attrs.get("Parent") or attrs.get("gene_id")

            if transcript_id:
                transcripts[transcript_id] = {
                    "transcript_id": transcript_id,
                    "gene_id": parent_gene or transcript_id,
                    "chrom": row["chrom"],
                    "start": row["start"],
                    "end": row["end"],
                    "strand": row["strand"],
                }

        elif feature_type == "exon":
            parent_text = attrs.get("Parent") or attrs.get("transcript_id") or attrs.get("gene_id")

            if not parent_text:
                continue

            parents = [p.strip() for p in parent_text.split(",") if p.strip()]

            for parent in parents:
                exon_name = attrs.get("ID") or f"{parent}.exon"

                exons_by_transcript[parent].append(
                    Interval(
                        chrom=row["chrom"],
                        start=row["start"],
                        end=row["end"],
                        name=exon_name,
                        score="0",
                        strand=row["strand"],
                    )
                )

    transcript_models: List[Dict[str, Any]] = []

    for transcript_id, exons in sorted(exons_by_transcript.items()):
        transcript_info = transcripts.get(transcript_id)

        if transcript_info:
            gene_id = transcript_info["gene_id"]
            chrom = transcript_info["chrom"]
            strand = transcript_info["strand"]
        else:
            # Fallback: exon Parent may directly point to a gene
            gene_info = genes.get(transcript_id, {})
            gene_id = transcript_id
            chrom = gene_info.get("chrom", exons[0].chrom)
            strand = gene_info.get("strand", exons[0].strand)

        exons_sorted = sorted(exons, key=lambda x: (x.chrom, x.start, x.end))

        transcript_models.append(
            {
                "gene_id": gene_id,
                "transcript_id": transcript_id,
                "chrom": chrom,
                "strand": strand,
                "exons": exons_sorted,
            }
        )

    return transcript_models


def read_bed(bed_path: str | Path) -> List[Interval]:
    """
    Read a BED-like file.

    Requires at least 3 columns:
      chrom, start, end

    Optional:
      name, score, strand
    """
    bed_path = Path(bed_path)
    intervals: List[Interval] = []

    with bed_path.open() as handle:
        for line in handle:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")

            if len(parts) < 3:
                continue

            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3] if len(parts) >= 4 else "."
            score = parts[4] if len(parts) >= 5 else "0"
            strand = parts[5] if len(parts) >= 6 and parts[5] in {"+", "-"} else "."

            intervals.append(
                Interval(
                    chrom=chrom,
                    start=start,
                    end=end,
                    name=name,
                    score=score,
                    strand=strand,
                )
            )

    return intervals


def write_bed(intervals: Iterable[Interval], output_path: str | Path) -> None:
    output_path = Path(output_path)

    with output_path.open("w") as handle:
        for interval in intervals:
            handle.write(
                "\t".join(
                    [
                        interval.chrom,
                        str(interval.start),
                        str(interval.end),
                        interval.name,
                        str(interval.score),
                        interval.strand,
                    ]
                )
                + "\n"
            )


def write_tsv(rows: Iterable[Dict[str, Any]], output_path: str | Path, fieldnames: List[str]) -> None:
    output_path = Path(output_path)

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for row in rows:
            writer.writerow(row)


def overlap_bp(a: Interval, b: Interval) -> int:
    """
    Return overlap length between two BED-like intervals.
    """
    if a.chrom != b.chrom:
        return 0

    return max(0, min(a.end, b.end) - max(a.start, b.start))


def merge_intervals(intervals: List[tuple[int, int]]) -> List[tuple[int, int]]:
    """
    Merge overlapping coordinate intervals.
    """
    if not intervals:
        return []

    intervals = sorted(intervals)
    merged = [intervals[0]]

    for start, end in intervals[1:]:
        last_start, last_end = merged[-1]

        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))

    return merged


def format_float(value: float) -> str:
    return f"{value:.6f}"
