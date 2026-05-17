from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

from mandragora.utils import ensure_outdir, read_bed, read_gff3_rows, write_tsv


def detect_file_type(path: str | Path) -> str:
    """
    Detect a simple file type from suffix.

    This is intentionally conservative for v0.1.
    """
    path = Path(path)
    suffixes = [suffix.lower() for suffix in path.suffixes]

    if ".gff3" in suffixes or ".gff" in suffixes:
        return "gff3"
    if ".gtf" in suffixes:
        return "gtf"
    if ".bed" in suffixes:
        return "bed"
    if ".fa" in suffixes or ".fasta" in suffixes or ".fna" in suffixes:
        return "fasta"

    return "unknown"


def count_fasta_records(fasta_path: str | Path) -> int:
    """
    Count FASTA records without requiring Biopython.

    This keeps prepare lightweight.
    """
    fasta_path = Path(fasta_path)
    count = 0

    with fasta_path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1

    return count


def run_prepare(
    annotation_path: Optional[str | Path] = None,
    repeats_path: Optional[str | Path] = None,
    genome_path: Optional[str | Path] = None,
    outdir: str | Path = "results/prepare",
) -> Dict[str, Path]:
    """
    Inspect input files and write a simple preparation report.

    Future versions will standardize GFF3/GTF/BED/RepeatMasker outputs
    into clean BED files for downstream MANDRAGORA modules.
    """
    outdir = ensure_outdir(outdir)

    rows = []

    if annotation_path is not None:
        annotation_path = Path(annotation_path)
        file_type = detect_file_type(annotation_path)

        if file_type in {"gff3", "gtf"}:
            gff_rows = read_gff3_rows(annotation_path)
            feature_types = {}

            for row in gff_rows:
                feature_type = row["type"]
                feature_types[feature_type] = feature_types.get(feature_type, 0) + 1

            rows.append(
                {
                    "input_name": "annotation",
                    "path": str(annotation_path),
                    "detected_type": file_type,
                    "records": len(gff_rows),
                    "notes": ",".join(
                        f"{key}:{value}" for key, value in sorted(feature_types.items())
                    ),
                }
            )

        else:
            rows.append(
                {
                    "input_name": "annotation",
                    "path": str(annotation_path),
                    "detected_type": file_type,
                    "records": "NA",
                    "notes": "Annotation inspection currently supports GFF3/GTF-like files.",
                }
            )

    if repeats_path is not None:
        repeats_path = Path(repeats_path)
        file_type = detect_file_type(repeats_path)

        if file_type == "bed":
            repeat_rows = read_bed(repeats_path)
            rows.append(
                {
                    "input_name": "repeats",
                    "path": str(repeats_path),
                    "detected_type": file_type,
                    "records": len(repeat_rows),
                    "notes": "BED repeat intervals loaded.",
                }
            )
        else:
            rows.append(
                {
                    "input_name": "repeats",
                    "path": str(repeats_path),
                    "detected_type": file_type,
                    "records": "NA",
                    "notes": "Repeat inspection currently supports BED. RepeatMasker .out/.gff support planned.",
                }
            )

    if genome_path is not None:
        genome_path = Path(genome_path)
        file_type = detect_file_type(genome_path)

        if file_type == "fasta":
            fasta_records = count_fasta_records(genome_path)
            rows.append(
                {
                    "input_name": "genome",
                    "path": str(genome_path),
                    "detected_type": file_type,
                    "records": fasta_records,
                    "notes": "FASTA records counted.",
                }
            )
        else:
            rows.append(
                {
                    "input_name": "genome",
                    "path": str(genome_path),
                    "detected_type": file_type,
                    "records": "NA",
                    "notes": "Genome inspection currently supports FASTA.",
                }
            )

    report_path = outdir / "prepare_report.tsv"

    write_tsv(
        rows,
        report_path,
        fieldnames=[
            "input_name",
            "path",
            "detected_type",
            "records",
            "notes",
        ],
    )

    return {
        "prepare_report_tsv": report_path,
    }
