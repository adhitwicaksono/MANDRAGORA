from pathlib import Path

from mandragora.promoter import (
    analyze_promoter_repeat_overlap,
    infer_promoters_from_annotation,
    read_fasta_lengths,
    run_promoter_analysis,
)


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"


def test_read_fasta_lengths():
    genome = EXAMPLES / "toy_genome.fa"

    lengths = read_fasta_lengths(genome)

    assert lengths["scaffold_1"] == 1000
    assert lengths["scaffold_2"] == 800


def test_infer_promoters_from_toy_annotation():
    annotation = EXAMPLES / "toy_annotation.gff3"
    genome = EXAMPLES / "toy_genome.fa"

    promoters = infer_promoters_from_annotation(
        annotation_path=annotation,
        upstream=100,
        genome_path=genome,
    )

    assert len(promoters) == 3

    by_name = {promoter.name: promoter for promoter in promoters}

    assert by_name["gene1.promoter"].chrom == "scaffold_1"
    assert by_name["gene1.promoter"].start == 0
    assert by_name["gene1.promoter"].end == 100
    assert by_name["gene1.promoter"].strand == "+"

    assert by_name["gene2.promoter"].chrom == "scaffold_1"
    assert by_name["gene2.promoter"].start == 900
    assert by_name["gene2.promoter"].end == 1000
    assert by_name["gene2.promoter"].strand == "-"

    assert by_name["gene3.promoter"].chrom == "scaffold_2"
    assert by_name["gene3.promoter"].start == 0
    assert by_name["gene3.promoter"].end == 100
    assert by_name["gene3.promoter"].strand == "+"


def test_promoter_repeat_overlap_from_toy_files():
    annotation = EXAMPLES / "toy_annotation.gff3"
    genome = EXAMPLES / "toy_genome.fa"
    repeats = EXAMPLES / "toy_promoter_repeats.bed"

    promoters = infer_promoters_from_annotation(
        annotation_path=annotation,
        upstream=100,
        genome_path=genome,
    )

    rows = analyze_promoter_repeat_overlap(
        promoters=promoters,
        repeat_bed_path=repeats,
    )

    by_promoter = {row["promoter_id"]: row for row in rows}

    assert by_promoter["gene1.promoter"]["repeat_overlap_bp"] == 40
    assert by_promoter["gene1.promoter"]["repeat_fraction"] == "0.400000"
    assert by_promoter["gene1.promoter"]["repeat_classes"] == "Promoter_repeat_A"

    assert by_promoter["gene2.promoter"]["repeat_overlap_bp"] == 30
    assert by_promoter["gene2.promoter"]["repeat_fraction"] == "0.300000"
    assert by_promoter["gene2.promoter"]["repeat_classes"] == "Promoter_repeat_B"

    assert by_promoter["gene3.promoter"]["repeat_overlap_bp"] == 10
    assert by_promoter["gene3.promoter"]["repeat_fraction"] == "0.100000"
    assert by_promoter["gene3.promoter"]["repeat_classes"] == "Promoter_repeat_C"


def test_run_promoter_analysis_writes_outputs(tmp_path):
    annotation = EXAMPLES / "toy_annotation.gff3"
    genome = EXAMPLES / "toy_genome.fa"
    repeats = EXAMPLES / "toy_promoter_repeats.bed"

    outputs = run_promoter_analysis(
        annotation_path=annotation,
        genome_path=genome,
        repeat_bed_path=repeats,
        upstream=100,
        outdir=tmp_path,
    )

    assert outputs["promoters_bed"].exists()
    assert outputs["promoters_fa"].exists()
    assert outputs["promoter_repeat_overlap_tsv"].exists()
    assert outputs["promoter_summary_tsv"].exists()

    promoters_bed = outputs["promoters_bed"].read_text()
    promoter_fasta = outputs["promoters_fa"].read_text()
    summary = outputs["promoter_summary_tsv"].read_text()

    assert "gene1.promoter" in promoters_bed
    assert "gene2.promoter" in promoters_bed
    assert "gene3.promoter" in promoters_bed

    assert ">gene1.promoter" in promoter_fasta
    assert ">gene2.promoter" in promoter_fasta
    assert ">gene3.promoter" in promoter_fasta

    assert "total_promoters" in summary
    assert "promoters_overlapping_repeats" in summary
