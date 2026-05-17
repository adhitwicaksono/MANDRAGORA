from pathlib import Path

from mandragora.intron import infer_introns_from_annotation, run_intron_analysis


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"
EXPECTED = EXAMPLES / "expected_outputs"


def normalize_text(path: Path) -> str:
    return path.read_text().strip()


def test_infer_introns_from_toy_annotation():
    annotation = EXAMPLES / "toy_annotation.gff3"

    introns, gene_summary, stats = infer_introns_from_annotation(annotation)

    assert len(introns) == 3

    intron_names = [intron.name for intron in introns]
    assert intron_names == [
        "gene1.intron1",
        "gene1.intron2",
        "gene2.intron1",
    ]

    intron_lengths = [intron.length for intron in introns]
    assert intron_lengths == [100, 100, 150]

    summary_by_gene = {row["gene_id"]: row for row in gene_summary}

    assert summary_by_gene["gene1"]["intron_count"] == 2
    assert summary_by_gene["gene2"]["intron_count"] == 1
    assert summary_by_gene["gene3"]["intron_count"] == 0

    stats_by_metric = {row["metric"]: row["value"] for row in stats}

    assert stats_by_metric["total_transcripts"] == 3
    assert stats_by_metric["total_introns"] == 3
    assert stats_by_metric["total_intron_bp"] == 350


def test_run_intron_analysis_outputs_match_expected(tmp_path):
    annotation = EXAMPLES / "toy_annotation.gff3"

    outputs = run_intron_analysis(annotation, tmp_path)

    assert normalize_text(outputs["introns_bed"]) == normalize_text(
        EXPECTED / "expected_introns.bed"
    )

    assert normalize_text(outputs["gene_intron_summary_tsv"]) == normalize_text(
        EXPECTED / "expected_gene_intron_summary.tsv"
    )

    assert normalize_text(outputs["intron_stats_tsv"]) == normalize_text(
        EXPECTED / "expected_intron_stats.tsv"
    )
