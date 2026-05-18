from pathlib import Path

from mandragora.omen import analyze_gene_omen, classify_omen_level, run_omen_analysis


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"


def test_classify_omen_level():
    assert classify_omen_level(0) == "NONE"
    assert classify_omen_level(1) == "LOW"
    assert classify_omen_level(2) == "MODERATE"
    assert classify_omen_level(3) == "MODERATE"
    assert classify_omen_level(4) == "HIGH"
    assert classify_omen_level(5) == "HIGH"
    assert classify_omen_level(6) == "SEVERE"


def test_analyze_gene_omen_from_toy_files():
    annotation = EXAMPLES / "toy_annotation.gff3"
    repeats = EXAMPLES / "toy_repeats.bed"

    omen_rows, summary_rows = analyze_gene_omen(
        annotation_path=annotation,
        repeat_bed_path=repeats,
    )

    assert len(omen_rows) == 3

    by_gene = {row["gene_id"]: row for row in omen_rows}

    assert by_gene["gene1"]["gene_length"] == 400
    assert by_gene["gene1"]["intron_count"] == 2
    assert by_gene["gene1"]["repeat_fraction"] == "0.300000"
    assert by_gene["gene1"]["omen_score"] == 2
    assert by_gene["gene1"]["omen_level"] == "MODERATE"
    assert "repeat_fraction_ge_0.25" in by_gene["gene1"]["flags"]

    assert by_gene["gene2"]["gene_length"] == 300
    assert by_gene["gene2"]["intron_count"] == 1
    assert by_gene["gene2"]["repeat_fraction"] == "0.133333"
    assert by_gene["gene2"]["omen_score"] == 0
    assert by_gene["gene2"]["omen_level"] == "NONE"
    assert by_gene["gene2"]["flags"] == "."

    assert by_gene["gene3"]["gene_length"] == 150
    assert by_gene["gene3"]["exon_count"] == 1
    assert by_gene["gene3"]["repeat_fraction"] == "0.266667"
    assert by_gene["gene3"]["omen_score"] == 4
    assert by_gene["gene3"]["omen_level"] == "HIGH"
    assert "short_gene" in by_gene["gene3"]["flags"]
    assert "repeat_fraction_ge_0.25" in by_gene["gene3"]["flags"]
    assert "single_exon_repeat_overlap" in by_gene["gene3"]["flags"]

    summary = {row["metric"]: row["value"] for row in summary_rows}

    assert summary["total_genes"] == 3
    assert summary["none_omen_genes"] == 1
    assert summary["moderate_omen_genes"] == 1
    assert summary["high_omen_genes"] == 1


def test_run_omen_analysis_writes_outputs(tmp_path):
    annotation = EXAMPLES / "toy_annotation.gff3"
    repeats = EXAMPLES / "toy_repeats.bed"

    outputs = run_omen_analysis(
        annotation_path=annotation,
        repeat_bed_path=repeats,
        outdir=tmp_path,
    )

    assert outputs["gene_omen_scores_tsv"].exists()
    assert outputs["gene_omen_summary_tsv"].exists()

    scores_text = outputs["gene_omen_scores_tsv"].read_text()
    summary_text = outputs["gene_omen_summary_tsv"].read_text()

    assert "gene_id" in scores_text
    assert "gene1" in scores_text
    assert "gene2" in scores_text
    assert "gene3" in scores_text
    assert "omen_score" in scores_text
    assert "gene_omen" not in scores_text

    assert "total_genes" in summary_text
    assert "high_omen_genes" in summary_text
