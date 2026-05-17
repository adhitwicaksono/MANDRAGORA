from pathlib import Path

from mandragora.repeat_overlap import (
    analyze_gene_repeat_overlap,
    run_repeat_overlap_analysis,
)


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"
EXPECTED = EXAMPLES / "expected_outputs"


def normalize_text(path: Path) -> str:
    return path.read_text().strip()


def test_analyze_gene_repeat_overlap_from_toy_files():
    annotation = EXAMPLES / "toy_annotation.gff3"
    repeats = EXAMPLES / "toy_repeats.bed"

    gene_rows, summary_rows, class_rows = analyze_gene_repeat_overlap(
        gene_path=annotation,
        repeat_bed_path=repeats,
    )

    assert len(gene_rows) == 3

    by_gene = {row["gene_id"]: row for row in gene_rows}

    assert by_gene["gene1"]["repeat_overlap_bp"] == 120
    assert by_gene["gene1"]["repeat_fraction"] == "0.300000"
    assert by_gene["gene1"]["repeat_count"] == 2
    assert by_gene["gene1"]["repeat_classes"] == "DNA/TIR,LTR/Gypsy"

    assert by_gene["gene2"]["repeat_overlap_bp"] == 40
    assert by_gene["gene2"]["repeat_fraction"] == "0.133333"
    assert by_gene["gene2"]["repeat_classes"] == "LINE/L1"

    assert by_gene["gene3"]["repeat_overlap_bp"] == 40
    assert by_gene["gene3"]["repeat_fraction"] == "0.266667"
    assert by_gene["gene3"]["repeat_classes"] == "LTR/Copia"

    summary = {row["metric"]: row["value"] for row in summary_rows}

    assert summary["total_genes"] == 3
    assert summary["total_repeats"] == 5
    assert summary["genes_overlapping_repeats"] == 3
    assert summary["percent_genes_overlapping_repeats"] == "100.000000"
    assert summary["total_gene_bp"] == 850
    assert summary["total_gene_repeat_overlap_bp"] == 200

    by_class = {row["repeat_class"]: row for row in class_rows}

    assert by_class["LTR/Gypsy"]["total_overlap_bp"] == 80
    assert by_class["DNA/TIR"]["total_overlap_bp"] == 40
    assert by_class["LINE/L1"]["total_overlap_bp"] == 40
    assert by_class["LTR/Copia"]["total_overlap_bp"] == 40


def test_run_repeat_overlap_analysis_outputs_match_expected(tmp_path):
    annotation = EXAMPLES / "toy_annotation.gff3"
    repeats = EXAMPLES / "toy_repeats.bed"

    outputs = run_repeat_overlap_analysis(
        gene_path=annotation,
        repeat_bed_path=repeats,
        outdir=tmp_path,
    )

    assert normalize_text(outputs["genes_with_repeat_overlap_tsv"]) == normalize_text(
        EXPECTED / "expected_gene_repeat_overlap.tsv"
    )

    assert normalize_text(outputs["gene_repeat_overlap_summary_tsv"]) == normalize_text(
        EXPECTED / "expected_repeat_overlap_summary.tsv"
    )

    assert normalize_text(outputs["repeat_class_summary_tsv"]) == normalize_text(
        EXPECTED / "expected_repeat_class_summary.tsv"
    )
