from pathlib import Path
from textwrap import wrap

base = Path("examples")
expected = base / "expected_outputs"
expected.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# 1. Toy genome
# ------------------------------------------------------------
# Mostly simple sequence, but we insert canonical-looking splice motifs
# for future splice motif testing.
scaffold_1 = list("A" * 1000)
scaffold_2 = list("C" * 800)

# Plus-strand gene1 intron 1: GFF 151-250, BED 150-250
scaffold_1[150:152] = list("GT")
scaffold_1[248:250] = list("AG")

# Plus-strand gene1 intron 2: GFF 301-400, BED 300-400
scaffold_1[300:302] = list("GT")
scaffold_1[398:400] = list("AG")

# Minus-strand gene2 intron: GFF 651-800, BED 650-800
# Designed so reverse-complement intron looks like GT...AG
scaffold_1[650:652] = list("CT")
scaffold_1[798:800] = list("AC")

def write_fasta(path, records):
    with open(path, "w") as f:
        for name, seq in records.items():
            f.write(f">{name}\n")
            for line in wrap(seq, 80):
                f.write(line + "\n")

write_fasta(
    base / "toy_genome.fa",
    {
        "scaffold_1": "".join(scaffold_1),
        "scaffold_2": "".join(scaffold_2),
    },
)

# ------------------------------------------------------------
# 2. Toy annotation GFF3
# ------------------------------------------------------------
gff3 = """##gff-version 3
##sequence-region scaffold_1 1 1000
##sequence-region scaffold_2 1 800
scaffold_1\ttoy\tgene\t101\t500\t.\t+\t.\tID=gene1;Name=Toy_gene_1
scaffold_1\ttoy\tmRNA\t101\t500\t.\t+\t.\tID=gene1.t1;Parent=gene1
scaffold_1\ttoy\texon\t101\t150\t.\t+\t.\tID=gene1.exon1;Parent=gene1.t1
scaffold_1\ttoy\texon\t251\t300\t.\t+\t.\tID=gene1.exon2;Parent=gene1.t1
scaffold_1\ttoy\texon\t401\t500\t.\t+\t.\tID=gene1.exon3;Parent=gene1.t1
scaffold_1\ttoy\tCDS\t101\t150\t.\t+\t0\tID=gene1.cds1;Parent=gene1.t1
scaffold_1\ttoy\tCDS\t251\t300\t.\t+\t0\tID=gene1.cds2;Parent=gene1.t1
scaffold_1\ttoy\tCDS\t401\t500\t.\t+\t0\tID=gene1.cds3;Parent=gene1.t1
scaffold_1\ttoy\tgene\t601\t900\t.\t-\t.\tID=gene2;Name=Toy_gene_2
scaffold_1\ttoy\tmRNA\t601\t900\t.\t-\t.\tID=gene2.t1;Parent=gene2
scaffold_1\ttoy\texon\t601\t650\t.\t-\t.\tID=gene2.exon1;Parent=gene2.t1
scaffold_1\ttoy\texon\t801\t900\t.\t-\t.\tID=gene2.exon2;Parent=gene2.t1
scaffold_1\ttoy\tCDS\t601\t650\t.\t-\t0\tID=gene2.cds1;Parent=gene2.t1
scaffold_1\ttoy\tCDS\t801\t900\t.\t-\t0\tID=gene2.cds2;Parent=gene2.t1
scaffold_2\ttoy\tgene\t101\t250\t.\t+\t.\tID=gene3;Name=Toy_gene_3_single_exon
scaffold_2\ttoy\tmRNA\t101\t250\t.\t+\t.\tID=gene3.t1;Parent=gene3
scaffold_2\ttoy\texon\t101\t250\t.\t+\t.\tID=gene3.exon1;Parent=gene3.t1
scaffold_2\ttoy\tCDS\t101\t250\t.\t+\t0\tID=gene3.cds1;Parent=gene3.t1
"""
(base / "toy_annotation.gff3").write_text(gff3)

# ------------------------------------------------------------
# 3. Toy repeat BED
# BED format: chrom, start, end, repeat_class, score, strand
# ------------------------------------------------------------
repeats_bed = """scaffold_1\t180\t260\tLTR/Gypsy\t0\t+
scaffold_1\t430\t470\tDNA/TIR\t0\t+
scaffold_1\t720\t760\tLINE/L1\t0\t-
scaffold_2\t700\t760\tSimple_repeat\t0\t+
scaffold_2\t120\t160\tLTR/Copia\t0\t+
"""
(base / "toy_repeats.bed").write_text(repeats_bed)

# ------------------------------------------------------------
# 4. Expected intron output
# ------------------------------------------------------------
expected_introns = """scaffold_1\t150\t250\tgene1.intron1\t0\t+
scaffold_1\t300\t400\tgene1.intron2\t0\t+
scaffold_1\t650\t800\tgene2.intron1\t0\t-
"""
(expected / "expected_introns.bed").write_text(expected_introns)

expected_gene_intron_summary = """gene_id\ttranscript_id\tchrom\tstrand\texon_count\tintron_count\ttotal_intron_bp\tmax_intron_length
gene1\tgene1.t1\tscaffold_1\t+\t3\t2\t200\t100
gene2\tgene2.t1\tscaffold_1\t-\t2\t1\t150\t150
gene3\tgene3.t1\tscaffold_2\t+\t1\t0\t0\t0
"""
(expected / "expected_gene_intron_summary.tsv").write_text(expected_gene_intron_summary)

expected_intron_stats = """metric\tvalue
total_transcripts\t3
total_introns\t3
total_intron_bp\t350
min_intron_length\t100
max_intron_length\t150
mean_intron_length\t116.666667
median_intron_length\t100
"""
(expected / "expected_intron_stats.tsv").write_text(expected_intron_stats)

# ------------------------------------------------------------
# 5. Expected gene-repeat overlap output
# ------------------------------------------------------------
expected_gene_repeat_overlap = """gene_id\tchrom\tstart\tend\tstrand\tgene_length\trepeat_overlap_bp\trepeat_fraction\trepeat_count\trepeat_classes
gene1\tscaffold_1\t100\t500\t+\t400\t120\t0.300000\t2\tDNA/TIR,LTR/Gypsy
gene2\tscaffold_1\t600\t900\t-\t300\t40\t0.133333\t1\tLINE/L1
gene3\tscaffold_2\t100\t250\t+\t150\t40\t0.266667\t1\tLTR/Copia
"""
(expected / "expected_gene_repeat_overlap.tsv").write_text(expected_gene_repeat_overlap)

expected_repeat_overlap_summary = """metric\tvalue
total_genes\t3
total_repeats\t5
genes_overlapping_repeats\t3
percent_genes_overlapping_repeats\t100.000000
total_gene_bp\t850
total_gene_repeat_overlap_bp\t200
"""
(expected / "expected_repeat_overlap_summary.tsv").write_text(expected_repeat_overlap_summary)

expected_repeat_class_summary = """repeat_class\trepeats_overlapping_genes\tgenes_overlapped\ttotal_overlap_bp
DNA/TIR\t1\t1\t40
LINE/L1\t1\t1\t40
LTR/Copia\t1\t1\t40
LTR/Gypsy\t1\t1\t80
"""
(expected / "expected_repeat_class_summary.tsv").write_text(expected_repeat_class_summary)

print("Toy example files created under examples/")
