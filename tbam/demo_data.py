from __future__ import annotations

import random
from pathlib import Path

import pysam

DNA = "ACGT"


def generate_demo_bam(
    out_dir: Path,
    *,
    reads: int = 250,
    prefix: str = "demo",
) -> tuple[Path, Path, Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    bam_path = out_dir / f"{prefix}.bam"
    bai_path = out_dir / f"{prefix}.bam.bai"
    fasta_path = out_dir / f"{prefix}.fa"
    fai_path = out_dir / f"{prefix}.fa.fai"

    rng = random.Random(42)
    references = [("chr1", 60_000), ("chr2", 40_000), ("chrM", 16_569)]
    reference_sequences = _build_reference_sequences(references, rng=rng)
    _write_fasta(fasta_path, reference_sequences)
    pysam.faidx(str(fasta_path))

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": length} for name, length in references],
    }

    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out_bam:
        read_id = 0
        for ref_index, (_, ref_length) in enumerate(references):
            count = _reads_for_reference(reads, ref_index, len(references))
            if count <= 0:
                continue
            hotspots = _reference_hotspots(ref_length)

            for i in range(count):
                seq_len = 75
                hotspot = hotspots[i % len(hotspots)]
                jitter = rng.randint(-120, 120)
                start = max(100, min(ref_length - seq_len - 1, hotspot + jitter))
                ref_name = references[ref_index][0]
                ref_seq = reference_sequences[ref_name][start : start + seq_len]

                read = pysam.AlignedSegment()
                read.query_name = f"read_{read_id:06d}"
                read.query_sequence, nm = _mutate_sequence(ref_seq, rng=rng, max_mutations=2)
                read.flag = 0
                read.reference_id = ref_index
                read.reference_start = start
                read.mapping_quality = rng.randint(18, 60)
                read.cigar = ((0, seq_len),)
                read.next_reference_id = -1
                read.next_reference_start = -1
                read.template_length = 0
                read.query_qualities = pysam.qualitystring_to_array("I" * seq_len)
                read.set_tag("NM", nm)
                read.set_tag("AS", seq_len - rng.randint(0, 6))
                out_bam.write(read)
                read_id += 1

    pysam.index(str(bam_path))
    return bam_path, bai_path, fasta_path, fai_path


def _reads_for_reference(total_reads: int, ref_index: int, ref_count: int) -> int:
    base = total_reads // ref_count
    remainder = total_reads % ref_count
    return base + (1 if ref_index < remainder else 0)


def _build_reference_sequences(
    references: list[tuple[str, int]],
    *,
    rng: random.Random,
) -> dict[str, str]:
    return {
        name: "".join(rng.choice(DNA) for _ in range(length))
        for name, length in references
    }


def _write_fasta(path: Path, references: dict[str, str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, seq in references.items():
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                handle.write(seq[i : i + 60] + "\n")


def _mutate_sequence(seq: str, *, rng: random.Random, max_mutations: int) -> tuple[str, int]:
    if not seq:
        return seq, 0
    out = list(seq)
    changes = rng.randint(0, max_mutations)
    for _ in range(changes):
        idx = rng.randrange(0, len(out))
        old = out[idx]
        choices = [base for base in DNA if base != old]
        out[idx] = rng.choice(choices)
    return "".join(out), changes


def _reference_hotspots(ref_length: int) -> list[int]:
    anchor = [2_000, 6_000, 10_000]
    max_start = max(500, ref_length - 500)
    return [min(pos, max_start) for pos in anchor]
