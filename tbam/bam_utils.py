from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pysam


@dataclass(slots=True)
class ReadRow:
    index: int
    qname: str
    rname: str
    pos: int
    mapq: int
    cigar: str
    flag: int
    length: int


@dataclass(slots=True)
class BamSummary:
    path: Path
    total_reads: int
    mapped_reads: int
    unmapped_reads: int
    references: list[tuple[str, int]]


def validate_bam(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"BAM file not found: {path}")
    if path.suffix.lower() != ".bam":
        raise ValueError(f"Expected a .bam file, got: {path.name}")


def validate_reference_fasta(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Reference FASTA file not found: {path}")


def fetch_rows(
    bam_path: Path,
    *,
    max_reads: int = 1000,
    contig: str | None = None,
    min_mapq: int = 0,
) -> list[ReadRow]:
    validate_bam(bam_path)
    rows: list[ReadRow] = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        reads = _iter_reads(bam, contig=contig)
        for idx, read in enumerate(reads):
            if read.mapq < min_mapq:
                continue
            if len(rows) >= max_reads:
                break
            rows.append(
                ReadRow(
                    index=idx,
                    qname=read.query_name or "",
                    rname=read.reference_name or "*",
                    pos=(read.reference_start + 1) if read.reference_start >= 0 else 0,
                    mapq=read.mapping_quality,
                    cigar=read.cigarstring or "*",
                    flag=read.flag,
                    length=read.query_length or 0,
                )
            )
    return rows


def get_header_text(bam_path: Path) -> str:
    validate_bam(bam_path)
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        return bam.text or ""


def summarize_bam(bam_path: Path) -> BamSummary:
    validate_bam(bam_path)
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        total = 0
        mapped = 0
        unmapped = 0
        for read in bam.fetch(until_eof=True):
            total += 1
            if read.is_unmapped:
                unmapped += 1
            else:
                mapped += 1
        references = list(zip(bam.references, bam.lengths))
    return BamSummary(
        path=bam_path,
        total_reads=total,
        mapped_reads=mapped,
        unmapped_reads=unmapped,
        references=references,
    )


def get_read_details(bam_path: Path, read_index: int, *, contig: str | None = None) -> str:
    validate_bam(bam_path)
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        read = _find_read_by_index(bam, read_index, contig=contig)
        if read is not None:
            tags = ", ".join(f"{tag}:{value}" for tag, value in read.get_tags()) or "(none)"
            return (
                f"QNAME: {read.query_name}\n"
                f"FLAG: {read.flag}\n"
                f"RNAME: {read.reference_name}\n"
                f"POS: {(read.reference_start + 1) if read.reference_start >= 0 else 0}\n"
                f"MAPQ: {read.mapping_quality}\n"
                f"CIGAR: {read.cigarstring}\n"
                f"IS REVERSE: {read.is_reverse}\n"
                f"QUERY LEN: {read.query_length}\n"
                f"TAGS: {tags}\n"
                f"SEQ: {read.query_sequence or ''}\n"
                f"QUAL: {read.qual or ''}"
            )
    return "Read details not found."


def render_reference_panel(
    bam_path: Path,
    read_index: int,
    *,
    reference_fasta: Path | None = None,
    contig: str | None = None,
    reference_name: str | None = None,
    center_pos: int | None = None,
    flank: int = 35,
    max_rows: int = 12,
    min_mapq: int = 0,
    pan_bp: int = 0,
) -> str:
    validate_bam(bam_path)
    if flank < 0:
        raise ValueError("flank must be >= 0.")
    if max_rows <= 0:
        raise ValueError("max_rows must be > 0.")
    if min_mapq < 0:
        raise ValueError("min_mapq must be >= 0.")

    resolved_reference = _discover_reference_fasta(
        bam_path=bam_path,
        reference_fasta=reference_fasta,
    )

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        selected_read_name: str | None = None
        if reference_name is not None and center_pos is not None:
            if reference_name not in bam.references:
                return f"Alignment panel unavailable: reference '{reference_name}' was not found in BAM."
            ref_name = reference_name
            focus_center = center_pos
            focus_label = f"reference-centered at {ref_name}:{focus_center + 1}"
        else:
            target = _find_read_by_index(bam, read_index, contig=contig)
            if target is None:
                return "Alignment panel unavailable: selected read was not found."
            if target.is_unmapped or target.reference_name is None:
                return "Alignment panel unavailable: selected read is unmapped."
            if target.reference_end is None:
                return "Alignment panel unavailable: selected read has no reference end."
            ref_name = target.reference_name
            focus_center = (target.reference_start + target.reference_end) // 2
            selected_read_name = target.query_name
            focus_label = f"selected read: {target.query_name}"

        ref_length = bam.get_reference_length(ref_name)
        effective_center = focus_center + pan_bp
        effective_center = max(0, min(effective_center, max(0, ref_length - 1)))
        window_width = max(10, (2 * flank) + 1)
        window_start = max(0, min(effective_center - flank, max(0, ref_length - window_width)))
        window_end = min(ref_length, window_start + window_width)

        try:
            overlapping = list(bam.fetch(ref_name, window_start, window_end))
        except ValueError:
            return (
                "Alignment panel needs an indexed BAM (.bai). "
                "Create one with: samtools index your_file.bam"
            )

        mapped_overlapping = [
            read
            for read in overlapping
            if not read.is_unmapped
            and read.reference_start < window_end
            and (read.reference_end or read.reference_start) > window_start
            and read.mapping_quality >= min_mapq
        ]
        mapped_overlapping.sort(key=lambda read: (read.reference_start, read.query_name or ""))

        if selected_read_name is not None:
            target_read = next(
                (
                    read
                    for read in mapped_overlapping
                    if (read.query_name or "") == (selected_read_name or "")
                    and read.reference_start <= effective_center
                    and (read.reference_end or read.reference_start) >= effective_center
                ),
                None,
            )
            if target_read is not None:
                reads = _pick_alignment_reads(target_read, mapped_overlapping, max_rows=max_rows)
            else:
                reads = mapped_overlapping[:max_rows]
        else:
            reads = mapped_overlapping[:max_rows]

        insertion_lengths = _collect_insertion_lengths(
            reads=reads,
            window_start=window_start,
            window_end=window_end,
        )
        layout = _build_alignment_layout(
            window_start=window_start,
            window_end=window_end,
            insertion_lengths=insertion_lengths,
        )

        ref_seq, reference_note = _load_reference_window(
            reference_fasta=resolved_reference,
            reference_name=ref_name,
            start=window_start,
            end=window_end,
        )
        if not reads:
            reference_note = (
                "No mapped reads in this window."
                if not reference_note
                else f"No mapped reads in this window. {reference_note}"
            )
        if reference_fasta is None and resolved_reference is not None:
            prefix = f"Auto-detected reference: {resolved_reference.name}."
            reference_note = f"{prefix} {reference_note}".strip()
        return _format_alignment_text(
            focus_label=focus_label,
            selected_read_name=selected_read_name,
            reads=reads,
            ref_name=ref_name,
            ref_seq=ref_seq,
            reference_note=reference_note,
            layout=layout,
            min_mapq=min_mapq,
        )


def _iter_reads(
    bam: pysam.AlignmentFile, *, contig: str | None = None
) -> Iterable[pysam.AlignedSegment]:
    if contig:
        if contig not in bam.references:
            raise ValueError(f"Contig '{contig}' was not found in BAM header.")
        try:
            return bam.fetch(contig)
        except ValueError as exc:
            raise RuntimeError(
                "Contig filtering requires an indexed BAM (.bai). "
                "Create one with: samtools index your_file.bam"
            ) from exc
    return bam.fetch(until_eof=True)


def _find_read_by_index(
    bam: pysam.AlignmentFile,
    read_index: int,
    *,
    contig: str | None = None,
) -> pysam.AlignedSegment | None:
    reads = _iter_reads(bam, contig=contig)
    for idx, read in enumerate(reads):
        if idx == read_index:
            return read
    return None


def _load_reference_window(
    *,
    reference_fasta: Path | None,
    reference_name: str,
    start: int,
    end: int,
) -> tuple[str, str]:
    window_len = max(0, end - start)
    if window_len == 0:
        return "", ""
    if reference_fasta is None:
        return (
            "N" * window_len,
            "Reference FASTA not provided. Use --reference <ref.fa> and a .fai index.",
        )

    validate_reference_fasta(reference_fasta)
    try:
        _ensure_fasta_index(reference_fasta)
        with pysam.FastaFile(str(reference_fasta)) as fasta:
            resolved_name = _resolve_reference_name(reference_name, fasta.references)
            if resolved_name is None:
                return (
                    "N" * window_len,
                    (
                        f"Reference '{reference_name}' is missing from {reference_fasta.name}. "
                        "Contig names must match (for example chr1 vs 1)."
                    ),
                )
            seq = fasta.fetch(resolved_name, start, end).upper()
            if len(seq) < window_len:
                seq = seq + ("N" * (window_len - len(seq)))
            note = ""
            if resolved_name != reference_name:
                note = f"Using FASTA contig '{resolved_name}' for BAM contig '{reference_name}'."
            return seq, note
    except (ValueError, OSError, RuntimeError):
        return (
            "N" * window_len,
            (
                f"Could not read reference FASTA {reference_fasta.name}. "
                "Ensure it is indexed (.fai), e.g. samtools faidx ref.fa."
            ),
        )


def _discover_reference_fasta(
    *,
    bam_path: Path,
    reference_fasta: Path | None,
) -> Path | None:
    if reference_fasta is not None:
        return reference_fasta

    bam_dir = bam_path.parent
    bam_stem = bam_path.stem
    candidates = [
        bam_dir / f"{bam_stem}.fa",
        bam_dir / f"{bam_stem}.fasta",
        bam_dir / f"{bam_stem}.fna",
        bam_dir / "reference.fa",
        bam_dir / "reference.fasta",
        bam_dir / "ref.fa",
        bam_dir / "ref.fasta",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def _ensure_fasta_index(reference_fasta: Path) -> None:
    fai_path = reference_fasta.with_suffix(reference_fasta.suffix + ".fai")
    if fai_path.exists():
        return
    try:
        pysam.faidx(str(reference_fasta))
    except (ValueError, OSError, RuntimeError) as exc:
        raise RuntimeError(
            f"Could not create FASTA index for {reference_fasta}. "
            "Run: samtools faidx <reference.fa>"
        ) from exc


def _resolve_reference_name(
    bam_ref_name: str,
    fasta_references: tuple[str, ...],
) -> str | None:
    if bam_ref_name in fasta_references:
        return bam_ref_name

    # Common synonyms between BAM and FASTA conventions.
    candidates = [bam_ref_name]
    if bam_ref_name.startswith("chr"):
        candidates.append(bam_ref_name[3:])
    else:
        candidates.append(f"chr{bam_ref_name}")

    mt_aliases = {"M", "MT", "chrM", "chrMT"}
    if bam_ref_name in mt_aliases:
        candidates.extend(mt_aliases)

    # Case-insensitive fallback map.
    lower_map = {name.lower(): name for name in fasta_references}
    for candidate in candidates:
        if candidate in fasta_references:
            return candidate
        mapped = lower_map.get(candidate.lower())
        if mapped is not None:
            return mapped

    return None


@dataclass(slots=True)
class AlignmentLayout:
    window_start: int
    window_end: int
    width: int
    base_columns: dict[int, int]
    insertion_starts: dict[int, int]
    insertion_lengths: dict[int, int]


def _pick_alignment_reads(
    target: pysam.AlignedSegment,
    overlapping: list[pysam.AlignedSegment],
    *,
    max_rows: int,
) -> list[pysam.AlignedSegment]:
    selected: list[pysam.AlignedSegment] = [target]
    for read in overlapping:
        if _same_alignment(read, target):
            continue
        selected.append(read)
        if len(selected) >= max_rows:
            break
    return selected


def _collect_insertion_lengths(
    *,
    reads: list[pysam.AlignedSegment],
    window_start: int,
    window_end: int,
) -> dict[int, int]:
    insertion_lengths: dict[int, int] = {window_start - 1: 0}
    for read in reads:
        ref_pos = read.reference_start
        for op, length in read.cigartuples or []:
            if op in (0, 7, 8):
                ref_pos += length
                continue
            if op == 1:
                anchor = ref_pos - 1
                if window_start - 1 <= anchor < window_end:
                    insertion_lengths[anchor] = max(insertion_lengths.get(anchor, 0), length)
                continue
            if op in (2, 3):
                ref_pos += length
                continue
            if op in (4, 5, 6):
                continue
    return insertion_lengths


def _build_alignment_layout(
    *,
    window_start: int,
    window_end: int,
    insertion_lengths: dict[int, int],
) -> AlignmentLayout:
    base_columns: dict[int, int] = {}
    insertion_starts: dict[int, int] = {}
    insertion_cols = {k: v for k, v in insertion_lengths.items() if v > 0}
    cursor = 0

    pre_window_key = window_start - 1
    pre_insertions = insertion_cols.get(pre_window_key, 0)
    if pre_insertions:
        insertion_starts[pre_window_key] = cursor
        cursor += pre_insertions

    for ref_pos in range(window_start, window_end):
        base_columns[ref_pos] = cursor
        cursor += 1
        ins_after = insertion_cols.get(ref_pos, 0)
        if ins_after:
            insertion_starts[ref_pos] = cursor
            cursor += ins_after

    return AlignmentLayout(
        window_start=window_start,
        window_end=window_end,
        width=cursor,
        base_columns=base_columns,
        insertion_starts=insertion_starts,
        insertion_lengths=insertion_cols,
    )


def _format_alignment_text(
    *,
    focus_label: str,
    selected_read_name: str | None,
    reads: list[pysam.AlignedSegment],
    ref_name: str,
    ref_seq: str,
    reference_note: str,
    layout: AlignmentLayout,
    min_mapq: int,
) -> str:
    labels = [
        f"Region: {ref_name}:{layout.window_start + 1}-{layout.window_end}",
        "Legend: uppercase=match lowercase=mismatch '-'=deletion '~'=skipped insertion=between columns",
        f"Focus: {focus_label}  Reads shown: {len(reads)} (MAPQ >= {min_mapq})",
    ]
    if reference_note:
        labels.append(f"Note: {reference_note}")

    coord_line, tick_line = _build_coordinate_lines(layout)
    ref_track = _render_reference_track(ref_seq=ref_seq, layout=layout)
    lines = labels.copy()

    # Reads stacked first (IGV-like), reference is pinned to the bottom of the block.
    for idx, read in enumerate(reads):
        aligned = _render_read_track(read=read, ref_seq=ref_seq, layout=layout)
        is_selected = (
            selected_read_name is not None and (read.query_name or "") == selected_read_name
        )
        prefix = "TARGET" if is_selected else "READ"
        name = read.query_name or "<unnamed>"
        lines.append(f"{prefix:<7}{name[:9]:<9}{aligned}  mq={read.mapping_quality}")

    lines.append(f"{'':<16}{coord_line}")
    lines.append(f"{'':<16}{tick_line}")
    lines.append(f"{'REF':<16}{ref_track}")

    return "\n".join(lines)


def _build_coordinate_lines(layout: AlignmentLayout) -> tuple[str, str]:
    coords = [" "] * layout.width
    ticks = [" "] * layout.width
    for pos in range(layout.window_start, layout.window_end):
        col = layout.base_columns[pos]
        one_based = pos + 1
        if one_based % 10 == 0:
            ticks[col] = "|"
            coord_text = str(one_based)
            left = max(0, col - len(coord_text) + 1)
            for j, char in enumerate(coord_text):
                idx = left + j
                if idx < layout.width:
                    coords[idx] = char
        elif one_based % 5 == 0:
            ticks[col] = ":"
    return "".join(coords), "".join(ticks)


def _render_reference_track(*, ref_seq: str, layout: AlignmentLayout) -> str:
    chars = [" "] * layout.width
    for pos in range(layout.window_start, layout.window_end):
        col = layout.base_columns[pos]
        idx = pos - layout.window_start
        base = ref_seq[idx] if 0 <= idx < len(ref_seq) else "N"
        chars[col] = base.upper()
    return "".join(chars)


def _render_read_track(
    *,
    read: pysam.AlignedSegment,
    ref_seq: str,
    layout: AlignmentLayout,
) -> str:
    chars = [" "] * layout.width
    seq = read.query_sequence or ""
    ref_pos = read.reference_start
    query_pos = 0
    for op, length in read.cigartuples or []:
        if op in (0, 7, 8):
            for i in range(length):
                if layout.window_start <= ref_pos < layout.window_end:
                    base_index = query_pos + i
                    base = seq[base_index] if base_index < len(seq) else "N"
                    ref_index = ref_pos - layout.window_start
                    ref_base = ref_seq[ref_index] if 0 <= ref_index < len(ref_seq) else "N"
                    col = layout.base_columns[ref_pos]
                    if ref_base.upper() == "N":
                        chars[col] = base.upper()
                    else:
                        chars[col] = (
                            base.upper() if base.upper() == ref_base.upper() else base.lower()
                        )
                ref_pos += 1
            query_pos += length
            continue

        if op == 1:
            ins_seq = seq[query_pos : query_pos + length]
            anchor = ref_pos - 1
            start_col = layout.insertion_starts.get(anchor)
            ins_slots = layout.insertion_lengths.get(anchor, 0)
            if start_col is not None and ins_slots > 0:
                for i in range(min(length, ins_slots)):
                    base = ins_seq[i] if i < len(ins_seq) else "N"
                    chars[start_col + i] = base.lower()
                if length > ins_slots:
                    chars[start_col + ins_slots - 1] = "+"
            query_pos += length
            continue

        if op in (2, 3):
            for _ in range(length):
                if layout.window_start <= ref_pos < layout.window_end:
                    col = layout.base_columns[ref_pos]
                    chars[col] = "-" if op == 2 else "~"
                ref_pos += 1
            continue

        if op == 4:
            query_pos += length
            continue

        if op in (5, 6):
            continue

    return "".join(chars)


def _same_alignment(left: pysam.AlignedSegment, right: pysam.AlignedSegment) -> bool:
    return (
        left.reference_name == right.reference_name
        and left.reference_start == right.reference_start
        and left.reference_end == right.reference_end
        and left.query_name == right.query_name
    )
