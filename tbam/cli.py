from __future__ import annotations

import argparse
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="tbam",
        description="CLI BAM viewer built with Textual.",
    )
    subparsers = parser.add_subparsers(dest="command")

    view = subparsers.add_parser("view", help="Open BAM file in Textual TUI")
    view.add_argument("bam", type=Path, help="Path to BAM file")
    view.add_argument("--max-reads", type=int, default=1000, help="Max reads to load")
    view.add_argument("--contig", type=str, default=None, help="Optional contig filter")
    view.add_argument("--min-mapq", type=int, default=0, help="Minimum MAPQ")
    view.add_argument(
        "--reference",
        type=Path,
        default=None,
        help="Reference FASTA for alignment panel (optional but recommended)",
    )
    view.add_argument(
        "--window-bp",
        type=int,
        default=35,
        help="Half-window size in bp for reference-centered pileup view",
    )
    view.add_argument(
        "--panel-reads",
        type=int,
        default=40,
        help="How many mapped reads to render in the alignment panel",
    )

    inspect = subparsers.add_parser("inspect", help="Print BAM summary in terminal")
    inspect.add_argument("bam", type=Path, help="Path to BAM file")

    panel = subparsers.add_parser(
        "panel-preview",
        help="Print IGV-like reference/read panel in terminal (no TUI)",
    )
    panel.add_argument("bam", type=Path, help="Path to BAM file")
    panel.add_argument("--reference", type=Path, default=None, help="Reference FASTA")
    panel.add_argument("--read-index", type=int, default=0, help="Read index to center on")
    panel.add_argument("--contig", type=str, default=None, help="Optional contig filter")
    panel.add_argument("--window-bp", type=int, default=35, help="Flanking bases")
    panel.add_argument("--panel-reads", type=int, default=40, help="Number of reads to render")
    panel.add_argument("--min-mapq", type=int, default=0, help="Minimum MAPQ")

    demo = subparsers.add_parser("demo-data", help="Generate demo BAM + BAI + FASTA")
    demo.add_argument(
        "--out-dir",
        type=Path,
        default=Path("sample_data"),
        help="Output directory for demo BAM",
    )
    demo.add_argument("--reads", type=int, default=300, help="Approximate number of reads")
    demo.add_argument("--prefix", type=str, default="demo", help="Filename prefix")
    demo.add_argument("--open", action="store_true", help="Open viewer after generation")
    demo.add_argument("--max-reads", type=int, default=1000, help="Max reads to load in TUI")
    demo.add_argument(
        "--window-bp",
        type=int,
        default=35,
        help="Half-window size in bp for reference-centered pileup view",
    )
    demo.add_argument(
        "--panel-reads",
        type=int,
        default=40,
        help="How many mapped reads to render in the alignment panel",
    )

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "view":
        _check_dependency("pysam")
        from tbam.bam_utils import validate_bam, validate_reference_fasta

        try:
            validate_bam(args.bam)
            if args.reference is not None:
                validate_reference_fasta(args.reference)
            _launch_view(
                args.bam,
                max_reads=args.max_reads,
                contig=args.contig,
                min_mapq=args.min_mapq,
                reference_fasta=args.reference,
                window_bp=args.window_bp,
                panel_reads=args.panel_reads,
            )
            return 0
        except (FileNotFoundError, ValueError, RuntimeError) as exc:
            _print_error(str(exc))
            return 2

    if args.command == "inspect":
        _check_dependency("pysam")
        from tbam.bam_utils import summarize_bam, validate_bam

        try:
            validate_bam(args.bam)
            summary = summarize_bam(args.bam)
            print(f"BAM: {summary.path}")
            print(f"Total reads: {summary.total_reads}")
            print(f"Mapped reads: {summary.mapped_reads}")
            print(f"Unmapped reads: {summary.unmapped_reads}")
            print("References:")
            for ref_name, length in summary.references:
                print(f"  - {ref_name}: {length}")
            return 0
        except (FileNotFoundError, ValueError, RuntimeError) as exc:
            _print_error(str(exc))
            return 2

    if args.command == "panel-preview":
        _check_dependency("pysam")
        from tbam.bam_utils import render_reference_panel, validate_bam, validate_reference_fasta

        try:
            validate_bam(args.bam)
            if args.reference is not None:
                validate_reference_fasta(args.reference)
            panel_text = render_reference_panel(
                args.bam,
                args.read_index,
                reference_fasta=args.reference,
                contig=args.contig,
                flank=args.window_bp,
                max_rows=args.panel_reads,
                min_mapq=args.min_mapq,
            )
            print(panel_text)
            return 0
        except (FileNotFoundError, ValueError, RuntimeError) as exc:
            _print_error(str(exc))
            return 2

    if args.command == "demo-data":
        _check_dependency("pysam")
        from tbam.demo_data import generate_demo_bam

        try:
            bam_path, bai_path, fasta_path, fai_path = generate_demo_bam(
                args.out_dir,
                reads=args.reads,
                prefix=args.prefix,
            )
            print(f"Generated BAM: {bam_path}")
            print(f"Generated index: {bai_path}")
            print(f"Generated reference FASTA: {fasta_path}")
            print(f"Generated FASTA index: {fai_path}")
            print("Use this command to open it:")
            print(f"  tbam view {bam_path} --reference {fasta_path}")
            if args.open:
                _launch_view(
                    bam_path,
                    max_reads=args.max_reads,
                    contig=None,
                    min_mapq=0,
                    reference_fasta=fasta_path,
                    window_bp=args.window_bp,
                    panel_reads=args.panel_reads,
                )
            return 0
        except (ValueError, RuntimeError) as exc:
            _print_error(str(exc))
            return 2

    parser.print_help()
    return 1


def _launch_view(
    bam_path: Path,
    *,
    max_reads: int,
    contig: str | None,
    min_mapq: int,
    reference_fasta: Path | None,
    window_bp: int,
    panel_reads: int,
) -> None:
    _check_dependency("textual")
    from tbam.app import BamViewerApp

    app = BamViewerApp(
        bam_path=bam_path,
        max_reads=max_reads,
        contig=contig,
        min_mapq=min_mapq,
        reference_fasta=reference_fasta,
        window_bp=window_bp,
        panel_reads=panel_reads,
    )
    app.run()


def _check_dependency(module_name: str) -> None:
    try:
        __import__(module_name)
    except ModuleNotFoundError as exc:
        raise SystemExit(
            f"Missing dependency '{module_name}'. Install with: pip install -e ."
        ) from exc


def _print_error(message: str) -> None:
    print(f"Error: {message}", file=sys.stderr)


if __name__ == "__main__":
    raise SystemExit(main())
