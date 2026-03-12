from __future__ import annotations

from pathlib import Path

from textual import on
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical
from textual.widgets import DataTable, Footer, Header, Static

from tbam.bam_utils import (
    fetch_rows,
    get_header_text,
    get_read_details,
    render_reference_panel,
)


class BamViewerApp(App[None]):
    BASE_BG_COLORS = {
        "A": "#b22222",   # dark red
        "C": "#1e66f5",   # blue
        "G": "#2e7d32",   # green
        "T": "#b36b00",   # dark orange
        "N": "#555555",   # gray
    }
    MISMATCH_BG_COLOR = "#666666"

    CSS = """
    Screen {
        layout: vertical;
    }

    #top {
        height: 8;
    }

    #table {
        width: 3fr;
        border: round #666666;
    }

    #detail {
        width: 2fr;
        border: round #666666;
        padding: 1 1;
    }

    #alignment_panel {
        height: 1fr;
        border: round #666666;
        padding: 1 1;
    }

    #align_meta {
        height: 2;
    }

    #align_reads {
        height: 1fr;
        border: round #666666;
        text-wrap: nowrap;
        overflow: auto auto;
        content-align: left bottom;
    }

    #align_ref {
        height: 5;
        border: round #666666;
        text-wrap: nowrap;
        overflow: auto hidden;
    }

    #header_view {
        height: 12;
        border: round #666666;
        padding: 1 1;
        display: none;
    }
    """

    BINDINGS = [
        Binding("q", "quit", "Quit"),
        Binding("h", "toggle_header", "Header"),
        Binding("r", "reload", "Reload"),
        Binding("[", "shrink_window", "Zoom In"),
        Binding("]", "expand_window", "Zoom Out"),
        Binding("a", "pan_left", "Pan Left"),
        Binding("d", "pan_right", "Pan Right"),
        Binding("0", "center_panel", "Center"),
    ]

    def __init__(
        self,
        bam_path: Path,
        *,
        max_reads: int = 1000,
        contig: str | None = None,
        min_mapq: int = 0,
        reference_fasta: Path | None = None,
        window_bp: int = 35,
        panel_reads: int = 40,
    ) -> None:
        super().__init__()
        self.bam_path = bam_path
        self.max_reads = max_reads
        self.contig = contig
        self.min_mapq = min_mapq
        self.reference_fasta = reference_fasta
        self.window_bp = window_bp
        self.panel_reads = panel_reads
        self.header_visible = False
        self.current_row_index: int | None = None
        self.rows_by_index: dict[int, tuple[str, int]] = {}
        self.panel_reference_name: str | None = contig
        self.panel_center_pos: int | None = None

    def compose(self) -> ComposeResult:
        yield Header()
        with Horizontal(id="top"):
            yield DataTable(id="table", zebra_stripes=True, cursor_type="row")
            yield Static("Select a row to see details.", id="detail")
        with Vertical(id="alignment_panel"):
            yield Static("Select a row to view reference alignment.", id="align_meta")
            yield Static("", id="align_reads")
            yield Static("", id="align_ref")
        yield Static("", id="header_view")
        yield Footer()

    def on_mount(self) -> None:
        table = self.query_one("#table", DataTable)
        table.add_columns("Idx", "QNAME", "RNAME", "POS", "MAPQ", "CIGAR", "FLAG", "LEN")
        self.load_rows()
        self.load_header()
        self.title = f"TBAM Viewer - {self.bam_path.name}"
        self._update_subtitle()

    def load_rows(self) -> None:
        table = self.query_one("#table", DataTable)
        table.clear(columns=False)
        self.rows_by_index = {}
        try:
            rows = fetch_rows(
                self.bam_path,
                max_reads=self.max_reads,
                contig=self.contig,
                min_mapq=self.min_mapq,
            )
        except (ValueError, RuntimeError) as exc:
            self.query_one("#detail", Static).update(f"Error: {exc}")
            self.query_one("#align_meta", Static).update(f"Error: {exc}")
            self.query_one("#align_reads", Static).update("")
            self.query_one("#align_ref", Static).update("")
            return

        for row in rows:
            self.rows_by_index[row.index] = (row.rname, row.pos)
            table.add_row(
                str(row.index),
                row.qname,
                row.rname,
                str(row.pos),
                str(row.mapq),
                row.cigar,
                str(row.flag),
                str(row.length),
                key=str(row.index),
            )
        if rows:
            table.cursor_coordinate = (0, 0)
            self.current_row_index = rows[0].index
            if self.panel_reference_name is None:
                self.panel_reference_name = rows[0].rname if rows[0].rname != "*" else None
            if self.panel_center_pos is None:
                self.panel_center_pos = max(0, rows[0].pos - 1)
            self.update_detail_for_row(rows[0].index)
            self.update_alignment_panel()
        else:
            self.current_row_index = None
            self.query_one("#detail", Static).update("No reads match current filters.")
            self.query_one("#align_meta", Static).update("No alignment view available.")
            self.query_one("#align_reads", Static).update("")
            self.query_one("#align_ref", Static).update("")

    def load_header(self) -> None:
        header_panel = self.query_one("#header_view", Static)
        try:
            header_text = get_header_text(self.bam_path).strip() or "(No header text available)"
            header_panel.update(header_text)
        except (FileNotFoundError, ValueError, RuntimeError) as exc:
            header_panel.update(f"Error: {exc}")

    @on(DataTable.RowHighlighted)
    def on_row_highlighted(self, event: DataTable.RowHighlighted) -> None:
        row_key = event.row_key.value if event.row_key else None
        if row_key is None:
            return
        row_index = int(row_key)
        self.current_row_index = row_index
        self.update_detail_for_row(row_index)

    def update_detail_for_row(self, row_index: int) -> None:
        detail = self.query_one("#detail", Static)
        try:
            detail.update(get_read_details(self.bam_path, row_index, contig=self.contig))
        except (ValueError, RuntimeError) as exc:
            detail.update(f"Error: {exc}")

    def update_alignment_panel(self) -> None:
        meta = self.query_one("#align_meta", Static)
        reads_stack = self.query_one("#align_reads", Static)
        ref = self.query_one("#align_ref", Static)
        if self.panel_reference_name is None or self.panel_center_pos is None:
            meta.update("Alignment view unavailable: no reference focus is set.")
            reads_stack.update("")
            ref.update("")
            return
        try:
            panel_text = render_reference_panel(
                self.bam_path,
                0,
                reference_fasta=self.reference_fasta,
                contig=self.contig,
                reference_name=self.panel_reference_name,
                center_pos=self.panel_center_pos,
                flank=self.window_bp,
                max_rows=self.panel_reads,
                min_mapq=self.min_mapq,
            )
            self._load_alignment_widgets(
                panel_text=panel_text,
                meta=meta,
                reads_stack=reads_stack,
                ref=ref,
            )
        except (ValueError, RuntimeError, FileNotFoundError) as exc:
            meta.update(f"Error: {exc}")
            reads_stack.update("")
            ref.update("")

    def action_toggle_header(self) -> None:
        self.header_visible = not self.header_visible
        panel = self.query_one("#header_view", Static)
        panel.styles.display = "block" if self.header_visible else "none"
        panel.refresh(layout=True)

    def action_reload(self) -> None:
        self.load_rows()
        self.load_header()

    def action_shrink_window(self) -> None:
        self.window_bp = max(5, self.window_bp - 5)
        self._update_subtitle()
        self.update_alignment_panel()

    def action_expand_window(self) -> None:
        self.window_bp = min(250, self.window_bp + 5)
        self._update_subtitle()
        self.update_alignment_panel()

    def action_pan_left(self) -> None:
        if self.panel_center_pos is not None:
            self.panel_center_pos = max(0, self.panel_center_pos - 10)
        self._update_subtitle()
        self.update_alignment_panel()

    def action_pan_right(self) -> None:
        if self.panel_center_pos is not None:
            self.panel_center_pos += 10
        self._update_subtitle()
        self.update_alignment_panel()

    def action_center_panel(self) -> None:
        if self.current_row_index is not None and self.current_row_index in self.rows_by_index:
            rname, pos = self.rows_by_index[self.current_row_index]
            if rname != "*":
                self.panel_reference_name = rname
                self.panel_center_pos = max(0, pos - 1)
        self._update_subtitle()
        self.update_alignment_panel()

    def _update_subtitle(self) -> None:
        reference_label = self.reference_fasta.name if self.reference_fasta else "none"
        active_contig = self.panel_reference_name or self.contig
        contig_label = f"contig={active_contig} " if active_contig else ""
        center_label = (
            f"center={self.panel_center_pos + 1}" if self.panel_center_pos is not None else "center=NA"
        )
        self.sub_title = (
            f"{contig_label}min_mapq={self.min_mapq} max_reads={self.max_reads} "
            f"window_bp={self.window_bp} panel_reads={self.panel_reads} "
            f"{center_label} reference={reference_label}"
        )

    def _load_alignment_widgets(
        self,
        *,
        panel_text: str,
        meta: Static,
        reads_stack: Static,
        ref: Static,
    ) -> None:
        lines = panel_text.splitlines()
        read_lines: list[str] = []
        header_lines: list[str] = []

        for line in lines:
            if line.startswith("TARGET") or line.startswith("READ"):
                read_lines.append(line)
            else:
                if line.strip():
                    header_lines.append(line)

        # Show read stack bottom-up so rows build upward from the reference track.
        reads_stack.update(self._colorize_read_block(list(reversed(read_lines))))
        reads_stack.scroll_end(animate=False)

        # Renderer always appends coordinate, tick, and REF as the last 3 lines.
        coord_line = ""
        tick_line = ""
        ref_line = ""
        if len(lines) >= 3:
            coord_line, tick_line, ref_line = lines[-3], lines[-2], lines[-1]
        if not ref_line.lstrip().startswith("REF"):
            # Fallback if formatting changes in the renderer.
            ref_line = next((line for line in reversed(lines) if line.lstrip().startswith("REF")), "")
        if ref_line:
            ref.update("\n".join([coord_line, tick_line, self._colorize_ref_line(ref_line)]).rstrip())
        else:
            ref.update("Reference track unavailable for current view.")

        if header_lines:
            meta.update("\n".join(header_lines[:3]))
        else:
            meta.update("Alignment view")

    def _colorize_read_block(self, read_lines: list[str]) -> str:
        out_lines: list[str] = []
        for line in read_lines:
            if len(line) <= 16:
                out_lines.append(line)
                continue
            prefix = line[:16]
            track_and_meta = line[16:]
            track = track_and_meta
            suffix = ""
            if "  mq=" in track_and_meta:
                track, rest = track_and_meta.rsplit("  mq=", 1)
                suffix = f"  mq={rest}"
            colored_track = self._colorize_track(track)
            out_lines.append(f"{prefix}{colored_track}{suffix}")
        return "\n".join(out_lines).rstrip()

    def _colorize_track(self, track: str) -> str:
        chars: list[str] = []
        for ch in track:
            if ch == " ":
                chars.append(ch)
                continue
            if ch in "-~+":
                chars.append(self._bg_style(ch, "#8b0000"))
                continue
            if ch.islower():
                chars.append(self._bg_style(ch.upper(), self.MISMATCH_BG_COLOR))
                continue
            base = ch.upper()
            bg = self.BASE_BG_COLORS.get(base)
            if bg is None:
                chars.append(ch)
            else:
                chars.append(self._bg_style(ch, bg))
        return "".join(chars)

    def _colorize_ref_line(self, ref_line: str) -> str:
        if len(ref_line) <= 16:
            return ref_line
        prefix = ref_line[:16]
        track = ref_line[16:]
        return f"{prefix}{self._colorize_track(track)}"

    def _bg_style(self, char: str, bg: str) -> str:
        return f"[white on {bg}]{char}[/white on {bg}]"
