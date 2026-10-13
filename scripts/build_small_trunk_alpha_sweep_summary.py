#!/usr/bin/env python3
"""Rebuild Small Trunk alpha-sweep summary CSV + LaTeX table from per-alpha stats CSVs.

Usage:
  python build_small_trunk_alpha_sweep_summary.py \\
      --input-dir path/to/EvoEngine_App \\
      [--expected-alphas 2500,900,400,100,25,9,4,1] \\
      [--output-dir path/to/out]

Rod / segment count is strands + subdivision events (or the CSV segment_count column when present).
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Optional


ALPHA_FILENAME_RE = re.compile(
    r"meshing_statistics_Small_Trunk_alpha_(?P<alpha>\d+)_(?P<sections>\d+)_(?P<stamp>.+)\.csv$",
    re.IGNORECASE,
)

# Matches the C++ summary CSV produced by DynamicStrandsDemo::RunSmallTrunkAlphaSweepMeshing.
SUMMARY_HEADER = [
    "alpha",
    "cutoff",
    "succeeded",
    "section_count",
    "runtime_s",
    "strand_count",
    "branch_count",
    "segment_count",
    "subdivision",
    "section",
    "separation",
    "flip",
    "radius",
    "crossing",
    "alpha_recorded",
    "triangle_count",
    "vertex_count",
    "failure",
]


@dataclass
class TotalsRow:
    alpha: float
    cutoff: float
    succeeded: bool
    section_count: int = 0
    runtime_s: float = 0.0
    strand_count: Optional[int] = None
    branch_count: Optional[int] = None
    segment_count: Optional[int] = None
    subdivision: int = 0
    section: int = 0
    separation: int = 0
    flip: int = 0
    radius: int = 0
    crossing: int = 0
    alpha_recorded: Optional[float] = None
    triangle_count: Optional[int] = None
    vertex_count: Optional[int] = None
    failure: str = ""
    source: str = ""

    def resolved_segment_count(self) -> Optional[int]:
        if self.segment_count is not None:
            return self.segment_count
        if self.strand_count is None:
            return None
        return self.strand_count + self.subdivision


def parse_optional_int(value: str) -> Optional[int]:
    value = (value or "").strip()
    if not value:
        return None
    return int(float(value))


def parse_optional_float(value: str) -> Optional[float]:
    value = (value or "").strip()
    if not value:
        return None
    return float(value)


def thin_space_int(n: int) -> str:
    s = f"{n:d}"
    parts = []
    while s:
        parts.append(s[-3:])
        s = s[:-3]
    return "\\,".join(reversed(parts))


def format_runtime(seconds: float) -> str:
    # Match the sample style: two decimals.
    return f"{seconds:.2f}"


def format_alpha(alpha: float) -> str:
    if abs(alpha - round(alpha)) < 1e-9:
        return f"{int(round(alpha))}.0"
    return f"{alpha:g}"


def load_totals_from_csv(path: Path) -> TotalsRow:
    match = ALPHA_FILENAME_RE.match(path.name)
    if not match:
        raise ValueError(f"Unexpected alpha stats filename: {path.name}")

    alpha_from_name = float(match.group("alpha"))
    sections_from_name = int(match.group("sections"))

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        total = None
        section_rows = 0
        for row in reader:
            section_id = (row.get("section_id") or "").strip()
            if section_id == "total":
                total = row
            elif section_id:
                section_rows += 1

    if total is None:
        raise ValueError(f"No totals row in {path}")

    alpha_recorded = parse_optional_float(total.get("alpha", ""))
    alpha = alpha_recorded if alpha_recorded is not None else alpha_from_name
    failure = (total.get("failure") or "").strip().strip('"')
    succeeded = not failure

    return TotalsRow(
        alpha=alpha,
        cutoff=math.sqrt(alpha),
        succeeded=succeeded,
        section_count=section_rows if section_rows > 0 else sections_from_name,
        runtime_s=float(total["runtime_s"]),
        strand_count=parse_optional_int(total.get("strand_count", "")),
        branch_count=parse_optional_int(total.get("branch_count", "")),
        segment_count=parse_optional_int(total.get("segment_count", "")),
        subdivision=parse_optional_int(total.get("subdivision", "")) or 0,
        section=parse_optional_int(total.get("section", "")) or 0,
        separation=parse_optional_int(total.get("separation", "")) or 0,
        flip=parse_optional_int(total.get("flip", "")) or 0,
        radius=parse_optional_int(total.get("radius", "")) or 0,
        crossing=parse_optional_int(total.get("crossing", "")) or 0,
        alpha_recorded=alpha_recorded if alpha_recorded is not None else alpha,
        triangle_count=parse_optional_int(total.get("triangle_count", "")),
        vertex_count=parse_optional_int(total.get("vertex_count", "")),
        failure=failure,
        source=path.name,
    )


def collect_rows(input_dir: Path, expected_alphas: list[float]) -> list[TotalsRow]:
    files = sorted(input_dir.glob("meshing_statistics_Small_Trunk_alpha_*.csv"))
    # Prefer per-alpha stats files only (skip summaries / event lists).
    files = [
        p
        for p in files
        if "event_list" not in p.name.lower()
        and "summary" not in p.name.lower()
        and ALPHA_FILENAME_RE.match(p.name)
    ]

    by_alpha: dict[float, TotalsRow] = {}
    for path in files:
        row = load_totals_from_csv(path)
        key = float(row.alpha)
        # Keep the newest file if duplicates exist for the same alpha.
        previous = by_alpha.get(key)
        if previous is None or path.stat().st_mtime >= (input_dir / previous.source).stat().st_mtime:
            by_alpha[key] = row

    rows: list[TotalsRow] = []
    for alpha in expected_alphas:
        key = float(alpha)
        if key in by_alpha:
            rows.append(by_alpha[key])
        else:
            rows.append(
                TotalsRow(
                    alpha=key,
                    cutoff=math.sqrt(key),
                    succeeded=False,
                    alpha_recorded=key,
                    failure="missing CSV (likely crash / nullptr before stats write)",
                    source="",
                )
            )
    return rows


def write_summary_csv(path: Path, rows: list[TotalsRow]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_HEADER)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "alpha": row.alpha,
                    "cutoff": row.cutoff,
                    "succeeded": 1 if row.succeeded else 0,
                    "section_count": row.section_count,
                    "runtime_s": row.runtime_s,
                    "strand_count": "" if row.strand_count is None else row.strand_count,
                    "branch_count": "" if row.branch_count is None else row.branch_count,
                    "segment_count": "" if row.resolved_segment_count() is None else row.resolved_segment_count(),
                    "subdivision": row.subdivision,
                    "section": row.section,
                    "separation": row.separation,
                    "flip": row.flip,
                    "radius": row.radius,
                    "crossing": row.crossing,
                    "alpha_recorded": "" if row.alpha_recorded is None else row.alpha_recorded,
                    "triangle_count": "" if row.triangle_count is None else row.triangle_count,
                    "vertex_count": "" if row.vertex_count is None else row.vertex_count,
                    "failure": row.failure,
                }
            )


def latex_cell_int(value: Optional[int]) -> str:
    if value is None:
        return "---"
    return thin_space_int(value)


def write_latex_table(path: Path, rows: list[TotalsRow], figure_placeholder: str) -> None:
    lines = [
        "% Auto-generated by build_small_trunk_alpha_sweep_summary.py",
        "% Requires: \\usepackage{booktabs,makecell}",
        "% Rod elements = strand_count + subdivision (or CSV segment_count when present).",
        "\\begin{tabular}{rrrrrrrrr}",
        "\\toprule",
        "\\textbf{Figure} & \\textbf{\\# Strands} & \\makecell[r]{\\textbf{\\# Rod}\\\\ \\textbf{elements}} & "
        "\\textbf{$\\alpha$} & \\textbf{\\# Sections} & \\makecell[r]{\\textbf{\\# Flip}\\\\ \\textbf{Events}} & "
        "\\makecell[r]{\\textbf{\\# Radius}\\\\ \\textbf{Events}} & "
        "\\makecell[r]{\\textbf{\\# Crossing}\\\\ \\textbf{Events}} & \\textbf{Runtime ($s$)} \\\\",
        "\\midrule",
    ]

    # Present ascending alpha for the paper table.
    for row in sorted(rows, key=lambda r: r.alpha):
        rods = latex_cell_int(row.resolved_segment_count())
        if row.succeeded:
            runtime = format_runtime(row.runtime_s)
            strands = latex_cell_int(row.strand_count)
            sections = latex_cell_int(row.section_count)
            flip = latex_cell_int(row.flip)
            radius = latex_cell_int(row.radius)
            crossing = latex_cell_int(row.crossing)
        else:
            runtime = "---"
            strands = latex_cell_int(row.strand_count) if row.strand_count is not None else "---"
            sections = latex_cell_int(row.section_count) if row.section_count else "---"
            flip = "---"
            radius = "---"
            crossing = "---"

        lines.append(
            f"{figure_placeholder} & {strands} & {rods} & {format_alpha(row.alpha)} & "
            f"{sections} & {flip} & {radius} & {crossing} & {runtime} \\\\"
        )

    lines.extend(["\\bottomrule", "\\end{tabular}", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path(__file__).resolve().parents[4]
        / "out"
        / "build"
        / "x64-Release"
        / "EvoEngine_App",
        help="Directory containing meshing_statistics_Small_Trunk_alpha_*.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Where to write summary CSV + LaTeX (default: input-dir)",
    )
    parser.add_argument(
        "--expected-alphas",
        type=str,
        default="2500,900,400,100,25,9,4,1",
        help="Comma-separated alpha list (missing CSVs become failure rows)",
    )
    parser.add_argument(
        "--figure-placeholder",
        type=str,
        default=r"Fig.~\ref{fig:TODO}",
        help="LaTeX figure cell shared by all rows",
    )
    args = parser.parse_args()

    input_dir: Path = args.input_dir
    output_dir: Path = args.output_dir or input_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    expected_alphas = [float(x.strip()) for x in args.expected_alphas.split(",") if x.strip()]
    rows = collect_rows(input_dir, expected_alphas)

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    summary_path = output_dir / f"meshing_statistics_Small_Trunk_alpha_sweep_summary_{stamp}.csv"
    latex_path = output_dir / f"meshing_statistics_Small_Trunk_alpha_sweep_table_{stamp}.tex"
    latest_summary = output_dir / "meshing_statistics_Small_Trunk_alpha_sweep_summary.csv"
    latest_latex = output_dir / "meshing_statistics_Small_Trunk_alpha_sweep_table.tex"

    write_summary_csv(summary_path, rows)
    write_summary_csv(latest_summary, rows)
    write_latex_table(latex_path, rows, args.figure_placeholder)
    write_latex_table(latest_latex, rows, args.figure_placeholder)

    print(f"Wrote {summary_path}")
    print(f"Wrote {latest_summary}")
    print(f"Wrote {latex_path}")
    print(f"Wrote {latest_latex}")
    print()
    print("Rows:")
    for row in rows:
        status = "ok" if row.succeeded else f"FAIL ({row.failure})"
        segs = row.resolved_segment_count()
        segs_s = "---" if segs is None else str(segs)
        print(
            f"  alpha={row.alpha:g}  strands={row.strand_count}  segments={segs_s}  "
            f"runtime={row.runtime_s:.2f}s  flip={row.flip} radius={row.radius} "
            f"crossing={row.crossing}  [{status}]"
        )


if __name__ == "__main__":
    main()
