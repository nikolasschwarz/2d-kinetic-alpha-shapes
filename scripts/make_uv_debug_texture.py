"""Generate a bilinear UV debug texture.

Corner colors (image / UV space, v=0 at top):
  top-left     (u=0, v=0) -> black
  top-right    (u=1, v=0) -> blue
  bottom-left  (u=0, v=1) -> red
  bottom-right (u=1, v=1) -> white

Everywhere else is bilinearly interpolated between those corners.

A dark-green grid is overlaid for UV orientation checks.
Line thicknesses cycle 6,2,4,2,... along each axis. Each line is shared by
adjacent cells (half counted toward each side), so every cell spans a fixed
pitch. Default: 32 cells x 64px = 2048px, with e.g. a 6|2 cell laid out as
3px line + 60px base + 1px line.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from PIL import Image

# RGB in [0, 1]
TOP_LEFT = np.array([0.0, 0.0, 0.0])  # black
TOP_RIGHT = np.array([0.0, 0.0, 1.0])  # blue
BOTTOM_LEFT = np.array([1.0, 0.0, 0.0])  # red
BOTTOM_RIGHT = np.array([1.0, 1.0, 1.0])  # white
GRID_COLOR = np.array([0.0, 0.35, 0.0])  # dark green

DEFAULT_GRID = 32
DEFAULT_CELL_PITCH = 64
DEFAULT_LINE_PATTERN = (6, 2, 4, 2)


def expected_size(grid_divisions: int, cell_pitch: int) -> int:
    return grid_divisions * cell_pitch


def make_texture(
    width: int,
    height: int,
    grid_divisions: int = DEFAULT_GRID,
    cell_pitch: int = DEFAULT_CELL_PITCH,
    line_pattern: tuple[int, ...] = DEFAULT_LINE_PATTERN,
) -> np.ndarray:
    if not line_pattern:
        raise ValueError("line pattern must not be empty")
    if any(w <= 0 or w % 2 != 0 for w in line_pattern):
        raise ValueError(f"line widths must be positive even integers so halves are integral: {line_pattern}")

    expected = expected_size(grid_divisions, cell_pitch)
    if width != expected or height != expected:
        raise ValueError(
            f"Size {width}x{height} does not match grid layout "
            f"({grid_divisions} cells x {cell_pitch}px pitch → expected {expected}px)."
        )

    # Pixel centers in [0, 1]
    u = (np.arange(width, dtype=np.float64) + 0.5) / width
    v = (np.arange(height, dtype=np.float64) + 0.5) / height
    uu, vv = np.meshgrid(u, v)

    # Bilinear blend of the four corners
    top = (1.0 - uu)[..., None] * TOP_LEFT + uu[..., None] * TOP_RIGHT
    bottom = (1.0 - uu)[..., None] * BOTTOM_LEFT + uu[..., None] * BOTTOM_RIGHT
    rgb = (1.0 - vv)[..., None] * top + vv[..., None] * bottom

    # One line per cell boundary. Line k (width pattern[k]) is centered on k*pitch.
    # The wrap line (k=0) is split: half at the start, half at the end → seamless tile.
    def paint_axis_lines(paint_vertical: bool) -> None:
        for boundary_index in range(grid_divisions):
            line_width = line_pattern[boundary_index % len(line_pattern)]
            half = line_width // 2
            center = boundary_index * cell_pitch

            if boundary_index == 0:
                # Right half at the left edge, left half at the right edge.
                segments = [(0, half), (width - half, width)] if paint_vertical else [
                    (0, half),
                    (height - half, height),
                ]
            else:
                segments = [(center - half, center + half)]

            for start, end in segments:
                if paint_vertical:
                    rgb[:, start:end, :] = GRID_COLOR
                else:
                    rgb[start:end, :, :] = GRID_COLOR

    paint_axis_lines(paint_vertical=True)
    paint_axis_lines(paint_vertical=False)

    return np.clip(np.rint(rgb * 255.0), 0, 255).astype(np.uint8)


def parse_pattern(text: str) -> tuple[int, ...]:
    parts = [p.strip() for p in text.split(",") if p.strip()]
    if not parts:
        raise argparse.ArgumentTypeError("line pattern must contain at least one width")
    return tuple(int(p) for p in parts)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path(__file__).with_name("uv_debug_texture.png"),
        help="Output PNG path",
    )
    parser.add_argument("--grid", type=int, default=DEFAULT_GRID, help="Number of grid cells per axis")
    parser.add_argument(
        "--cell-pitch",
        type=int,
        default=DEFAULT_CELL_PITCH,
        help="Full cell size in pixels including half of each bordering line",
    )
    parser.add_argument(
        "--line-pattern",
        type=parse_pattern,
        default=DEFAULT_LINE_PATTERN,
        help="Comma-separated repeating line widths (even), e.g. 6,2,4,2",
    )
    parser.add_argument(
        "-s",
        "--size",
        type=int,
        default=None,
        help="Square texture size (optional; defaults to grid * cell-pitch)",
    )
    parser.add_argument("--width", type=int, default=None, help="Override width")
    parser.add_argument("--height", type=int, default=None, help="Override height")
    args = parser.parse_args()

    computed = expected_size(args.grid, args.cell_pitch)
    width = args.width or args.size or computed
    height = args.height or args.size or computed
    image = Image.fromarray(
        make_texture(
            width,
            height,
            grid_divisions=args.grid,
            cell_pitch=args.cell_pitch,
            line_pattern=args.line_pattern,
        ),
        mode="RGB",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    image.save(args.output)
    print(
        f"Wrote {args.output.resolve()} ({width}x{height}), "
        f"{args.grid} cells x {args.cell_pitch}px pitch, lines {args.line_pattern}"
    )


if __name__ == "__main__":
    main()
