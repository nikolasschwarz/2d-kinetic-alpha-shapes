"""Generate a bilinear UV debug texture.

Corner colors (image / UV space, v=0 at top):
  top-left     (u=0, v=0) -> black
  top-right    (u=1, v=0) -> blue
  bottom-left  (u=0, v=1) -> red
  bottom-right (u=1, v=1) -> white

Everywhere else is bilinearly interpolated between those corners.

A dark-green grid is overlaid for UV orientation checks:
  outer lines = 1px, inner lines = 2px, cells = 62x62.
When tiled, adjacent outer 1px edges meet to form a 2px stroke.
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


def make_texture(
    width: int,
    height: int,
    grid_divisions: int = 8,
    cell_size: int = 62,
    outer_line_width: int = 1,
    inner_line_width: int = 2,
) -> np.ndarray:
    expected = (
        2 * outer_line_width
        + grid_divisions * cell_size
        + (grid_divisions - 1) * inner_line_width
    )
    if width != expected or height != expected:
        raise ValueError(
            f"Size {width}x{height} does not match grid layout "
            f"({grid_divisions} cells of {cell_size}px, outer={outer_line_width}px, "
            f"inner={inner_line_width}px → expected {expected}px)."
        )

    # Pixel centers in [0, 1]
    u = (np.arange(width, dtype=np.float64) + 0.5) / width
    v = (np.arange(height, dtype=np.float64) + 0.5) / height
    uu, vv = np.meshgrid(u, v)

    # Bilinear blend of the four corners
    top = (1.0 - uu)[..., None] * TOP_LEFT + uu[..., None] * TOP_RIGHT
    bottom = (1.0 - uu)[..., None] * BOTTOM_LEFT + uu[..., None] * BOTTOM_RIGHT
    rgb = (1.0 - vv)[..., None] * top + vv[..., None] * bottom

    # Layout: [outer 1px][cell][inner 2px][cell]...[cell][outer 1px]
    # When tiled, adjacent outer 1px lines meet to form a 2px stroke like the inners.
    def paint_axis_lines(paint_vertical: bool) -> None:
        # Outer start
        if paint_vertical:
            rgb[:, 0:outer_line_width, :] = GRID_COLOR
        else:
            rgb[0:outer_line_width, :, :] = GRID_COLOR

        cursor = outer_line_width
        for cell_index in range(grid_divisions):
            cursor += cell_size
            if cell_index < grid_divisions - 1:
                if paint_vertical:
                    rgb[:, cursor : cursor + inner_line_width, :] = GRID_COLOR
                else:
                    rgb[cursor : cursor + inner_line_width, :, :] = GRID_COLOR
                cursor += inner_line_width

        # Outer end
        if paint_vertical:
            rgb[:, width - outer_line_width : width, :] = GRID_COLOR
        else:
            rgb[height - outer_line_width : height, :, :] = GRID_COLOR

    paint_axis_lines(paint_vertical=True)
    paint_axis_lines(paint_vertical=False)

    return np.clip(np.rint(rgb * 255.0), 0, 255).astype(np.uint8)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path(__file__).with_name("uv_debug_texture.png"),
        help="Output PNG path",
    )
    parser.add_argument(
        "-s",
        "--size",
        type=int,
        default=512,
        help="Square texture size in pixels (must match grid layout)",
    )
    parser.add_argument("--width", type=int, default=None, help="Override width")
    parser.add_argument("--height", type=int, default=None, help="Override height")
    parser.add_argument("--grid", type=int, default=8, help="Number of grid cells per axis")
    parser.add_argument("--cell-size", type=int, default=62, help="Interior cell size in pixels")
    parser.add_argument("--outer-width", type=int, default=1, help="Outer grid line thickness")
    parser.add_argument("--inner-width", type=int, default=2, help="Inner grid line thickness")
    args = parser.parse_args()

    width = args.width or args.size
    height = args.height or args.size
    image = Image.fromarray(
        make_texture(
            width,
            height,
            grid_divisions=args.grid,
            cell_size=args.cell_size,
            outer_line_width=args.outer_width,
            inner_line_width=args.inner_width,
        ),
        mode="RGB",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    image.save(args.output)
    print(f"Wrote {args.output.resolve()} ({width}x{height})")


if __name__ == "__main__":
    main()
