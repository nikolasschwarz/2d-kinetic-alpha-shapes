#!/usr/bin/env python3
"""
Find boundary_meshlet_id for a triangle in a kinDS OBJ export.

Provide three corner positions (as strings to preserve decimal precision).
Tolerance per axis is half of the least significant digit you typed for that
coordinate (e.g. 7.21593 -> 5e-6 on each axis for that point's matching).

Example:
  python find_boundary_meshlet_from_obj.py combined_mesh.obj \\
    --point 1.23 4.567 7.21593 \\
    --point 1.24 4.568 7.21601 \\
    --point 1.22 4.566 7.21580
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

Vec3 = Tuple[float, float, float]
Tol3 = Tuple[float, float, float]
QueryPoint = Tuple[Vec3, Tol3]


def tolerance_from_coord_string(s: str) -> float:
    """Half of the value of the last decimal digit shown in *s*."""
    s = s.strip()
    if not s:
        return 0.5

    lower = s.lower()
    if "e" in lower:
        mantissa, exp_str = lower.split("e", 1)
        exponent = int(exp_str)
        if "." in mantissa:
            decimals = len(mantissa.split(".", 1)[1])
        else:
            decimals = 0
        return 0.5 * 10 ** (exponent - decimals)

    if "." in s:
        decimals = len(s.split(".", 1)[1])
        return 0.5 * 10 ** (-decimals)

    return 0.5


def parse_query_point(coords: Sequence[str]) -> QueryPoint:
    if len(coords) != 3:
        raise ValueError("each --point requires exactly three coordinates (x y z)")
    x_s, y_s, z_s = (c.strip() for c in coords)
    point = (float(x_s), float(y_s), float(z_s))
    tol = (
        tolerance_from_coord_string(x_s),
        tolerance_from_coord_string(y_s),
        tolerance_from_coord_string(z_s),
    )
    return point, tol


def points_match(a: Vec3, b: Vec3, tol: Tol3) -> bool:
    return (
        abs(a[0] - b[0]) <= tol[0]
        and abs(a[1] - b[1]) <= tol[1]
        and abs(a[2] - b[2]) <= tol[2]
    )


def triangle_matches_unordered(face: Tuple[Vec3, Vec3, Vec3], query: List[QueryPoint]) -> bool:
    remaining = list(query)
    for vertex in face:
        matched = False
        for i, (q_point, q_tol) in enumerate(remaining):
            if points_match(vertex, q_point, q_tol):
                remaining.pop(i)
                matched = True
                break
        if not matched:
            return False
    return True


def extract_boundary_meshlet_id(metadata: Optional[str]) -> Optional[int]:
    if not metadata:
        return None
    try:
        parsed = json.loads(metadata)
        value = parsed.get("boundary_meshlet_id")
        return int(value) if value is not None else None
    except json.JSONDecodeError:
        match = re.search(r'"boundary_meshlet_id"\s*:\s*(\d+)', metadata)
        return int(match.group(1)) if match else None


def parse_face_vertex_index(token: str) -> int:
    return int(token.split("/", 1)[0]) - 1


def triangulate_face_indices(indices: Sequence[int]) -> Iterable[List[int]]:
    if len(indices) < 3:
        return
    if len(indices) == 3:
        yield list(indices)
        return
    v0 = indices[0]
    for i in range(1, len(indices) - 1):
        yield [v0, indices[i], indices[i + 1]]


def parse_obj(path: Path) -> Tuple[List[Vec3], List[Tuple[List[int], Optional[str], int]]]:
    vertices: List[Vec3] = []
    faces: List[Tuple[List[int], Optional[str], int]] = []

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_no, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            metadata: Optional[str] = None
            if " #" in line:
                face_part, metadata = line.split(" #", 1)
                metadata = metadata.strip()
                line = face_part.strip()

            parts = line.split()
            if not parts:
                continue

            tag = parts[0]
            if tag == "v" and len(parts) >= 4:
                vertices.append((float(parts[1]), float(parts[2]), float(parts[3])))
                continue

            if tag != "f":
                continue

            indices = [parse_face_vertex_index(token) for token in parts[1:]]
            for tri in triangulate_face_indices(indices):
                faces.append((tri, metadata, line_no))

    return vertices, faces


def face_positions(vertices: Sequence[Vec3], tri: Sequence[int]) -> Tuple[Vec3, Vec3, Vec3]:
    return vertices[tri[0]], vertices[tri[1]], vertices[tri[2]]


def format_point(point: Vec3) -> str:
    return f"({point[0]:.9g}, {point[1]:.9g}, {point[2]:.9g})"


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Locate a triangle in an OBJ by three vertex positions and print boundary_meshlet_id."
    )
    parser.add_argument("obj", type=Path, help="Path to OBJ file (with kinDS face metadata comments)")
    parser.add_argument(
        "--point",
        action="append",
        nargs=3,
        metavar=("X", "Y", "Z"),
        required=True,
        help="One triangle corner as three coordinate strings (repeat three times)",
    )
    parser.add_argument(
        "--from-blender",
        action="store_true",
        help="Remap coordinates from Blender import orientation to kinDS OBJ (x,z,-y)",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=None,
        help="Override per-axis tolerance (default: inferred from coordinate strings)",
    )
    parser.add_argument(
        "--max-hits",
        type=int,
        default=20,
        help="Maximum number of matching faces to print (default: 20)",
    )
    args = parser.parse_args(argv)

    if len(args.point) != 3:
        parser.error("provide exactly three --point arguments (triangle corners)")

    obj_path: Path = args.obj
    if not obj_path.is_file():
        print(f"error: file not found: {obj_path}", file=sys.stderr)
        return 1

    query: List[QueryPoint] = [parse_query_point(p) for p in args.point]
    if args.from_blender:
        remapped: List[QueryPoint] = []
        for (point, tol) in query:
            bx, by, bz = point
            remapped.append(((bx, bz, -by), tol))
        query = remapped

    if args.tolerance is not None:
        query = [(point, (args.tolerance, args.tolerance, args.tolerance)) for point, _ in query]

    print("Query points and per-axis tolerances:")
    for i, (point, tol) in enumerate(query):
        print(f"  p{i}: {format_point(point)}  tol=({tol[0]:.6g}, {tol[1]:.6g}, {tol[2]:.6g})")

    vertices, faces = parse_obj(obj_path)
    matches: List[dict] = []

    for tri, metadata, line_no in faces:
        if any(idx < 0 or idx >= len(vertices) for idx in tri):
            continue
        positions = face_positions(vertices, tri)
        if not triangle_matches_unordered(positions, query):
            continue
        meshlet_id = extract_boundary_meshlet_id(metadata)
        matches.append(
            {
                "line": line_no,
                "vertex_indices_obj": [i + 1 for i in tri],
                "positions": positions,
                "boundary_meshlet_id": meshlet_id,
                "metadata": metadata,
            }
        )

    if not matches:
        print(f"\nNo matching triangle found in {obj_path}")
        return 2

    print(f"\nFound {len(matches)} matching triangle(s) in {obj_path}:")
    for i, hit in enumerate(matches[: args.max_hits]):
        v_idx = hit["vertex_indices_obj"]
        p0, p1, p2 = hit["positions"]
        meshlet = hit["boundary_meshlet_id"]
        meshlet_text = str(meshlet) if meshlet is not None else "(not present)"
        print(f"\n[{i + 1}] line {hit['line']}")
        print(f"    f {v_idx[0]} {v_idx[1]} {v_idx[2]}")
        print(f"    v[{v_idx[0]}] = {format_point(p0)}")
        print(f"    v[{v_idx[1]}] = {format_point(p1)}")
        print(f"    v[{v_idx[2]}] = {format_point(p2)}")
        print(f"    boundary_meshlet_id = {meshlet_text}")
        if hit["metadata"]:
            print(f"    metadata = {hit['metadata']}")

    if len(matches) > args.max_hits:
        print(f"\n... {len(matches) - args.max_hits} more match(es) omitted")

    unique_meshlets = sorted({m["boundary_meshlet_id"] for m in matches if m["boundary_meshlet_id"] is not None})
    if unique_meshlets:
        print(f"\nUnique boundary_meshlet_id values among matches: {unique_meshlets}")
    else:
        print("\nNo boundary_meshlet_id found on matching faces (yellow/regular meshlets, or metadata disabled).")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
