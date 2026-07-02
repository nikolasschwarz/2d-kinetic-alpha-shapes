#!/usr/bin/env python3
"""Extract input branches 1 and 2 at global sections 15-16 from test_oak_2.txt."""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_ints(line: str) -> list[int]:
    return [int(x) for x in line.split()]


def read_nonempty_lines(path: Path) -> list[str]:
    lines: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        stripped = raw.split("#", 1)[0].strip()
        if stripped:
            lines.append(stripped)
    return lines


def branch_at(per_strand: list[int], height: int) -> int:
    idx = min(height, len(per_strand) - 1)
    return per_strand[idx]


def support_at(per_strand: list[tuple[float, float]], height: int) -> tuple[float, float]:
    if height >= len(per_strand):
        raise IndexError(f"support height {height} out of range (len={len(per_strand)})")
    return per_strand[height]


def find_section(lines: list[str], name: str) -> int:
    header = f"[Section] {name}"
    for index, line in enumerate(lines):
        if line == header or line == name:
            return index + 1
    raise ValueError(f"section not found: {name}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("input_examples/test_oak_2.txt"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tests/fixtures/test_oak_2_branches_1_2_sections_15_16.txt"),
    )
    parser.add_argument("--section-start", type=int, default=15)
    parser.add_argument("--section-end", type=int, default=16)
    parser.add_argument("--branches", type=str, default="1,2")
    args = parser.parse_args()

    section_indices = list(range(args.section_start, args.section_end + 1))
    wanted_branches = {int(x.strip()) for x in args.branches.split(",") if x.strip()}

    lines = read_nonempty_lines(args.input)
    assert lines[0] == "StrandTree 1"
    tree_height = int(lines[1])
    strand_count = int(lines[2])
    _ = tree_height  # validated against sliced output

    support_start = find_section(lines, "support_points")
    subdiv_start = find_section(lines, "subdivisions_by_strand")
    physics_start = find_section(lines, "physics_strand_to_segment_indices")
    transforms_start = find_section(lines, "transforms_by_height_and_branch")
    branch_indices_start = find_section(lines, "branch_indices")
    strands_by_branch_start = find_section(lines, "strands_by_branch_id")

    support_points: list[list[tuple[float, float]]] = []
    cursor = support_start
    for _ in range(strand_count):
        count = int(lines[cursor])
        cursor += 1
        points = []
        for _ in range(count):
            x, y = map(float, lines[cursor].split())
            points.append((x, y))
            cursor += 1
        support_points.append(points)

    branch_indices: list[list[int]] = []
    cursor = branch_indices_start
    for _ in range(strand_count):
        values = parse_ints(lines[cursor])
        branch_indices.append(values[1:])
        cursor += 1

    selected_old_strands = [
        strand_id
        for strand_id in range(strand_count)
        if all(height < len(support_points[strand_id]) for height in section_indices)
        and any(
            branch_at(branch_indices[strand_id], height) in wanted_branches for height in section_indices
        )
    ]
    old_to_new = {old: new for new, old in enumerate(selected_old_strands)}

    sliced_support = [
        [support_at(support_points[old], height) for height in section_indices] for old in selected_old_strands
    ]

    sliced_branch_indices = [
        [branch_at(branch_indices[old], height) for height in section_indices] for old in selected_old_strands
    ]

    # Parse transforms for the requested global heights.
    cursor = transforms_start
    global_transform_height_count = int(lines[cursor])
    cursor += 1
    global_transforms: list[list[list[float]]] = []
    for _ in range(global_transform_height_count):
        branch_count = int(lines[cursor])
        cursor += 1
        by_branch = []
        for _ in range(branch_count):
            matrix = list(map(float, lines[cursor].split()))
            cursor += 1
            by_branch.append(matrix)
        global_transforms.append(by_branch)

    max_branch_id = max(
        branch_at(branch_indices[old], height)
        for old in selected_old_strands
        for height in section_indices
    )
    sliced_transforms = []
    for height in section_indices:
        source = global_transforms[height]
        by_branch = []
        for branch_id in range(max_branch_id + 1):
            if branch_id < len(source):
                by_branch.append(source[branch_id])
            else:
                by_branch.append([1.0 if row == col else 0.0 for col in range(4) for row in range(4)])
        sliced_transforms.append(by_branch)

    # Parse global strands_by_branch_id and slice membership for the selected strands.
    cursor = strands_by_branch_start
    global_height_count = int(lines[cursor])
    cursor += 1
    global_membership: list[list[list[int]]] = []
    for _ in range(global_height_count):
        branch_count = int(lines[cursor])
        cursor += 1
        by_branch = []
        for _ in range(branch_count):
            values = parse_ints(lines[cursor])
            by_branch.append(values[1:])
            cursor += 1
        global_membership.append(by_branch)

    sliced_membership: list[list[list[int]]] = []
    for height in section_indices:
        by_branch: list[list[int]] = []
        for branch_id in range(max_branch_id + 1):
            old_members = (
                global_membership[height][branch_id]
                if branch_id < len(global_membership[height])
                else []
            )
            remapped = [old_to_new[old] for old in old_members if old in old_to_new]
            by_branch.append(remapped)
        sliced_membership.append(by_branch)

    local_height = len(section_indices) - 1
    out_lines = [
        "# Extracted from test_oak_2: input branches 1 and 2 at global sections "
        f"{args.section_start}-{args.section_end}.",
        "StrandTree 1",
        str(local_height),
        str(len(selected_old_strands)),
        "[Section] support_points",
    ]
    for points in sliced_support:
        out_lines.append(str(len(points)))
        for x, y in points:
            out_lines.append(f"{x} {y}")

    out_lines.append("[Section] subdivisions_by_strand")
    for _ in selected_old_strands:
        out_lines.append("0")

    out_lines.append("[Section] physics_strand_to_segment_indices")
    for _ in selected_old_strands:
        out_lines.append("0")

    out_lines.append("[Section] transforms_by_height_and_branch")
    out_lines.append(str(len(sliced_transforms)))
    for by_branch in sliced_transforms:
        out_lines.append(str(len(by_branch)))
        for matrix in by_branch:
            out_lines.append(" ".join(f"{value:.17g}" for value in matrix))

    out_lines.append("[Section] branch_indices")
    for per_strand in sliced_branch_indices:
        out_lines.append(str(len(per_strand)) + " " + " ".join(map(str, per_strand)))

    out_lines.append("[Section] strands_by_branch_id")
    out_lines.append(str(len(sliced_membership)))
    for by_branch in sliced_membership:
        out_lines.append(str(len(by_branch)))
        for members in by_branch:
            out_lines.append(str(len(members)) + " " + " ".join(map(str, members)))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output} with {len(selected_old_strands)} strands at local height {local_height}")
    print("Old strand ids:", selected_old_strands)


if __name__ == "__main__":
    main()
