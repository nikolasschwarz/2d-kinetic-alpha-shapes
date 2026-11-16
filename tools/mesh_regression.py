"""Mesh regression harness for kinDS-demo.

Runs kinDS-demo on a set of StrandTree files in serial and/or parallel scheduling, fingerprints the combined export
(OBJ positions / UVs / normals plus the GPU-attribute JSON sidecar) independent of vertex and triangle order, and
records meshing time and logged errors. Parallel scheduling emits meshlets in nondeterministic order, so only an
order-independent fingerprint is comparable between runs.

  python tools/mesh_regression.py --exe build/bin/kinDS-demo.exe --trees input_examples/*.txt --out base.json
  python tools/mesh_regression.py --exe ... --trees ... --out new.json --baseline base.json

Exit code 1 when a fingerprint or error count differs from the baseline.
"""
import argparse
import concurrent.futures
import glob
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time

COMMON_ARGS = ['--mesh-cap-at-start', '--cutoff', '10', '--export-mode', 'combined', '--export-gpu-attributes']
MODES = {'serial': [], 'parallel': ['--parallelize-branches']}


def fingerprint(obj_path, json_path):
    """Order-independent fingerprint. 'hash': sorted triangles whose corners carry position, UV, normal and per-vertex
    JSON attributes (strict). 'corners' + 'area': distinct corners and total area, which ignore how a polygon was split
    into triangles (parallel scheduling may fan-triangulate a cap from a different start vertex)."""
    attrs = json.load(open(json_path, encoding='utf-8')) if os.path.exists(json_path) else {}
    face_attrs = {k: v for k, v in attrs.items() if k == 'has_neighbor'}
    vert_attrs = {k: v for k, v in attrs.items() if k != 'has_neighbor'}
    verts, uvs, normals, tris = [], [], [], []
    face_index = 0
    with open(obj_path, encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.split('#', 1)[0]
            if line.startswith('v '):
                verts.append(tuple(line.split()[1:4]))
            elif line.startswith('vt '):
                uvs.append(tuple(line.split()[1:]))
            elif line.startswith('vn '):
                normals.append(tuple(line.split()[1:4]))
            elif line.startswith('f '):
                corners = []
                for tok in line.split()[1:]:
                    parts = (tok.split('/') + ['', ''])[:3]
                    vi = int(parts[0]) - 1
                    corner = (verts[vi],
                              uvs[int(parts[1]) - 1] if parts[1] else None,
                              normals[int(parts[2]) - 1] if parts[2] else None,
                              tuple(repr(v[vi]) for _, v in sorted(vert_attrs.items())))
                    corners.append(corner)
                extra = tuple(repr(v[face_index]) for _, v in sorted(face_attrs.items()))
                face_index += 1
                for k in range(1, len(corners) - 1):
                    tri = [corners[0], corners[k], corners[k + 1]]
                    r = tri.index(min(tri))
                    tris.append((tuple(tri[r:] + tri[:r]), extra))
    tris.sort()
    area = 0.0
    for (a, b, c), _ in tris:
        p = [[float(x) for x in corner[0]] for corner in (a, b, c)]
        u = [p[1][i] - p[0][i] for i in range(3)]
        v = [p[2][i] - p[0][i] for i in range(3)]
        n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
        area += 0.5 * (n[0] ** 2 + n[1] ** 2 + n[2] ** 2) ** 0.5
    corners = sorted({corner for (tri, _) in tris for corner in tri})
    return {'triangles': len(tris), 'vertices': len(verts), 'hash': digest(tris), 'corners': digest(corners),
            'area': area}


def digest(items):
    h = hashlib.sha256()
    for item in items:
        h.update(repr(item).encode())
    return h.hexdigest()[:20]


def same_mesh(old, new, level):
    """'strict': identical triangles. 'corners': identical corners and triangle count, area within 1e-9 relative.
    'loose': triangle / vertex counts within 0.01 % and area within 1e-5 relative. Parallel scheduling already varies
    the combined mesh between runs of the same build (order-dependent meshlet post-processing; up to ~5e-6 relative
    area on large trees), so it is 'loose'."""
    if old is None or new is None:
        return old == new
    if level == 'strict':
        return old['hash'] == new['hash']
    rel_area = abs(old['area'] - new['area']) / max(abs(old['area']), 1e-30)
    if level == 'corners':
        return old['corners'] == new['corners'] and old['triangles'] == new['triangles'] and rel_area <= 1e-9
    return all(abs(old[k] - new[k]) <= max(4, 1e-4 * old[k]) for k in ('triangles', 'vertices')) and rel_area <= 1e-5


def mesh_level(mode, part):
    if mode == 'serial':
        return 'strict'
    return 'corners' if part == 'boundary' else 'loose'


def run_case(exe, tree, mode, repeat, work, keep):
    name = os.path.splitext(os.path.basename(tree))[0] + '_' + mode
    out_dir = os.path.join(work, f'{name}_r{repeat}')
    shutil.rmtree(out_dir, ignore_errors=True)
    os.makedirs(out_dir)
    cmd = [os.path.abspath(exe), '--log-file', 'run.log'] + COMMON_ARGS + MODES[mode] + ['--mesh', os.path.abspath(tree)]
    start = time.perf_counter()
    proc = subprocess.run(cmd, cwd=out_dir, capture_output=True, text=True, errors='replace')
    wall = time.perf_counter() - start
    m = re.search(r'Generated (\d+) meshlets in ([0-9.eE+-]+) s', proc.stdout)
    log_path = os.path.join(out_dir, 'run.log')
    log = open(log_path, encoding='utf-8', errors='replace').read() if os.path.exists(log_path) else ''
    result = {'tree': os.path.basename(tree), 'mode': mode, 'exit_code': proc.returncode,
              'meshlets': int(m.group(1)) if m else None, 'meshing_seconds': float(m.group(2)) if m else None,
              'wall_seconds': round(wall, 3), 'errors': len(re.findall(r'\] ERROR:', log))}
    combined = os.path.join(out_dir, 'combined_mesh.obj')
    if os.path.exists(combined):
        result['combined'] = fingerprint(combined, os.path.join(out_dir, 'combined_mesh.json'))
        result['boundary'] = fingerprint(os.path.join(out_dir, 'boundary_mesh.obj'), '')
    if not keep:
        for sub in ('raw_meshlets', 'segment_meshlets'):
            shutil.rmtree(os.path.join(out_dir, sub), ignore_errors=True)
    return name, result


def compare(baseline, results):
    failures = []
    for name, new in sorted(results.items()):
        old = baseline.get(name)
        if old is None:
            print(f'  {name}: not in baseline')
            continue
        diffs = [k for k in ('exit_code', 'errors', 'meshlets') if old.get(k) != new.get(k)]
        diffs += [k for k in ('combined', 'boundary')
                  if not same_mesh(old.get(k), new.get(k), mesh_level(new['mode'], k))]
        t_old, t_new = old.get('meshing_seconds'), new.get('meshing_seconds')
        speed = f'{t_old:.3f}s -> {t_new:.3f}s ({t_new / t_old:.2f}x)' if t_old and t_new else 'n/a'
        status = 'OK  ' if not diffs else 'DIFF'
        print(f'  {status} {name}: {speed}' + (f'  changed: {", ".join(diffs)}' if diffs else ''))
        if diffs:
            failures.append(name)
    return failures


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--exe', required=True)
    ap.add_argument('--trees', nargs='+', required=True, help='StrandTree files (globs allowed)')
    ap.add_argument('--modes', nargs='+', default=['serial', 'parallel'], choices=sorted(MODES))
    ap.add_argument('--work', default='mesh_regression_runs')
    ap.add_argument('--out', required=True, help='results JSON')
    ap.add_argument('--baseline', help='compare against this results JSON')
    ap.add_argument('--jobs', type=int, default=1, help='concurrent runs (use 1 for timing)')
    ap.add_argument('--repeat', type=int, default=1, help='runs per case; keeps the fastest meshing time')
    ap.add_argument('--keep', action='store_true', help='keep raw/segment meshlet exports')
    args = ap.parse_args()

    trees = sorted({p for pattern in args.trees for p in glob.glob(pattern)})
    cases = [(t, m, r) for t in trees for m in args.modes for r in range(args.repeat)]
    results = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for name, result in pool.map(lambda c: run_case(args.exe, *c, args.work, args.keep), cases):
            prev = results.get(name)
            if prev is not None and not all(same_mesh(prev.get(k), result.get(k), mesh_level(result['mode'], k))
                                       for k in ('combined', 'boundary')):
                print(f'WARNING: {name} differs between repeats', flush=True)
            if prev is None or (result['meshing_seconds'] or 1e9) < (prev['meshing_seconds'] or 1e9):
                results[name] = result
            print(f'{name}: exit={result["exit_code"]} errors={result["errors"]} '
                  f't={result["meshing_seconds"]} {result.get("combined", {}).get("hash")}', flush=True)
    json.dump(results, open(args.out, 'w', encoding='utf-8'), indent=1, sort_keys=True)

    if args.baseline:
        print('compare against', args.baseline)
        failures = compare(json.load(open(args.baseline, encoding='utf-8')), results)
        print(f'{len(failures)} case(s) differ')
        sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main()
