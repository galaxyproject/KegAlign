#!/usr/bin/env python
"""Group a KegAlign run's diagonal-partition splits into pair files -- one per chromosome pair.

A **pair file** is the unit Growler maps over: every anchor segment for ONE (target, query)
chromosome pair, BOTH strands, in a single gzipped segments file. `--output-type collection`
emits one pair file per element instead of tarring all the splits into `data_package.tgz`.

    build_pairs.py --commands lastz-commands.txt --segments-dir . --out pairs/
    build_pairs.py --verify --bundle /path/to/extracted --lastz /path/to/lastz --pair EH23a.chr9,EH23b.chrX
    build_pairs.py --self-test

⛔ WHY THIS IS SAFE -- the lastz 1.04.52 manual, "Segment File", is the whole specification:

    Query sequence names must appear in the same order as they do in the query file. For each
    query sequence, normally all positive strand intervals must appear before any negative
    strand intervals. Sequence names for the target may appear in any order, and are only
    meaningful if the `multiple` action is used.

A pair file holds exactly ONE query sequence, so the first rule is trivially satisfied. The second is
why plus splits are written before minus. The third is why target order inside a pair file is free --
and why the 35-of-5177 splits that span more than one target chromosome can be cut apart by
target with no sorting: their target names arrive in contiguous runs.

⚠ MEASURED, NOT ASSUMED. Over all 5,177 splits of a real bundle: 0 span more than one QUERY
sequence, 0 have a query name recurring after another, 35 span more than one TARGET. An earlier
150-file sample found zero exceptions and was reported as "zero" -- at a 0.68% rate a sample
that size expects one. A structural property a design rests on gets the exhaustive pass.

⚠ THE ELEMENTS ARE GZIPPED, AND THAT IS NOT OPTIONAL. Segments are 97.1% of a bundle and the
tarball compresses the whole thing ~3.6x (14.58 GB -> 4.10 GB measured). Emitting pair files as plain
datasets would roughly triple the object store. But lastz CANNOT read a gzipped segments file --
`--segments=x.gz` gives `FAILURE: bad field (x.gz: line 1, ...)` and zero output -- so the
consumer must decompress its own pair file before invoking lastz. That is ~140 MB per job, against
the current runner inflating all 14.58 GB into the job directory before anything runs.
"""
from __future__ import annotations

import argparse
import collections
import gzip
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile
import typing

#: Separates target and query in a pair file's element identifier. Two underscores, because a single
#: one occurs inside real sequence names; Galaxy element identifiers are matched by regex in
#: `__APPLY_RULES__`, and a separator that appears in the payload makes that regex ambiguous.
PAIR_SEP = "__"

SEGMENT_COLUMNS = 8


class PairKey(typing.NamedTuple):
    target: str
    query: str

    @property
    def identifier(self) -> str:
        return f"{self.target}{PAIR_SEP}{self.query}"


def arg_value(args: list[str], prefix: str) -> str | None:
    for arg in args:
        if arg.startswith(prefix):
            return arg[len(prefix):]
    return None


def strand_of(command_args: list[str]) -> str:
    """The strand a command searches, as the '+'/'-' that appears in its segment lines.

    ⚠ `--strand` is ABSENT from a command only when the caller wants both; lastz's default is
    `both`. A pair file is written plus-then-minus precisely so that the merged file is legal under
    the default, and the consumer can then drop `--strand` altogether.
    """
    value = arg_value(command_args, "--strand=")
    if value == "minus":
        return "-"
    if value == "plus":
        return "+"
    return "both"


def split_by_target(lines: list[str]) -> list[tuple[str, list[str]]]:
    """Cut a segment file into contiguous runs sharing a target name.

    Returns [(target, lines), ...] preserving input order. For the 5,142 splits that hold a
    single target this is one run; for the 35 that hold two or three it is that many, and no
    sorting is needed because the runs are already contiguous.
    """
    runs: list[tuple[str, list[str]]] = []
    for line in lines:
        fields = line.split("\t")
        if len(fields) < SEGMENT_COLUMNS - 1:
            raise ValueError(f"segment line has {len(fields)} fields: {line!r}")
        target = fields[0]
        if runs and runs[-1][0] == target:
            runs[-1][1].append(line)
        else:
            runs.append((target, [line]))
    return runs


def query_of(lines: list[str]) -> str:
    """The single query name in a segment file.

    ⛔ Raises if there is more than one. A pair file spanning two query sequences would violate the
    manual's ordering rule, and this is the one invariant that cannot be repaired by reordering
    here -- it would need a sort against the 2bit's sequence order, which this script does not
    have. Failing loudly is correct; silently emitting an illegal pair file is not.
    """
    names = {line.split("\t")[3] for line in lines}
    if len(names) != 1:
        raise ValueError(f"split spans {len(names)} query sequences: {sorted(names)}")
    return names.pop()


def assign(commands: list[dict], read_segments: typing.Callable[[str], list[str]]) -> dict[PairKey, dict[str, list[str]]]:
    """Group every segment line of every command into pair_files, keyed (target, query).

    The returned mapping is pair file -> {"+": [...], "-": [...]}, each list in the order the
    commands were given. Writing is a separate step so this can be unit-tested without a
    filesystem, and so the ordering rule is visible in one place.
    """
    pair_files: dict[PairKey, dict[str, list[str]]] = collections.defaultdict(lambda: {"+": [], "-": []})
    for command in commands:
        segments = arg_value(command.get("args", []), "--segments=")
        if segments is None:
            continue
        lines = read_segments(segments)
        if not lines:
            continue
        query = query_of(lines)
        strand = strand_of(command["args"])
        for target, run in split_by_target(lines):
            if strand == "both":
                # the file itself carries strand per line; trust it rather than the flag
                for line in run:
                    pair_files[PairKey(target, query)][line.split("\t")[6]].append(line)
            else:
                pair_files[PairKey(target, query)][strand].extend(run)
    return pair_files


def pair_lines(strands: dict[str, list[str]]) -> list[str]:
    """Plus strand first, then minus -- the manual's rule, applied in exactly one place."""
    return strands["+"] + strands["-"]


def write_pairs(pair_files: dict[PairKey, dict[str, list[str]]], out_dir: pathlib.Path, compresslevel: int = 6) -> list[tuple[str, int, int]]:
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest: list[tuple[str, int, int]] = []
    for key in sorted(pair_files, key=lambda k: k.identifier):
        lines = pair_lines(pair_files[key])
        path = out_dir / f"{key.identifier}.segments.gz"
        payload = "".join(line + "\n" for line in lines).encode()
        with gzip.open(path, "wb", compresslevel=compresslevel) as fh:
            fh.write(payload)
        manifest.append((key.identifier, len(lines), path.stat().st_size))
    return manifest


# --------------------------------------------------------------------------------------- verify

MAF_BLOCK_START = re.compile(r"^a score=")


def maf_blocks(path: pathlib.Path) -> collections.Counter:
    """MAF blocks as a multiset, so two runs can be compared regardless of block order."""
    blocks: collections.Counter = collections.Counter()
    current: list[str] = []
    for line in path.read_text().splitlines():
        if MAF_BLOCK_START.match(line):
            if current:
                blocks[tuple(current)] += 1
            current = [line]
        elif line.strip():
            current.append(line)
    if current:
        blocks[tuple(current)] += 1
    return blocks


def run_lastz(lastz: str, target: str, query: str, segments: pathlib.Path, out: pathlib.Path,
              extra: list[str], workdir: pathlib.Path) -> float:
    import time
    argv = [lastz, target, query, "--allocate:traceback=1.99G", "--format=maf-",
            f"--segments={segments}", f"--output={out}", *extra]
    begin = time.perf_counter()
    proc = subprocess.run(argv, cwd=workdir, capture_output=True, text=True)
    elapsed = time.perf_counter() - begin
    if proc.returncode != 0 and proc.stderr.strip():
        sys.exit(f"lastz failed ({proc.returncode}): {proc.stderr[:600]}")
    return elapsed


def verify(args: argparse.Namespace) -> int:
    """Run lastz once on a pair file, and once per constituent split, and compare the alignments.

    ⛔ THIS IS THE CLAIM THE WHOLE DESIGN RESTS ON, so it is tested against real data and real
    lastz rather than argued from the manual. A pair file is only legitimate if merging splits changes
    nothing about what comes out.
    """
    import json

    bundle = pathlib.Path(args.bundle).resolve()
    workdir = bundle / "galaxy" / "files"
    commands = [json.loads(line) for line in (bundle / "galaxy" / "commands.json").read_text().splitlines() if line.strip()]

    def read_segments(name: str) -> list[str]:
        return (workdir / name).read_text().splitlines()

    want_target, want_query = (args.pair.split(",") + [""])[:2]
    selected = []
    for command in commands:
        segments = arg_value(command["args"], "--segments=")
        if segments is None:
            continue
        head = (workdir / segments).open().readline().split("\t")
        if head[0] == want_target and head[3] == want_query:
            selected.append(command)
    if not selected:
        sys.exit(f"no commands for pair {args.pair}")
    selected.sort(key=lambda c: arg_value(c["args"], "--output=") or "")
    if args.smallest:
        # ⚠ Selecting the SMALLEST splits, not the first N. A whole chromosome pair is 86-256 MB
        # of segments and takes hours in one lastz; the merge property is the same code path at
        # any size, so a size-bounded selection answers it in minutes. The selection is still
        # returned in command order, so the pair file is built exactly as it would be in production.
        by_size = sorted(selected, key=lambda c: (workdir / arg_value(c["args"], "--segments=")).stat().st_size)
        keep = {id(c) for c in by_size[: args.smallest]}
        selected = [c for c in selected if id(c) in keep]
    elif args.max_splits:
        selected = selected[: args.max_splits]
    print(f"pair {want_target} x {want_query}: {len(selected)} split(s)", file=sys.stderr)

    pair_files = assign(selected, read_segments)
    if len(pair_files) != 1:
        print(f"  note: {len(pair_files)} pair files in this selection: {[k.identifier for k in pair_files]}", file=sys.stderr)

    tmp = pathlib.Path(tempfile.mkdtemp(prefix="pair_verify_", dir=args.tmpdir))
    manifest = write_pairs(pair_files, tmp / "pairs")
    for identifier, n_lines, n_bytes in manifest:
        print(f"  pair file {identifier}: {n_lines:,} segments, {n_bytes:,} bytes gzipped", file=sys.stderr)

    target_spec = arg_value(selected[0]["args"], "--target=")
    query_spec = arg_value(selected[0]["args"], "--query=")
    shared = [a for a in selected[0]["args"]
              if not a.startswith(("--target=", "--query=", "--segments=", "--output=", "--strand=", "--format="))]

    # --- arm A: one lastz per pair file, no --strand (default is both; measured identical) ---
    per_pair: collections.Counter = collections.Counter()
    pair_seconds = 0.0
    for identifier, _, _ in manifest:
        plain = tmp / f"{identifier}.segments"
        with gzip.open(tmp / "pairs" / f"{identifier}.segments.gz", "rb") as src, plain.open("wb") as dst:
            shutil.copyfileobj(src, dst)          # the consumer's gunzip step, in miniature
        out = tmp / f"{identifier}.pairfile.maf-"
        pair_seconds += run_lastz(args.lastz, target_spec, query_spec, plain, out, shared, workdir)
        per_pair += maf_blocks(out)

    # --- arm B: the original splits, each with its own --strand ---
    per_split: collections.Counter = collections.Counter()
    split_seconds = 0.0
    for command in selected:
        segments = workdir / arg_value(command["args"], "--segments=")
        strand = [a for a in command["args"] if a.startswith("--strand=")]
        out = tmp / (pathlib.Path(arg_value(command["args"], "--output=")).name + ".split")
        split_seconds += run_lastz(args.lastz, target_spec, query_spec, segments, out, shared + strand, workdir)
        per_split += maf_blocks(out)

    print(f"\n{'':<22}{'blocks':>10}{'seconds':>10}", file=sys.stderr)
    print(f"{'pair file (1 lastz)':<22}{sum(per_pair.values()):>10,}{pair_seconds:>10.1f}", file=sys.stderr)
    print(f"{'splits (%d lastz)' % len(selected):<22}{sum(per_split.values()):>10,}{split_seconds:>10.1f}", file=sys.stderr)

    # ⛔ COMPARE BY LOCUS, NOT BY BLOCK TEXT. lastz's own manual warns that the same alignment
    # may come back with "minor variations such as shifting of equally-scoring gap placements".
    # Two runs can therefore emit the same alignment -- same score, same target and query
    # coordinates, same length -- as different STRINGS. Comparing full block text calls that a
    # difference, which on a 238,076-segment pair file reported three alignments as "missing from the
    # pair file" when every one of them was present. The locus is the thing that has to match.
    locus = lambda b: (b[0], b[1].split()[2], b[2].split()[2])
    pair_loci = collections.Counter(locus(b) for b in per_pair.elements())
    split_loci = collections.Counter(locus(b) for b in per_split.elements())

    lost = set(split_loci) - set(pair_loci)
    gained = set(pair_loci) - set(split_loci)
    # ▶ AND THE DUPLICATES ARE A FINDING ABOUT THE SPLIT PATH, NOT ABOUT THE KEG. When anchors
    # for one alignment fall in two different splits, each process finds it independently and
    # both report it. One lastz sees them together and reports it once, so a pair file DEDUPLICATES
    # what the current production path emits twice.
    dup_splits = {k: v for k, v in split_loci.items() if v > 1}
    dup_pair = {k: v for k, v in pair_loci.items() if v > 1}

    print(f"\n{'':<22}{'blocks':>10}{'loci':>10}", file=sys.stderr)
    print(f"{'pair file':<22}{sum(pair_loci.values()):>10,}{len(pair_loci):>10,}", file=sys.stderr)
    print(f"{'splits':<22}{sum(split_loci.values()):>10,}{len(split_loci):>10,}", file=sys.stderr)
    if dup_splits or dup_pair:
        print(f"\nloci reported more than once: splits {len(dup_splits)} "
              f"(+{sum(v - 1 for v in dup_splits.values())} blocks), pair file {len(dup_pair)}",
              file=sys.stderr)

    if not lost and not gained:
        extra = sum(v - 1 for v in dup_splits.values())
        print(f"\nok - identical alignment loci; the pair file de-duplicates {extra} block(s) the "
              f"splits report twice", file=sys.stderr)
        return 0
    print(f"\nnot ok - {len(lost)} locus/loci only in the splits, {len(gained)} only in the pair file",
          file=sys.stderr)
    for loc in list(lost)[:3]:
        print("  split-only locus:", loc, file=sys.stderr)
    for loc in list(gained)[:3]:
        print("  pair file-only locus  :", loc, file=sys.stderr)
    return 1


def self_test() -> int:
    failures: list[str] = []

    def check(label: str, got: typing.Any, want: typing.Any) -> None:
        if got != want:
            print(f"not ok - {label}\n     got: {got}\nexpected: {want}")
            failures.append(label)
        else:
            print(f"ok - {label}")

    def seg(t: str, q: str, s: str, start: int) -> str:
        return f"{t}\t{start}\t{start+9}\t{q}\t{start}\t{start+9}\t{s}\t3000"

    # a split holding one target, and one holding three contiguous runs
    single = [seg("T1", "Q1", "+", i) for i in range(3)]
    multi = ([seg("T1", "Q1", "+", 0)] + [seg("T2", "Q1", "+", 1)] * 2 + [seg("T3", "Q1", "+", 2)])
    check("single-target split is one run", [t for t, _ in split_by_target(single)], ["T1"])
    check("multi-target split cuts into contiguous runs", [t for t, _ in split_by_target(multi)], ["T1", "T2", "T3"])
    check("cutting loses no lines", sum(len(r) for _, r in split_by_target(multi)), len(multi))

    check("query_of finds the single query", query_of(single), "Q1")
    try:
        query_of(single + [seg("T1", "Q2", "+", 9)])
        check("two queries must raise", "no raise", "ValueError")
    except ValueError:
        check("two queries must raise", "ValueError", "ValueError")

    check("--strand=minus reads as '-'", strand_of(["--strand=minus"]), "-")
    check("absent --strand means both", strand_of(["--ydrop=1"]), "both")

    # ▶ THE ORDERING RULE. Build a pair file from a minus command given BEFORE a plus command and
    # assert the written order is still plus-then-minus -- this is the invariant the lastz
    # manual imposes, and the one a naive concatenation in command order would break.
    store = {"m.segments": [seg("T1", "Q1", "-", 5)], "p.segments": [seg("T1", "Q1", "+", 1)]}
    cmds = [{"args": ["--segments=m.segments", "--strand=minus"]},
            {"args": ["--segments=p.segments", "--strand=plus"]}]
    pair_files = assign(cmds, lambda n: store[n])
    check("one pair file from one chromosome pair", [k.identifier for k in pair_files], ["T1__Q1"])
    written = pair_lines(pair_files[PairKey("T1", "Q1")])
    check("plus is written before minus, whatever the command order",
          [ln.split("\t")[6] for ln in written], ["+", "-"])

    # a multi-target split lands in as many pair files as it has targets, losing nothing
    pair_files2 = assign([{"args": ["--segments=x", "--strand=plus"]}], lambda n: multi)
    check("multi-target split fans into 3 pair files", sorted(k.identifier for k in pair_files2),
          ["T1__Q1", "T2__Q1", "T3__Q1"])
    check("  and every line survives the fan-out",
          sum(len(pair_lines(v)) for v in pair_files2.values()), len(multi))

    # ▶ THE DIRECTORY-SCAN PATH: no --strand in the synthesised command, so the strand has to
    # come from column 7 of each line. A pair file built this way must match one built from commands.
    mixed = [seg("T1", "Q1", "-", 5), seg("T1", "Q1", "+", 1)]
    from_file = assign([{"args": ["--segments=x"]}], lambda n: mixed)
    check("strand read from column 7 when --strand is absent",
          [ln.split("\t")[6] for ln in pair_lines(from_file[PairKey("T1", "Q1")])], ["+", "-"])
    from_cmds = assign([{"args": ["--segments=m", "--strand=minus"]},
                        {"args": ["--segments=p", "--strand=plus"]}],
                       lambda n: [mixed[0]] if n == "m" else [mixed[1]])
    check("  and agrees with the command-driven path",
          pair_lines(from_file[PairKey("T1", "Q1")]), pair_lines(from_cmds[PairKey("T1", "Q1")]))

    # ⚠ the separator must not occur in the payload it separates
    check("separator is two underscores", PAIR_SEP, "__")
    check("a single-underscore name still round-trips",
          PairKey("chr_1", "chr_2").identifier.split(PAIR_SEP), ["chr_1", "chr_2"])

    if failures:
        print(f"\n{len(failures)} test(s) failed")
        return 1
    print("\nall tests passed")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--commands", help="commands.json to group (optional; segments carry "
                                          "target, query and strand themselves)")
    parser.add_argument("--segments-dir", default=".", help="directory holding the .segments files")
    parser.add_argument("--out", default="pairs", help="directory to write pair files into")
    parser.add_argument("--compresslevel", type=int, default=6)
    parser.add_argument("--verify", action="store_true", help="compare a pair file against its splits with real lastz")
    parser.add_argument("--bundle", help="extracted bundle, for --verify")
    parser.add_argument("--pair", help="TARGET,QUERY to verify, e.g. EH23a.chr9,EH23b.chrX")
    parser.add_argument("--max-splits", type=int, default=0, help="first N splits in --verify (0 = all)")
    parser.add_argument("--smallest", type=int, default=0, help="the N SMALLEST splits, for a quick --verify")
    parser.add_argument("--lastz", default="lastz")
    parser.add_argument("--tmpdir", default=None)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    if args.self_test:
        return self_test()
    if args.verify:
        if not (args.bundle and args.pair):
            parser.error("--verify needs --bundle and --pair")
        return verify(args)

    seg_dir = pathlib.Path(args.segments_dir)
    if args.commands:
        import json
        text = pathlib.Path(args.commands).read_text().splitlines()
        commands = [json.loads(line) for line in text if line.strip()]
    else:
        # ▶ NO COMMAND MANIFEST IS NEEDED. Everything the grouping depends on -- target, query
        # and strand -- is in the segment lines themselves (columns 1, 4 and 7). Reading
        # `commands.json` would mean either parsing `lastz-commands.txt` with bashlex a second
        # time or running after package_output.py, and this repository already carries FOUR
        # drifted copies of two scripts. A second parser is the last thing it needs.
        commands = [{"args": [f"--segments={p.name}"]} for p in sorted(seg_dir.glob("*.segments"))]
        if not commands:
            sys.exit(f"ERROR: no *.segments files under {seg_dir}")
    manifest = write_pairs(assign(commands, lambda n: (seg_dir / n).read_text().splitlines()),
                          pathlib.Path(args.out), args.compresslevel)
    total_lines = sum(n for _, n, _ in manifest)
    total_bytes = sum(b for _, _, b in manifest)
    print(f"{len(manifest)} pair_files, {total_lines:,} segments, {total_bytes/1e9:.2f} GB gzipped")
    return 0


if __name__ == "__main__":
    sys.exit(main())
