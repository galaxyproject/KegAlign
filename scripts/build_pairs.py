#!/usr/bin/env python
"""Group a KegAlign run's diagonal-partition splits into pair files -- one per chromosome pair.

A **pair file** is the unit Growler maps over: every anchor segment for ONE (target, query)
chromosome pair, BOTH strands, in a single gzipped segments file. `--output-type collection`
emits one pair file per element instead of tarring all the splits into `data_package.tgz`.

    build_pairs.py --commands lastz-commands.txt --segments-dir . --out pairs/
    build_pairs.py --commands lastz-commands.txt --out pairs/ --batch-lines 2000000 --query-2bit query.2bit
    build_pairs.py --verify --bundle /path/to/extracted --lastz /path/to/lastz --pair EH23a.chr9,EH23b.chrX
    build_pairs.py --verify-batches --bundle /path/to/extracted --smallest 24
    build_pairs.py --self-test

⛔ WHY THIS IS SAFE -- the lastz 1.04.52 manual, "Segment File", is the whole specification:

    Query sequence names must appear in the same order as they do in the query file. For each
    query sequence, normally all positive strand intervals must appear before any negative
    strand intervals. Sequence names for the target may appear in any order, and are only
    meaningful if the `multiple` action is used.

A pair file holds exactly ONE query sequence, so the first rule is trivially satisfied. The second is
why plus splits are written before minus. The third is why target order inside a pair file is free --
and why the 35-of-5177 splits that span more than one target chromosome can be cut apart by target
with no sorting.

⚠ CONTIGUITY IS NOT UNIVERSAL, AND THIS DOCSTRING CLAIMED IT WAS. Measured over all 35: target
names arrive in contiguous runs in 33 of them. The two exceptions --
`tmp46.block0.r454826725.{plus,minus}.segments`, the only 2 of 5,177 files with no `.splitN` suffix,
i.e. the diagonal partition's unsplit remainder -- interleave FOUR target names across 22,464 and
22,445 runs rather than four blocks.

Correctness does not depend on contiguity and never did: `split_by_target` compares each line to the
previous one, so an interleaved file yields one run per line and every line still routes to the
right writer. What was wrong is the stated reason, which is the more dangerous half -- a reader who
believed it would think a cheaper grouping was safe.

▶ `--batch-lines` PACKS SEVERAL QUERIES INTO ONE FILE, and the same three rules are what make it
legal. One file per chromosome pair makes the element count scale with the QUERY'S SEQUENCE COUNT,
which is fine for ten chromosomes and ruinous for a contig-level assembly: a 300,000-contig query
would emit a Galaxy element per contig per target chromosome. Batching emits queries in query-file
order (rule 1), each query's plus lines before its minus lines (rule 2), and stops splitting by
target (rule 3) -- so the count is chosen rather than inherited. `read_2bit_order` supplies the
order. Nothing about a single-query pair file changes; `--batch-lines 0` is the default and its
output is unchanged.

⚠ WHICH OF THOSE RULES LASTZ ACTUALLY ENFORCES, measured on lastz 1.04.52 against a 6-target x
25-query fixture carrying 256 anchors on both strands, each perturbation run through
`growler_lastz`'s own command line (2bit inputs, `[multiple]` target, gzipped pair file, no
`--strand`):

    query order correct          identical MAF blocks       0 errors
    query order REVERSED         no output at all           FAILURE: extra segments in file
                                                            "(for this usage segments must appear
                                                             in the same order as the query file...)"
    minus written before plus    identical MAF blocks       0 errors
    strands fully interleaved    identical MAF blocks       0 errors

So the QUERY-ORDER rule is load-bearing and fails LOUDLY -- an out-of-order batch produces no
alignments and a non-zero exit, not a quiet subset. The STRAND-ORDER rule was not enforced at all
here, which is what the manual's "normally" is doing. ⛔ IT IS STILL OBEYED: it costs nothing to
emit, one fixture on one build is not grounds to overturn a documented contract, and the whole
reason this file writes plus-then-minus is so the consumer can drop `--strand` entirely.

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
import hashlib
import json
import pathlib
import re
import shutil
import struct
import subprocess
import sys
import tempfile
import time
import typing

#: Separates target and query in a pair file's element identifier. Two underscores, because a single
#: one occurs inside real sequence names; Galaxy element identifiers are matched by regex in
#: `__APPLY_RULES__`, and a separator that appears in the payload makes that regex ambiguous.
PAIR_SEP = "__"

SEGMENT_COLUMNS = 8

#: gzip level for the pair files.
#:
#: ⚠ LEVEL 1, DELIBERATELY. Measured on real Cannabis intermediates: level 6 writes at 22-26 MB/s,
#: level 1 at 141-151 MB/s, for about +24% bytes. This write is serial and sits inside the KegAlign
#: job, so the ~6x throughput is job wall clock; the bytes are object store, which is not the
#: constraint.
#:
#: ▶ IT MOVES A PUBLISHED FIGURE. At level 6 the whole 4.10 GB bundle came to 3.74 GB of pair files,
#: 4.16 GB with the two 2bits -- a 1.5% wash against the tarball. At level 1 expect ~4.6 GB of pair
#: files, ~5.0 GB with the 2bits: about 22% more object store than the tarball it replaces, still
#: against a tarball that must then be inflated whole into a job directory.
COMPRESSLEVEL = 1


class PairKey(typing.NamedTuple):
    target: str
    query: str

    @property
    def identifier(self) -> str:
        return f"{self.target}{PAIR_SEP}{self.query}"


def arg_value(args: list[str], prefix: str) -> str | None:
    for arg in args:
        if arg.startswith(prefix):
            return arg[len(prefix) :]
    return None


def require_arg(command_args: list[str], prefix: str) -> str:
    """`arg_value`, but for the arguments a command cannot legally be missing.

    Returning `str | None` everywhere pushes an `is None` check onto every caller and mypy is
    right to insist on it. Where the absence would be a malformed command rather than an option
    left unset, say so here once.
    """
    value = arg_value(command_args, prefix)
    if value is None:
        raise SystemExit(f"ERROR: command has no {prefix}...: {command_args}")
    return value


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

    Returns [(target, lines), ...] preserving input order. For the 5,142 splits that hold a single
    target this is one run. For the other 35 it is one run per CONTIGUOUS BLOCK, which is not the
    same as one run per target name.

    ⚠ MEASURED OVER ALL 35, because the earlier wording ("two or three") was a quantifier asserted
    over a population and checked on the ones looked at: 23 hold two target names, 8 hold three and
    **4 hold four**. In 33 the names arrive contiguously, so runs == names; in the remaining 2 they
    interleave, and those yield ~22,450 single-line runs apiece. Comparing against the previous line
    is what makes both cases correct without a sort.
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

    ⛔ Raises if there is more than one, and that is about the INPUT split, not the output file.
    KegAlign emits one query per split (measured: 0 of 5,177 span two), so a split that holds more
    is a malformed bundle rather than a case to accommodate -- there is no evidence about what
    order its lines would be in, and inventing one is how a silent wrong answer starts.

    ⚠ THIS IS NO LONGER WHAT LIMITS A PAIR FILE TO ONE QUERY. It used to be: the reason given here
    was that repairing the order "would need a sort against the 2bit's sequence order, which this
    script does not have". It does now -- `read_2bit_order` -- and `batch_pairs` uses it to put
    many queries in one file. This check stays because a multi-query SPLIT is still a bug upstream.
    """
    names = {line.split("\t")[3] for line in lines}
    if len(names) != 1:
        raise ValueError(f"split spans {len(names)} query sequences: {sorted(names)}")
    return names.pop()


def read_2bit_order(path: pathlib.Path) -> list[str]:
    """Sequence names in the order the `.2bit` stores them -- the order lastz will read them in.

    ⛔ THIS IS PARSED, NOT SHELLED OUT TO `twoBitInfo`. The order is a CORRECTNESS input to
    batching (see `batch_pairs`), and a batch written against the wrong order is silently wrong
    rather than loudly broken. Adding a runtime dependency on a UCSC binary to obtain it would
    also put the one number that matters behind a tool the conda package does not install.

    The format is fixed and small: a 16-byte header (magic, version, sequenceCount, reserved)
    then `sequenceCount` index records of `nameSize:uint8, name, offset:uint32`. Byte order is
    whichever way the magic reads, which is how `twoBitInfo` itself decides.
    """
    blob = path.read_bytes()
    if len(blob) < 16:
        raise ValueError(f"{path}: too short to be a 2bit file ({len(blob)} bytes)")
    for endian in ("<", ">"):
        magic, version, count, _reserved = struct.unpack_from(endian + "4I", blob, 0)
        if magic == 0x1A412743:
            break
    else:
        raise ValueError(f"{path}: not a 2bit file (magic {blob[:4]!r})")
    # ⚠ version 1 moves the index offsets to 64 bits. Nothing in this project writes one, and
    # reading it as 32-bit would yield plausible-looking garbage, so refuse rather than guess.
    if version != 0:
        raise ValueError(f"{path}: 2bit version {version} is not supported (only 0)")
    names: list[str] = []
    pos = 16
    for _ in range(count):
        size = blob[pos]
        pos += 1
        names.append(blob[pos : pos + size].decode())
        pos += size + 4
    return names


def batch_pairs(
    pair_files: dict[PairKey, dict[str, list[str]]], query_order: list[str], max_lines: int
) -> list[tuple[str, list[str]]]:
    """Pack whole query sequences into as few legal segment files as possible.

    Returns `[(identifier, lines), ...]`. One query is never split across two batches, so a query
    whose own segments exceed `max_lines` simply becomes an oversized batch of one.

    ⚠ THIS IS THE BUFFERED ORACLE, NOT THE PRODUCTION WRITER -- `stream_batches` is what runs, for
    the same reason `stream_pairs` replaced `write_pairs`. It is kept because `--verify-batches`
    and the self-test compare the two, which is what makes the streaming one checkable at all.

    ▶ WHY THIS IS LEGAL, and why it is the SAME specification the one-query form rests on. The
    lastz 1.04.52 manual, "Segment File", imposes exactly three things:

        Query sequence names must appear in the same order as they do in the query file. For each
        query sequence, normally all positive strand intervals must appear before any negative
        strand intervals. Sequence names for the target may appear in any order, and are only
        meaningful if the `multiple` action is used.

    A batch emits queries in `query_order` (rule 1), and each query's plus lines before its minus
    lines (rule 2). Rule 3 is why the per-target split is DROPPED here rather than preserved: it
    was never required, and `growler_lastz` already passes `target.2bit[multiple]` with the whole
    2bit and no `subset=`. Splitting by target multiplied the element count by the target's
    chromosome count for no gain.

    ⛔ EVERY QUERY MUST APPEAR IN `query_order`. A query the order does not name cannot be placed,
    and guessing its position is precisely the silent-wrongness this function exists to avoid.
    """
    by_query: dict[str, dict[str, list[str]]] = collections.defaultdict(lambda: {"+": [], "-": []})
    for key in sorted(pair_files, key=lambda k: k.identifier):
        for strand in ("+", "-"):
            by_query[key.query][strand].extend(pair_files[key][strand])
    counts = collections.Counter({name: len(pair_lines(strands)) for name, strands in by_query.items()})
    return [
        (identifier, [line for name in names for line in pair_lines(by_query[name])])
        for identifier, names in batch_plan(query_order, counts, max_lines)
    ]


def batch_plan(
    query_order: list[str], query_lines: collections.Counter[str], max_lines: int
) -> list[tuple[str, list[str]]]:
    """Which queries share a batch, and what each batch is called. `[(identifier, [query, ...])]`.

    Pure, and separated from both writers for exactly that reason: the boundaries and the names
    are the part that has to be identical between the buffered oracle and the streaming writer,
    and the only way to be sure of that is for there to be one copy of it.

    ⛔ A BATCH IS NAMED FOR ITS FIRST QUERY, NOT ITS POSITION. `batch00000`, `batch00001`, ... look
    tidier and are a trap: add one query sequence anywhere and every boundary after it shifts, so
    the same name means a different thing between two runs. That breaks job caching, breaks resume,
    and makes a diff of two runs meaningless -- silently, because the names still line up. A first
    query name is stable under any change elsewhere in the genome, and it is unique by construction
    because each query lands in exactly one batch.

    ⚠ IT ALSO CHANGES THE DISCOVERY ORDER, which Galaxy takes as ASCII over identifiers. Batches
    then arrive in ASCII order of their first query rather than in query-file order. Nothing
    downstream depends on that -- the lastz ordering rules are WITHIN a file, and MAF/AXT blocks
    from different batches are independent -- but a consumer that assumed element order was genome
    order would be wrong, and would have been wrong about `batch00000` too, just less visibly.

    ⛔ EVERY QUERY MUST APPEAR IN `query_order`. A query the order does not name cannot be placed,
    and guessing its position is precisely the silent-wrongness batching exists to avoid.
    """
    unplaceable = sorted(set(query_lines) - set(query_order))
    if unplaceable:
        raise ValueError(
            f"{len(unplaceable)} query sequence(s) are absent from the query order and cannot be "
            f"placed: {unplaceable[:5]}"
        )
    plan: list[tuple[str, list[str]]] = []
    current: list[str] = []
    current_lines = 0
    for name in query_order:
        n = query_lines.get(name, 0)
        if n == 0:
            continue
        # ⚠ The test is on the batch SO FAR, so a single query larger than max_lines becomes an
        # oversized batch of one rather than being split. Splitting it would put half a query's
        # plus intervals after the other half's minus intervals, which lastz rejects.
        if current and current_lines + n > max_lines:
            plan.append((current[0], current))
            current, current_lines = [], 0
        current.append(name)
        current_lines += n
    if current:
        plan.append((current[0], current))
    return plan


def assign(
    commands: list[dict[str, typing.Any]], read_segments: typing.Callable[[str], list[str]]
) -> dict[PairKey, dict[str, list[str]]]:
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


def write_elements(
    elements: list[tuple[str, list[str]]], out_dir: pathlib.Path, compresslevel: int = COMPRESSLEVEL
) -> list[tuple[str, int, int]]:
    """Write `[(identifier, lines)]` as one gzipped segments file each.

    Both shapes land here -- one file per chromosome pair, and one file per query batch -- so the
    gzip level, the `.segments.gz` suffix and the manifest have a single definition.

    ⚠ IT IS THE BUFFERED REFERENCE IMPLEMENTATION, NOT THE PRODUCTION WRITER. On the per-pair
    shape `stream_pairs` is what runs; see its docstring. This one is kept because it is the
    ORACLE -- `--compare` runs both over the same input and asserts the pair files agree -- and
    because batching has no streaming form (a batch is defined by whole queries, so it cannot be
    emitted until the query is complete). Do not call it on a real bundle with `assign`: that
    needs roughly twice the bundle's segment bytes in RAM, and was OOM-killed on 15 GB.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest: list[tuple[str, int, int]] = []
    for identifier, lines in elements:
        path = out_dir / f"{identifier}.segments.gz"
        payload = "".join(line + "\n" for line in lines).encode()
        with gzip.open(path, "wb", compresslevel=compresslevel) as fh:
            fh.write(payload)
        manifest.append((identifier, len(lines), path.stat().st_size))
    return manifest


# ---------------------------------------------------------------------------------- streaming

#: Plus before minus: the pass order IS the ordering rule (see `stream_pairs`).
STRAND_ORDER = ("+", "-")

#: How many gzip writers may be open at once. Each costs ~270 KB of zlib deflate state, so the
#: default is ~70 MB at worst -- and a chromosome-level pair has only ~100 pair files, so nothing
#: is ever evicted in the case this is built for. The cap exists so a scaffold-level assembly hits
#: an LRU rather than EMFILE.
MAX_OPEN_WRITERS = 256


def _line_strand(command_args: list[str]) -> typing.Callable[[str], str]:
    """Which strand a line belongs to, resolved the way `assign` resolves it.

    ⚠ NOT ALWAYS COLUMN 7. When a command declares `--strand`, `assign` files every one of its
    lines under that strand and never looks at the column; only for `both` does the column decide.
    The directory-scan path synthesises commands with no `--strand` at all, so it is ALWAYS
    `both` there -- which is why this cannot be hoisted into a sort over commands, and why
    `stream_pairs` makes two passes instead.
    """
    strand = strand_of(command_args)
    if strand == "both":
        return lambda line: line.split("\t")[6]
    return lambda _line: strand


def stream_pairs(
    commands: list[dict[str, typing.Any]],
    read_segments: typing.Callable[[str], list[str]],
    out_dir: pathlib.Path,
    compresslevel: int = COMPRESSLEVEL,
    max_open: int = MAX_OPEN_WRITERS,
) -> list[tuple[str, int, int]]:
    """Group splits into pair files WITHOUT holding the bundle in memory.

    ⛔ WHY THIS REPLACED `assign` + `write_pairs`. Those hold every segment line of the whole
    bundle in one dict of `str`, then join each pair's lines into one more full copy before
    gzipping. Measured on real *Cannabis* splits: peak RSS = 52.4 MB + 1.96 x (input MB), a tight
    linear fit over 124-2160 MB, which extrapolates to ~28 GB for a 14.15 GB bundle. It was
    OOM-killed on a 15 GB machine. The cost is O(WHOLE BUNDLE) for an output that is per-pair --
    and the tarball path it competes with is O(1), because `package_output.py` just streams files
    into a tar.

    ▶ THE ORDERING RULE IS THE PASS ORDER. A pair file must hold its plus lines before its minus
    lines. Rather than buffer one strand to achieve that, this makes one pass per strand and
    appends straight into one gzip writer per pair, held open across both. Memory becomes
    O(number of pair files) -- tens of MB -- at the cost of reading the input twice, which is
    seconds against the hours of lastz that follow.

    ⚠ A command with a definite `--strand` is skipped entirely on the other strand's pass, so the
    command-manifest path reads each split once, not twice. Only `both` commands are read twice.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    writers: collections.OrderedDict[PairKey, gzip.GzipFile] = collections.OrderedDict()
    counts: collections.Counter[PairKey] = collections.Counter()
    started: set[PairKey] = set()

    def path_of(key: PairKey) -> pathlib.Path:
        return out_dir / f"{key.identifier}.segments.gz"

    def writer_for(key: PairKey) -> gzip.GzipFile:
        found = writers.get(key)
        if found is not None:
            writers.move_to_end(key)
            return found
        while len(writers) >= max_open:
            _, evicted = writers.popitem(last=False)
            evicted.close()
        # ▶ `mtime=0` makes the bytes a function of the content alone, so two runs of this script
        # produce identical files. The buffered implementation stamped the current time, which
        # meant its output could never be compared byte-for-byte across runs.
        # ⚠ Reopening after an eviction appends a second gzip MEMBER. That is valid gzip and every
        # reader concatenates members transparently, but it means the raw bytes of an evicted
        # pair file differ from a single-member write -- which is why `--compare` asserts on the
        # DECOMPRESSED payload rather than on the file's bytes.
        mode = "ab" if key in started else "wb"
        started.add(key)
        opened = gzip.GzipFile(filename=path_of(key), mode=mode, compresslevel=compresslevel, mtime=0)
        writers[key] = opened
        return opened

    try:
        for wanted in STRAND_ORDER:
            for command in commands:
                args = command.get("args", [])
                segments = arg_value(args, "--segments=")
                if segments is None:
                    continue
                declared = strand_of(args)
                if declared not in ("both", wanted):
                    continue  # nothing of this command belongs to this pass
                lines = read_segments(segments)
                if not lines:
                    continue
                query = query_of(lines)
                strand_of_line = _line_strand(args)
                for target, run in split_by_target(lines):
                    key = PairKey(target, query)
                    for line in run:
                        if strand_of_line(line) != wanted:
                            continue
                        handle = writer_for(key)
                        handle.write(line.encode())
                        handle.write(b"\n")
                        counts[key] += 1
    finally:
        for handle in writers.values():
            handle.close()

    return [
        (key.identifier, counts[key], path_of(key).stat().st_size)
        for key in sorted(started, key=lambda k: k.identifier)
    ]


def stream_batches(
    commands: list[dict[str, typing.Any]],
    read_segments: typing.Callable[[str], list[str]],
    out_dir: pathlib.Path,
    query_order: list[str],
    max_lines: int,
    *,
    compresslevel: int = COMPRESSLEVEL,
) -> list[tuple[str, int, int]]:
    """Write query batches WITHOUT holding the bundle in memory -- `stream_pairs` for the batched shape.

    ⛔ WHY THIS EXISTS AT ALL. `batch_pairs` needs `assign`, which holds every segment line of the
    whole bundle in RAM: measured peak RSS = 52.4 MB + 1.96 x (input MB), extrapolating to ~28 GB
    for a 14.15 GB bundle, and OOM-killed on a 15 GB machine. That is the defect `stream_pairs`
    was written to remove, and a batched path built on `assign` would reintroduce it on exactly
    the shape intended for production.

    ▶ THE PASS ORDER IS THE ORDERING RULE, ONE LEVEL DOWN FROM `stream_pairs`. There, the outer
    loop is the strand. Here it CANNOT be: a batch holds several queries, and a strand-outermost
    pass would write Q1+ Q2+ Q1- Q2-, so each query name recurs non-contiguously and rule 1 is
    broken. lastz enforces rule 1 loudly -- a query out of order exits 1 -- so the loops are
    query-outermost, strand-inner: Q1+ Q1- Q2+ Q2-.

    ⚠ IT READS THE INPUT THREE TIMES, not twice. Batch boundaries depend on each query's total
    line count, and a boundary has to be known before the first line of that batch is written, so
    an indexing pass comes first. It keeps counts, never lines, so its memory is one split. Three
    passes over local disk is minutes against the hours of lastz that follow, and the alternative
    is the 28 GB.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- pass 1: which query each split belongs to, and how many lines it contributes ---
    by_query: dict[str, list[dict[str, typing.Any]]] = collections.defaultdict(list)
    query_lines: collections.Counter[str] = collections.Counter()
    for command in commands:
        segments = arg_value(command.get("args", []), "--segments=")
        if segments is None:
            continue
        lines = read_segments(segments)
        if not lines:
            continue
        # ⚠ `query_of` raises on a split spanning two queries. Measured over all 5,177 splits of a
        # real bundle: none do. It is a bug upstream if one ever does, not a case to handle here.
        by_query[query_of(lines)].append(command)
        query_lines[query_of(lines)] += len(lines)

    manifest: list[tuple[str, int, int]] = []
    for identifier, names in batch_plan(query_order, query_lines, max_lines):
        path = out_dir / f"{identifier}.segments.gz"
        written = 0
        # `mtime=0` for the same reason as in `stream_pairs`: the bytes are then a function of the
        # content alone, so two runs produce identical files and can be compared as bytes.
        with gzip.GzipFile(filename=path, mode="wb", compresslevel=compresslevel, mtime=0) as handle:
            for name in names:
                for wanted in STRAND_ORDER:
                    for command in by_query[name]:
                        args = command["args"]
                        if strand_of(args) not in ("both", wanted):
                            continue  # nothing of this command belongs to this pass
                        strand_of_line = _line_strand(args)
                        for line in read_segments(require_arg(args, "--segments=")):
                            if strand_of_line(line) != wanted:
                                continue
                            handle.write(line.encode())
                            handle.write(b"\n")
                            written += 1
        manifest.append((identifier, written, path.stat().st_size))
    return manifest


def payload_digests(out_dir: pathlib.Path) -> dict[str, tuple[str, int]]:
    """`identifier -> (md5 of the DECOMPRESSED payload, line count)` for every pair file.

    ⚠ The digest is of the payload, not of the file. A gzip header carries an mtime and a member
    boundary is invisible after decompression, so two byte-different files can be the same pair
    file -- and byte equality would fail for reasons that have nothing to do with correctness.
    """
    digests: dict[str, tuple[str, int]] = {}
    for path in sorted(out_dir.glob("*.segments.gz")):
        digest = hashlib.md5()
        lines = 0
        with gzip.open(path, "rb") as handle:
            for line in handle:
                digest.update(line)
                lines += 1
        digests[path.name[: -len(".segments.gz")]] = (digest.hexdigest(), lines)
    return digests


def coverage_complaint(splits_seen: int, splits_present: int, bounded: bool) -> str | None:
    """Why `--compare`'s verdict does not cover the input, or `None` when it does.

    ⛔ WHY A GATE THAT PROVES AGREEMENT IS NOT ENOUGH. Both implementations read the SAME command
    list, so a truncated list makes them agree on the truncation. Measured: restricting the
    directory glob to `*.plus.*` drops 2,593 of 5,177 real splits -- half the bundle -- and the
    comparison passes, as does the self-test. Equivalence and coverage are different properties and
    only one of them was being checked.

    `bounded` is the caller saying it truncated ON PURPOSE (`--max-splits`), which is legitimate --
    the buffered reference implementation cannot survive a whole bundle, which is the entire reason
    the streaming one exists. An UNBOUNDED run that still misses splits is the failure.
    """
    if splits_seen == splits_present:
        return None
    missing = splits_present - splits_seen
    if bounded:
        return (
            f"bounded by --max-splits: {splits_seen} of {splits_present} split(s) compared, "
            f"{missing} not covered by this verdict"
        )
    return (
        f"COVERAGE: {splits_seen} of {splits_present} split(s) reached the comparison and "
        f"{missing} did not, with no --max-splits to explain it. Both implementations read the "
        f"same list, so they agree on whatever it omits -- this verdict does not cover the "
        f"input it was pointed at."
    )


def compare_implementations(
    commands: list[dict[str, typing.Any]],
    read_segments: typing.Callable[[str], list[str]],
    tmp: pathlib.Path,
    compresslevel: int = COMPRESSLEVEL,
    max_open: int = MAX_OPEN_WRITERS,
) -> tuple[bool, list[str]]:
    """THE GATE. Run both implementations over one input and report every disagreement.

    Returns `(ok, complaints)`. This is what licenses replacing the buffered implementation: not
    an argument that the streaming one is equivalent, but a run of both on the same bytes.
    """
    buffered_dir, streamed_dir = tmp / "buffered", tmp / "streamed"
    write_pairs(assign(commands, read_segments), buffered_dir, compresslevel)
    stream_pairs(commands, read_segments, streamed_dir, compresslevel, max_open)

    buffered, streamed = payload_digests(buffered_dir), payload_digests(streamed_dir)
    complaints: list[str] = []
    for missing in sorted(set(buffered) - set(streamed)):
        complaints.append(f"{missing}: written by the buffered path, absent from the streamed one")
    for extra in sorted(set(streamed) - set(buffered)):
        complaints.append(f"{extra}: written by the streamed path, absent from the buffered one")
    for name in sorted(set(buffered) & set(streamed)):
        if buffered[name] != streamed[name]:
            complaints.append(
                f"{name}: buffered md5={buffered[name][0]} lines={buffered[name][1]}, "
                f"streamed md5={streamed[name][0]} lines={streamed[name][1]}"
            )
    return not complaints, complaints


def write_pairs(
    pair_files: dict[PairKey, dict[str, list[str]]], out_dir: pathlib.Path, compresslevel: int = COMPRESSLEVEL
) -> list[tuple[str, int, int]]:
    return write_elements(
        [(key.identifier, pair_lines(pair_files[key])) for key in sorted(pair_files, key=lambda k: k.identifier)],
        out_dir,
        compresslevel,
    )


# --------------------------------------------------------------------------------------- verify

MAF_BLOCK_START = re.compile(r"^a score=")

#: A lastz sequence-file ACTION restricting a 2bit to the names in a file, e.g.
#: `query.2bit[subset=query_block0.name]`. `--verify-batches` strips it; see there for why.
SUBSET_ACTION = re.compile(r"\[subset=[^\]]*\]")


def maf_blocks(path: pathlib.Path) -> collections.Counter[tuple[str, ...]]:
    """MAF blocks as a multiset, so two runs can be compared regardless of block order."""
    blocks: collections.Counter[tuple[str, ...]] = collections.Counter()
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


def run_lastz(
    lastz: str,
    target: str,
    query: str,
    segments: pathlib.Path,
    *,
    out: pathlib.Path,
    extra: list[str],
    workdir: pathlib.Path,
) -> float:
    argv = [
        lastz,
        target,
        query,
        "--allocate:traceback=1.99G",
        "--format=maf-",
        f"--segments={segments}",
        f"--output={out}",
        *extra,
    ]
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
    bundle = pathlib.Path(args.bundle).resolve()
    workdir = bundle / "galaxy" / "files"
    commands = [
        json.loads(line) for line in (bundle / "galaxy" / "commands.json").read_text().splitlines() if line.strip()
    ]

    def read_segments(name: str) -> list[str]:
        return (workdir / name).read_text().splitlines()

    want_target, want_query = [*args.pair.split(","), ""][:2]
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
        by_size = sorted(selected, key=lambda c: (workdir / require_arg(c["args"], "--segments=")).stat().st_size)
        keep = {id(c) for c in by_size[: args.smallest]}
        selected = [c for c in selected if id(c) in keep]
    elif args.max_splits:
        selected = selected[: args.max_splits]
    print(f"pair {want_target} x {want_query}: {len(selected)} split(s)", file=sys.stderr)

    # ⛔ THROUGH `stream_pairs`, NOT the buffered pair. This arm is what compares a pair file
    # against its splits with real lastz, so it has to build that pair file the way production
    # does -- verifying an implementation that no longer runs would be worse than not verifying.
    tmp = pathlib.Path(tempfile.mkdtemp(prefix="pair_verify_", dir=args.tmpdir))
    manifest = stream_pairs(selected, read_segments, tmp / "pairs")
    if len(manifest) != 1:
        print(
            f"  note: {len(manifest)} pair files in this selection: {[i for i, _, _ in manifest]}",
            file=sys.stderr,
        )
    for identifier, n_lines, n_bytes in manifest:
        print(f"  pair file {identifier}: {n_lines:,} segments, {n_bytes:,} bytes gzipped", file=sys.stderr)

    target_spec = require_arg(selected[0]["args"], "--target=")
    query_spec = require_arg(selected[0]["args"], "--query=")
    shared = [
        a
        for a in selected[0]["args"]
        if not a.startswith(("--target=", "--query=", "--segments=", "--output=", "--strand=", "--format="))
    ]

    # --- arm A: one lastz per pair file, no --strand (default is both; measured identical) ---
    per_pair: collections.Counter[tuple[str, ...]] = collections.Counter()
    pair_seconds = 0.0
    for identifier, _, _ in manifest:
        plain = tmp / f"{identifier}.segments"
        with gzip.open(tmp / "pairs" / f"{identifier}.segments.gz", "rb") as src, plain.open("wb") as dst:
            shutil.copyfileobj(src, dst)  # the consumer's gunzip step, in miniature
        out = tmp / f"{identifier}.pairfile.maf-"
        pair_seconds += run_lastz(args.lastz, target_spec, query_spec, plain, out=out, extra=shared, workdir=workdir)
        per_pair += maf_blocks(out)

    # --- arm B: the original splits, each with its own --strand ---
    per_split: collections.Counter[tuple[str, ...]] = collections.Counter()
    split_seconds = 0.0
    for command in selected:
        segments_path = workdir / require_arg(command["args"], "--segments=")
        strand = [a for a in command["args"] if a.startswith("--strand=")]
        out = tmp / (pathlib.Path(require_arg(command["args"], "--output=")).name + ".split")
        split_seconds += run_lastz(
            args.lastz, target_spec, query_spec, segments_path, out=out, extra=[*shared, *strand], workdir=workdir
        )
        per_split += maf_blocks(out)

    print(f"\n{'':<22}{'blocks':>10}{'seconds':>10}", file=sys.stderr)
    print(f"{'pair file (1 lastz)':<22}{sum(per_pair.values()):>10,}{pair_seconds:>10.1f}", file=sys.stderr)
    print(
        f"{f'splits ({len(selected)} lastz)':<22}{sum(per_split.values()):>10,}{split_seconds:>10.1f}", file=sys.stderr
    )

    return compare_loci("pair file", per_pair, f"splits ({len(selected)})", per_split)


def verify_batches(args: argparse.Namespace) -> int:
    """Run lastz once per pair file and once per BATCH, and compare the alignments.

    ⛔ THE CLAIM `--batch-lines` RESTS ON, tested the way `--verify` tests the pair file: against
    real data and real lastz rather than argued from the manual. Packing several queries into one
    segment file is only legitimate if it changes nothing about what comes out.

    ▶ WHY IT LIVES HERE RATHER THAN IN A SCRATCH DIRECTORY. The evidence for batching was a
    harness nobody else had, in a directory that does not survive the sandbox, quoted in a commit
    message on an unmerged branch. A claim that cannot be re-run by the next person is a claim
    with no evidence a month later.

    ⚠ It needs a real lastz and a real bundle, so it SKIPS when the binary is absent rather than
    failing -- a CI run without lastz should not go red over a gate it cannot execute. The
    self-test's own batched gate (`stream_batches` against `batch_pairs`) is what runs everywhere.
    """
    lastz_path = shutil.which(args.lastz)
    if lastz_path is None:
        print(f"skip - no lastz on PATH as {args.lastz!r}; nothing to verify", file=sys.stderr)
        return 0

    bundle = pathlib.Path(args.bundle).resolve()
    workdir = bundle / "galaxy" / "files"
    commands = [
        json.loads(line) for line in (bundle / "galaxy" / "commands.json").read_text().splitlines() if line.strip()
    ]

    def read_segments(name: str) -> list[str]:
        return (workdir / name).read_text().splitlines()

    selected = [c for c in commands if arg_value(c["args"], "--segments=") is not None]
    selected.sort(key=lambda c: arg_value(c["args"], "--output=") or "")
    if args.smallest:
        # ⚠ The SMALLEST splits, not the first N -- and across every query, not one pair, because
        # a selection confined to one query cannot produce a batch that holds two. Whether merging
        # changes the answer is the same code path at any size.
        by_size = sorted(selected, key=lambda c: (workdir / require_arg(c["args"], "--segments=")).stat().st_size)
        keep = {id(c) for c in by_size[: args.smallest]}
        selected = [c for c in selected if id(c) in keep]
    if not selected:
        sys.exit("no commands with segments in this bundle")

    # ⛔ THE `subset=` ACTIONS MUST GO, ON BOTH ARMS. A production split names a block file --
    # `query.2bit[subset=query_block0.name]` -- and a batch deliberately spans blocks, so the
    # subset of whichever split happened to be first would hide part of every batch. Dropping it
    # is the measured-safe move: with the whole 2bit the output is byte-identical (md5
    # 8e060801bfbb, 1,308,851 bytes) at +888 MB peak RSS and +25% time, and it is what
    # `growler_lastz` does in production. Both arms get the same specs, or the comparison is
    # between two different experiments.
    def whole_2bit(spec: str) -> str:
        return SUBSET_ACTION.sub("", spec)

    query_spec = whole_2bit(require_arg(selected[0]["args"], "--query="))
    query_2bit = pathlib.Path(args.query_2bit) if args.query_2bit else workdir / query_spec.split("[")[0]
    order = read_2bit_order(query_2bit)

    tmp = pathlib.Path(tempfile.mkdtemp(prefix="batch_verify_", dir=args.tmpdir))
    pairs = stream_pairs(selected, read_segments, tmp / "pairs")
    batch_lines = args.batch_lines or sum(n for _, n, _ in pairs)
    batches = stream_batches(selected, read_segments, tmp / "batches", order, batch_lines)

    # ⛔ THE TEST MUST BE ABLE TO SAY SOMETHING. If every batch holds one query, the two arms are
    # the same files under different names and the comparison passes without testing anything.
    # This is the failure mode a green run would otherwise hide, so it is checked before lastz is
    # paid for, not after.
    def queries_in(identifier: str) -> set[str]:
        with gzip.open(tmp / "batches" / f"{identifier}.segments.gz", "rt") as handle:
            return {line.split("\t")[3] for line in handle if line.strip()}

    packed = max((len(queries_in(i)) for i, _, _ in batches), default=0)
    print(
        f"{len(selected)} split(s) -> {len(pairs)} pair file(s) -> {len(batches)} batch(es); "
        f"largest batch holds {packed} query sequence(s)",
        file=sys.stderr,
    )
    if packed < 2:
        print(
            "not ok - no batch holds more than one query, so this run does not test batching. "
            "Raise --smallest, or lower --batch-lines.",
            file=sys.stderr,
        )
        return 1

    target_spec = whole_2bit(require_arg(selected[0]["args"], "--target="))
    shared = [
        a
        for a in selected[0]["args"]
        if not a.startswith(("--target=", "--query=", "--segments=", "--output=", "--strand=", "--format="))
    ]

    def arm(manifest: list[tuple[str, int, int]], subdir: str) -> collections.Counter[tuple[str, ...]]:
        blocks: collections.Counter[tuple[str, ...]] = collections.Counter()
        for identifier, _, _ in manifest:
            plain = tmp / f"{subdir}.{identifier}.segments"
            with gzip.open(tmp / subdir / f"{identifier}.segments.gz", "rb") as src, plain.open("wb") as dst:
                shutil.copyfileobj(src, dst)  # the consumer's gunzip step, in miniature
            out = tmp / f"{subdir}.{identifier}.maf-"
            run_lastz(args.lastz, target_spec, query_spec, plain, out=out, extra=shared, workdir=workdir)
            blocks += maf_blocks(out)
        return blocks

    return compare_loci(
        f"batches ({len(batches)})", arm(batches, "batches"), f"pair files ({len(pairs)})", arm(pairs, "pairs")
    )


def locus_of(block: tuple[str, ...]) -> tuple[str, str, str]:
    """A MAF block reduced to what must match: (score line, target start, query start).

    ⛔ THIS IS WHY THE COMPARISON IS NOT ON BLOCK TEXT. lastz's own manual warns that the same
    alignment may come back with "minor variations such as shifting of equally-scoring gap
    placements", so two runs can emit one alignment -- same score, same coordinates, same length --
    as two different STRINGS. Comparing full block text calls that a difference, which on a
    238,076-segment pair file reported three alignments as "missing from the pair file" when every
    one of them was present.
    """
    return (block[0], block[1].split()[2], block[2].split()[2])


def compare_loci(
    label_a: str,
    blocks_a: collections.Counter[tuple[str, ...]],
    label_b: str,
    blocks_b: collections.Counter[tuple[str, ...]],
) -> int:
    """Report whether two lastz arms found the same alignments. 0 if they did.

    Shared by every `--verify` mode, because the thing each of them is really asserting is the
    same one: regrouping the SAME anchor segments into different files changes which process
    finds an alignment, and must not change whether one is found.
    """
    loci_a = collections.Counter(locus_of(b) for b in blocks_a.elements())
    loci_b = collections.Counter(locus_of(b) for b in blocks_b.elements())

    lost = set(loci_b) - set(loci_a)
    gained = set(loci_a) - set(loci_b)
    # ▶ DUPLICATES ARE A FINDING, NOT AN ERROR. When the anchors for one alignment fall in two
    # different files, each lastz finds it independently and both report it; one lastz that sees
    # them together reports it once. So the side with fewer files DE-DUPLICATES what the other
    # emits twice, and the locus sets still agree.
    dup_a = {k: v for k, v in loci_a.items() if v > 1}
    dup_b = {k: v for k, v in loci_b.items() if v > 1}

    print(f"\n{'':<22}{'blocks':>10}{'loci':>10}", file=sys.stderr)
    print(f"{label_a:<22}{sum(loci_a.values()):>10,}{len(loci_a):>10,}", file=sys.stderr)
    print(f"{label_b:<22}{sum(loci_b.values()):>10,}{len(loci_b):>10,}", file=sys.stderr)

    if not lost and not gained:
        extra_a = sum(v - 1 for v in dup_a.values())
        extra_b = sum(v - 1 for v in dup_b.values())
        note = f"; duplicate blocks: {label_a} {extra_a}, {label_b} {extra_b}" if (extra_a or extra_b) else ""
        print(f"\nok - identical alignment loci{note}", file=sys.stderr)
        return 0
    print(f"\nnot ok - {len(lost)} locus/loci only in {label_b}, {len(gained)} only in {label_a}", file=sys.stderr)
    for loc in list(lost)[:3]:
        print(f"  {label_b}-only locus:", loc, file=sys.stderr)
    for loc in list(gained)[:3]:
        print(f"  {label_a}-only locus:", loc, file=sys.stderr)
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
        return f"{t}\t{start}\t{start + 9}\t{q}\t{start}\t{start + 9}\t{s}\t3000"

    # a split holding one target, and one holding three contiguous runs
    single = [seg("T1", "Q1", "+", i) for i in range(3)]
    multi = [seg("T1", "Q1", "+", 0)] + [seg("T2", "Q1", "+", 1)] * 2 + [seg("T3", "Q1", "+", 2)]
    check("single-target split is one run", [t for t, _ in split_by_target(single)], ["T1"])
    check("multi-target split cuts into contiguous runs", [t for t, _ in split_by_target(multi)], ["T1", "T2", "T3"])
    check("cutting loses no lines", sum(len(r) for _, r in split_by_target(multi)), len(multi))

    check("query_of finds the single query", query_of(single), "Q1")
    try:
        query_of([*single, seg("T1", "Q2", "+", 9)])
        check("two queries must raise", "no raise", "ValueError")
    except ValueError:
        check("two queries must raise", "ValueError", "ValueError")

    check("--strand=minus reads as '-'", strand_of(["--strand=minus"]), "-")
    check("absent --strand means both", strand_of(["--ydrop=1"]), "both")

    # ▶ THE ORDERING RULE. Build a pair file from a minus command given BEFORE a plus command and
    # assert the written order is still plus-then-minus -- this is the invariant the lastz
    # manual imposes, and the one a naive concatenation in command order would break.
    store = {"m.segments": [seg("T1", "Q1", "-", 5)], "p.segments": [seg("T1", "Q1", "+", 1)]}
    cmds = [{"args": ["--segments=m.segments", "--strand=minus"]}, {"args": ["--segments=p.segments", "--strand=plus"]}]
    pair_files = assign(cmds, lambda n: store[n])
    check("one pair file from one chromosome pair", [k.identifier for k in pair_files], ["T1__Q1"])
    written = pair_lines(pair_files[PairKey("T1", "Q1")])
    check("plus is written before minus, whatever the command order", [ln.split("\t")[6] for ln in written], ["+", "-"])

    # a multi-target split lands in as many pair files as it has targets, losing nothing
    pair_files2 = assign([{"args": ["--segments=x", "--strand=plus"]}], lambda n: multi)
    check(
        "multi-target split fans into 3 pair files",
        sorted(k.identifier for k in pair_files2),
        ["T1__Q1", "T2__Q1", "T3__Q1"],
    )
    check("  and every line survives the fan-out", sum(len(pair_lines(v)) for v in pair_files2.values()), len(multi))

    # ▶ THE DIRECTORY-SCAN PATH: no --strand in the synthesised command, so the strand has to
    # come from column 7 of each line. A pair file built this way must match one built from commands.
    mixed = [seg("T1", "Q1", "-", 5), seg("T1", "Q1", "+", 1)]
    from_file = assign([{"args": ["--segments=x"]}], lambda n: mixed)
    check(
        "strand read from column 7 when --strand is absent",
        [ln.split("\t")[6] for ln in pair_lines(from_file[PairKey("T1", "Q1")])],
        ["+", "-"],
    )
    from_cmds = assign(
        [{"args": ["--segments=m", "--strand=minus"]}, {"args": ["--segments=p", "--strand=plus"]}],
        lambda n: [mixed[0]] if n == "m" else [mixed[1]],
    )
    check(
        "  and agrees with the command-driven path",
        pair_lines(from_file[PairKey("T1", "Q1")]),
        pair_lines(from_cmds[PairKey("T1", "Q1")]),
    )

    # ▶ QUERY BATCHING. The invariant under test is the manual's first rule: queries appear in
    # QUERY-FILE order, whatever order the pair files were keyed in. `assign` keys them
    # alphabetically, so an order that is deliberately NOT alphabetical is the only honest fixture.
    #
    # ⚠ ONE SPLIT STILL HOLDS ONE QUERY, and `query_of` still refuses otherwise. That invariant is
    # about KegAlign's INPUT splits (measured: 0 of 5,177 span two queries) and is untouched. What
    # batching relaxes is the OUTPUT file, which the manual never restricted to one query.
    order = ["Qc", "Qa", "Qb"]
    splits = {
        "a": [seg("T1", "Qa", "+", 1), seg("T1", "Qa", "-", 2)],
        "b": [seg("T2", "Qb", "+", 3)],
        "c": [seg("T1", "Qc", "+", 4)],
    }
    keyed = assign([{"args": [f"--segments={n}"]} for n in ("a", "b", "c")], lambda n: splits[n])
    one = batch_pairs(keyed, order, max_lines=100)
    check("a generous batch size yields one element", [i for i, _ in one], ["Qc"])
    check(
        "queries are emitted in QUERY-FILE order, not key order",
        [ln.split("\t")[3] for ln in one[0][1]],
        ["Qc", "Qa", "Qa", "Qb"],
    )
    check(
        "  and within a query, plus still precedes minus",
        [ln.split("\t")[6] for ln in one[0][1] if ln.split("\t")[3] == "Qa"],
        ["+", "-"],
    )
    check("  and a batch may span several targets", sorted({ln.split("\t")[0] for ln in one[0][1]}), ["T1", "T2"])

    # a batch boundary never falls inside a query: Qa owns two lines and they stay together
    split = batch_pairs(keyed, order, max_lines=2)
    check("a tight batch size splits into several elements", len(split) > 1, True)
    check("no line is lost across batches", sum(len(v) for _, v in split), 4)
    for identifier, lines in split:
        names = [ln.split("\t")[3] for ln in lines]
        check(f"  {identifier} holds whole queries only", names == sorted(names, key=order.index), True)
    owners: dict[str, set[int]] = collections.defaultdict(set)
    for i, (_, lines) in enumerate(split):
        for ln in lines:
            owners[ln.split("\t")[3]].add(i)
    check(
        "every query lands in exactly ONE batch",
        {q: len(v) for q, v in sorted(owners.items())},
        {"Qa": 1, "Qb": 1, "Qc": 1},
    )
    check(
        "  concatenating the batches reproduces the single-batch order",
        [ln for _, lines in split for ln in lines],
        one[0][1],
    )

    # ⛔ a query the order cannot place must raise rather than be guessed at
    try:
        batch_pairs(keyed, ["Qa", "Qb"], max_lines=100)
        check("an unplaceable query must raise", "no raise", "ValueError")
    except ValueError:
        check("an unplaceable query must raise", "ValueError", "ValueError")

    # ⛔ THE NAMES MUST BE STABLE UNDER AN INSERTION ELSEWHERE. This is the whole reason a batch is
    # named for its first query rather than for its position: add a query at the FRONT and every
    # positional name would shift by one, so `batch00001` would mean a different set of queries in
    # two runs of the same pipeline -- with nothing to notice it, because the names still line up.
    counts = collections.Counter({"Qc": 1, "Qa": 2, "Qb": 1})
    before = batch_plan(order, counts, max_lines=2)
    after = batch_plan(["Qd", *order], collections.Counter({"Qd": 1, **counts}), max_lines=2)
    check("batches are named for their first query", [i for i, _ in before], ["Qc", "Qa", "Qb"])
    # Positionally these would be batch00000/1/2 before and 00000/1/2 after, with DIFFERENT
    # contents: the insertion pushes Qc into the first batch. By first query, the two batches that
    # did not change still answer to the same name and still hold the same queries.
    survivors = sorted(set(dict(before)) & set(dict(after)))
    check("  an insertion at the front leaves the later batches named as they were", survivors, ["Qa", "Qb"])
    check(
        "  and each surviving name still means the same queries",
        [dict(after)[n] for n in survivors],
        [dict(before)[n] for n in survivors],
    )
    check(
        "  while the batch the insertion landed in is a new name", sorted(set(dict(after)) - set(dict(before))), ["Qd"]
    )

    # ------------------------------------------------------------- streaming == buffered, batched
    # ▶ THE SECOND GATE. `stream_batches` is what production runs; `batch_pairs` is the oracle it
    # is checked against. Without this, the streaming batched writer has no test at all -- and it
    # is the one that must not fall back to `assign`, whose memory was the reason for #65.
    batched_commands = [{"args": [f"--segments={n}"]} for n in ("a", "b", "c")]
    for label, limit in (("one batch", 100), ("several batches", 2)):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            streamed = stream_batches(batched_commands, lambda n: splits[n], root / "s", order, limit)
            write_elements(batch_pairs(keyed, order, limit), root / "b")
            check(
                f"{label}: streamed batches match the buffered oracle",
                payload_digests(root / "s"),
                payload_digests(root / "b"),
            )
            check(f"  {label}: and the manifest counts the lines it wrote", sum(n for _, n, _ in streamed), 4)

    # The batch writer honours its compresslevel, pinned on the BYTES like the pair-file checks
    # above: nothing else asserted what level a batch is written at.
    #
    # ⚠ THIS DOES NOT COVER main()'s CALL SITE, and saying so is the point. `stream_batches` made
    # `compresslevel` keyword-only to satisfy a lint, main() still passed it positionally, and CI
    # failed on `mypy --strict` while this suite stayed green. Adding this check does not change
    # that: it calls `stream_batches` directly, so it would pass with main() still broken --
    # mutation-tested, and it does. Nothing here exercises main(), so mypy is the only thing
    # standing between a signature change and its callers. Do not read a green self-test as
    # covering the CLI.
    with tempfile.TemporaryDirectory() as tmp:
        root = pathlib.Path(tmp)
        stream_batches(batched_commands, lambda n: splits[n], root, order, 100, compresslevel=COMPRESSLEVEL)
        xfl = sorted({path.read_bytes()[8] for path in root.glob("*.segments.gz")})
        check("batch files carry the fastest-deflate XFL byte too", xfl, [4])

    # ⛔ AND THAT GATE MUST BE ABLE TO FAIL. A strand-outermost pass is the plausible wrong way to
    # write a multi-query batch: it is what `stream_pairs` does one level up, and here it yields
    # Qc+ Qa+ Qb+ Qa-, in which Qa recurs after Qb. That breaks rule 1, which lastz enforces with
    # rc=1. Plant exactly that file and assert the digests diverge.
    with tempfile.TemporaryDirectory() as tmp:
        root = pathlib.Path(tmp)
        stream_batches(batched_commands, lambda n: splits[n], root / "good", order, 100)
        strand_outermost = [splits["c"][0], splits["a"][0], splits["b"][0], splits["a"][1]]
        write_elements([("Qc", strand_outermost)], root / "bad")
        check(
            "a strand-outermost batch IS caught",
            payload_digests(root / "good") != payload_digests(root / "bad"),
            True,
        )

    # ▶ THE 2BIT ORDER READER, against bytes laid out by hand from the format description.
    def fake_2bit(names: list[str], endian: str = "<", version: int = 0) -> bytes:
        blob = struct.pack(endian + "4I", 0x1A412743, version, len(names), 0)
        for n in names:
            blob += bytes([len(n)]) + n.encode() + struct.pack(endian + "I", 0)
        return blob

    with tempfile.TemporaryDirectory() as tmp:
        p = pathlib.Path(tmp) / "q.2bit"
        p.write_bytes(fake_2bit(["chr3", "chr1", "chr10"]))
        check("2bit order is file order, not sorted", read_2bit_order(p), ["chr3", "chr1", "chr10"])
        p.write_bytes(fake_2bit(["chr3", "chr1"], endian=">"))
        check("  and a big-endian 2bit reads the same", read_2bit_order(p), ["chr3", "chr1"])
        p.write_bytes(fake_2bit(["chr1"], version=1))
        try:
            read_2bit_order(p)
            check("  version 1 must refuse", "no raise", "ValueError")
        except ValueError:
            check("  version 1 must refuse", "ValueError", "ValueError")
        p.write_bytes(b"not a 2bit file at all")
        try:
            read_2bit_order(p)
            check("  a non-2bit must refuse", "no raise", "ValueError")
        except ValueError:
            check("  a non-2bit must refuse", "ValueError", "ValueError")

    # ⚠ the separator must not occur in the payload it separates
    check("separator is two underscores", PAIR_SEP, "__")
    check(
        "a single-underscore name still round-trips",
        PairKey("chr_1", "chr_2").identifier.split(PAIR_SEP),
        ["chr_1", "chr_2"],
    )

    # ------------------------------------------------------------------ streaming == buffered
    # ▶ THE GATE, run on synthetic input so CI exercises it without a bundle. Each case is one
    # the streaming rewrite could plausibly break: the ordering rule, a multi-target fan-out, a
    # command whose declared strand disagrees with column 7, and the LRU eviction path.
    splits = {
        "p1.segments": [seg("T1", "Q1", "+", 1), seg("T2", "Q1", "+", 2)],
        "m1.segments": [seg("T1", "Q1", "-", 3)],
        "p2.segments": [seg("T1", "Q1", "+", 4)],
        "b1.segments": [seg("T3", "Q2", "-", 5), seg("T3", "Q2", "+", 6)],
    }
    mixed_commands = [
        {"args": ["--segments=m1.segments", "--strand=minus"]},
        {"args": ["--segments=p1.segments", "--strand=plus"]},
        {"args": ["--segments=b1.segments"]},
        {"args": ["--segments=p2.segments", "--strand=plus"]},
    ]
    with tempfile.TemporaryDirectory() as tmp:
        ok, complaints = compare_implementations(mixed_commands, lambda n: splits[n], pathlib.Path(tmp))
        check("streaming matches buffered on mixed strands and fan-out", (ok, complaints), (True, []))

        # ⛔ EVICTION. With one writer allowed, every pair file is reopened and appended to, so
        # each becomes multi-member gzip. The PAYLOAD must be unchanged -- this is the assertion
        # that makes `--max-open` safe rather than merely present.
        ok_lru, complaints_lru = compare_implementations(
            mixed_commands, lambda n: splits[n], pathlib.Path(tmp), max_open=1
        )
        check("eviction to multi-member gzip preserves the payload", (ok_lru, complaints_lru), (True, []))

    # ⛔ COVERAGE, THE HOLE THE FIRST VERSION OF THIS GATE HAD. Restricting the input glob to
    # `*.plus.*` drops 2,593 of 5,177 real splits and BOTH implementations agree on the remainder,
    # so `--compare` passed and so did this self-test. Agreement and coverage are separate
    # properties; these assert the second one is now checked and can fail.
    check("full coverage is silent", coverage_complaint(5177, 5177, False), None)
    check(
        "a bounded run says what it did not cover",
        "bounded by --max-splits" in (coverage_complaint(400, 5177, True) or ""),
        True,
    )
    check(
        "an UNBOUNDED shortfall is a COVERAGE failure",
        (coverage_complaint(2584, 5177, False) or "").startswith("COVERAGE:"),
        True,
    )
    # ▶ and the two must not be confused: the same shortfall is a warning when asked for and a
    # failure when not, which is the whole distinction `--max-splits` encodes.
    check(
        "the same shortfall reads differently bounded vs not",
        coverage_complaint(2584, 5177, True) != coverage_complaint(2584, 5177, False),
        True,
    )

    # ⛔ AND THE GATE MUST BE ABLE TO FAIL. A comparison that cannot detect a planted difference
    # would pass forever and license anything. Break the ordering rule in a copy of the streaming
    # output and assert the digests diverge.
    with tempfile.TemporaryDirectory() as tmp:
        root = pathlib.Path(tmp)
        good, bad = root / "good", root / "bad"
        stream_pairs(mixed_commands, lambda n: splits[n], good)
        reversed_pairs = {k: {"+": v["-"], "-": v["+"]} for k, v in assign(mixed_commands, lambda n: splits[n]).items()}
        write_pairs(reversed_pairs, bad)
        check(
            "a planted strand-order swap IS caught",
            payload_digests(good) != payload_digests(bad),
            True,
        )

    # ⛔ THE GZIP LEVEL, ASSERTED ON THE BYTES. Checking `COMPRESSLEVEL == 1` would only restate
    # the constant; what matters is that the level reaches the deflate call in BOTH writers. The
    # gzip header's XFL byte (offset 8) records it: CPython writes 4 for level 1, 2 for level 9
    # and 0 for everything between, so this distinguishes 1 from the default 6 that it keeps
    # getting "corrected" back to. See COMPRESSLEVEL for why the level is 1.
    with tempfile.TemporaryDirectory() as tmp:
        root = pathlib.Path(tmp)
        stream_pairs(mixed_commands, lambda n: splits[n], root / "streamed")
        write_pairs(assign(mixed_commands, lambda n: splits[n]), root / "buffered")
        for impl in ("streamed", "buffered"):
            xfl = sorted({path.read_bytes()[8] for path in (root / impl).glob("*.segments.gz")})
            check(f"{impl} pair files carry the fastest-deflate XFL byte", xfl, [4])

    if failures:
        print(f"\n{len(failures)} test(s) failed")
        return 1
    print("\nall tests passed")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--commands", help="commands.json to group (optional; segments carry target, query and strand themselves)"
    )
    parser.add_argument("--segments-dir", default=".", help="directory holding the .segments files")
    parser.add_argument("--out", default="pairs", help="directory to write pair files into")
    parser.add_argument("--compresslevel", type=int, default=COMPRESSLEVEL)
    parser.add_argument(
        "--batch-lines",
        type=int,
        default=0,
        help="pack whole query sequences into batches of about this many segment lines "
        "(0 = off, one file per chromosome pair). Needs --query-2bit. A single query is never "
        "split, so one larger than this becomes a batch of its own.",
    )
    parser.add_argument(
        "--query-2bit",
        help="the query .2bit, read for its sequence ORDER only. Required by --batch-lines, "
        "because a batch holding several queries must list them in query-file order.",
    )
    parser.add_argument("--verify", action="store_true", help="compare a pair file against its splits with real lastz")
    parser.add_argument(
        "--verify-batches",
        action="store_true",
        help="compare query BATCHES against the per-pair files with real lastz. Needs --bundle; "
        "skips when no lastz is on PATH.",
    )
    parser.add_argument("--bundle", help="extracted bundle, for --verify")
    parser.add_argument("--pair", help="TARGET,QUERY to verify, e.g. EH23a.chr9,EH23b.chrX")
    parser.add_argument(
        "--compare",
        action="store_true",
        help="THE GATE: run the buffered and streaming implementations over the same splits and "
        "assert the pair files agree. Bound it with --max-splits; the buffered one needs ~2x the "
        "input in RAM.",
    )
    parser.add_argument(
        "--max-open",
        type=int,
        default=MAX_OPEN_WRITERS,
        help=f"concurrent gzip writers before an LRU eviction (default {MAX_OPEN_WRITERS})",
    )
    parser.add_argument("--max-splits", type=int, default=0, help="first N splits in --verify/--compare (0 = all)")
    parser.add_argument("--smallest", type=int, default=0, help="the N SMALLEST splits, for a quick --verify")
    parser.add_argument("--lastz", default="lastz")
    parser.add_argument("--tmpdir", default=None)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    if args.self_test:
        return self_test()
    if args.verify_batches:
        if not args.bundle:
            parser.error("--verify-batches needs --bundle")
        return verify_batches(args)
    if args.verify:
        if not (args.bundle and args.pair):
            parser.error("--verify needs --bundle and --pair")
        return verify(args)

    seg_dir = pathlib.Path(args.segments_dir)
    if args.commands:
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

    def read_segments(name: str) -> list[str]:
        return (seg_dir / name).read_text().splitlines()

    if args.compare:
        # ⚠ Deliberately bounded. The gate runs the BUFFERED implementation too, so it can only be
        # run over an input small enough for that one to survive -- which is the whole reason the
        # streaming one exists. Equivalence is established where both fit and then relied on.
        chosen = commands[: args.max_splits] if args.max_splits else commands
        # ⛔ COVERAGE IS CHECKED AGAINST THE DIRECTORY, NOT AGAINST `commands`. Deriving both from
        # the same glob would make the check agree with whatever the glob omitted -- which is the
        # exact defect it exists to catch.
        present = len(list(seg_dir.glob("*.segments")))
        shortfall = coverage_complaint(len(chosen), present, bool(args.max_splits))
        print(f"comparing both implementations over {len(chosen)} of {present} split(s)")
        with tempfile.TemporaryDirectory(dir=args.tmpdir) as tmp:
            ok, complaints = compare_implementations(
                chosen, read_segments, pathlib.Path(tmp), args.compresslevel, args.max_open
            )
        for complaint in complaints:
            print(f"  ⛔ {complaint}")
        print("pair files agree" if ok else f"{len(complaints)} disagreement(s)")
        if shortfall:
            print(f"  {'⚠' if args.max_splits else '⛔'} {shortfall}")
        # An unbounded run that did not reach every split fails even when the two agree: agreement
        # over part of the input is not the claim `--compare` is asked for.
        return 0 if ok and not (shortfall and not args.max_splits) else 1

    if args.batch_lines:
        if not args.query_2bit:
            parser.error("--batch-lines needs --query-2bit: the batch order is the query file's order")
        order = read_2bit_order(pathlib.Path(args.query_2bit))
        manifest = stream_batches(
            commands, read_segments, pathlib.Path(args.out), order, args.batch_lines, compresslevel=args.compresslevel
        )
        shape = "batches"
    else:
        manifest = stream_pairs(commands, read_segments, pathlib.Path(args.out), args.compresslevel, args.max_open)
        shape = "pair_files"
    total_lines = sum(n for _, n, _ in manifest)
    total_bytes = sum(b for _, _, b in manifest)
    print(f"{len(manifest)} {shape}, {total_lines:,} segments, {total_bytes / 1e9:.2f} GB gzipped")
    return 0


if __name__ == "__main__":
    sys.exit(main())
