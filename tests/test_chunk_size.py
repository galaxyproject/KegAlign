#!/usr/bin/env python3
"""estimate_chunk_size() must survive a sparse alignment.

Run with:  python3 tests/test_chunk_size.py

statistics.quantiles() needs at least two data points on Python 3.10-3.12 and at
least one on 3.13+. runner.py guarded only the "<7 data points" case, so a run
whose alignment produced a single .segments file reached it with 1 and died with
an unhandled StatisticsError -- after the GPU work was already finished:

    File "runner.py", line 378, in estimate_chunk_size
      chunk_size = int(statistics.quantiles(fdict.values())[1] // line_size)
    statistics.StatisticsError: must have at least two data points

Observed on a real run: two sequences with a literal X every 500 bases, sparse
enough to yield one segments file. The conda-forge kegalign recipe pins
python 3.12.*, so that is the interpreter the shipped package uses.

This imports the real estimate_chunk_size() and drives it against temporary
directories of .segments files. It deliberately does NOT re-implement the
arithmetic: a copy would keep passing while runner.py drifted away from it.
"""

import argparse
import importlib.util
import os
import pathlib
import sys
import tempfile

HERE = pathlib.Path(__file__).resolve().parent
RUNNER = HERE.parent / "scripts" / "runner.py"

failures: list[str] = []


def load_runner():
    """Import scripts/runner.py without running its __main__ block."""
    spec = importlib.util.spec_from_file_location("kegalign_runner", RUNNER)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def check(label, got, want):
    if got != want:
        print(f"not ok - {label}\n     got: {got}\nexpected: {want}")
        failures.append(label)
    else:
        print(f"ok - {label}")


def run_in(files, line):
    """Call the real estimate_chunk_size() in a temp dir holding `files`.

    files: {name: number_of_lines}. `line` is the text of each line, so the
    caller controls line_size, which the function infers from the first line of
    the first file it opens.
    """
    with tempfile.TemporaryDirectory() as tmp:
        for name, count in files.items():
            (pathlib.Path(tmp) / name).write_text(line * count)
        cwd = os.getcwd()
        try:
            os.chdir(tmp)
            return runner.estimate_chunk_size(argparse.Namespace(debug=False))
        finally:
            os.chdir(cwd)


runner = load_runner()
MAX_CHUNK_SIZE = 50000
LINE = "chr1\t1\t100\tchr2\t1\t100\t+\t1000\n"  # 27 bytes

check("no segments files at all", run_in({}, LINE), MAX_CHUNK_SIZE)
check("one segments file (the observed crash)", run_in({"a.segments": 100}, LINE), 100)
check("two segments files", run_in({"a.segments": 100, "b.segments": 200}, LINE), 150)
check(
    "eight segments files",
    run_in({f"{c}.segments": 10 * (i + 1) for i, c in enumerate("abcdefgh")}, LINE),
    67,
)
check(
    "a lone empty file yields MAX_CHUNK_SIZE, not a division by zero", run_in({"a.segments": 0}, LINE), MAX_CHUNK_SIZE
)
check("a one-character line is a valid divisor", run_in({"a.segments": 3}, "\n"), 3)
# An empty file must not poison line_size for the files that follow it. Ten good
# files plus one empty one used to give either the right answer or MAX_CHUNK_SIZE
# depending on os.scandir() order.
check(
    "one empty file among good ones does not poison the estimate",
    run_in({"empty.segments": 0, **{f"g{i}.segments": 100 for i in range(10)}}, LINE),
    100,
)
# 0 is diagonal_partition.py's documented "skip partitioning" value -- it prints
# the command unchanged and exits. Flooring the estimate at 1 would turn that
# no-op into one output file and one LASTZ command per segment line, for every
# file in the run. Reaching 0 through estimate_chunk_size() needs the quantile to
# fall below line_size, which depends on which file os.scandir() latches
# line_size from and so cannot be constructed deterministically. Assert the
# absence of the clamp instead, and pin the contract at the other end.
source = RUNNER.read_text()
check("no lower clamp on the estimate", "max(1, min(chunk_size" in source, False)
check(
    "diagonal_partition.py still treats 0 as skip-partitioning",
    "if chunk_size == 0:" in (RUNNER.parent / "diagonal_partition.py").read_text(),
    True,
)
check("non-.segments files are ignored", run_in({"a.txt": 100, "b.log": 200}, LINE), MAX_CHUNK_SIZE)
check(
    "result is clamped to MAX_CHUNK_SIZE",
    run_in({"a.segments": 10**6, "b.segments": 10**6}, "x\n"),
    MAX_CHUNK_SIZE,
)

if failures:
    print(f"\n{len(failures)} test(s) failed")
    sys.exit(1)
print("\nall tests passed")
