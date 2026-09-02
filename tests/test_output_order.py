#!/usr/bin/env python3
"""MAF output order must not depend on which partitioner finished first.

Run with:  python3 tests/test_output_order.py

The shell implementation KegAlign replaced ran lastz in parallel and then
concatenated deterministically, plus before minus, version-sorted within each:

    for i in tmp*.plus.*;  do echo $i; done | sort -V
    for i in tmp*.minus.*; do echo $i; done | sort -V

The Python port dropped that: main() iterated lastz_commands.commands.values(),
which is insertion order -- whatever order the diagonal partitioners emitted.
That was masked for as long as the partitioners accidentally ran serially in the
parent process; once they actually run in parallel, the same input produces MAF
blocks in a different order run to run.

KegAlignSegment.__lt__ already encoded the intended order and was never called.
"""

import importlib.util
import pathlib
import random
import sys

RUNNER = pathlib.Path(__file__).resolve().parent.parent / "scripts" / "runner.py"
failures: list[str] = []


def load_runner():
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


runner = load_runner()

TEMPLATE = (
    "lastz ref.2bit[nameparse=darkspace][multiple][subset=ref_block{r}.name] "
    "query.2bit[nameparse=darkspace][subset=query_block{b}.name] --format=maf "
    "--ydrop=9430 --gappedthresh=3000 --strand={s} "
    "--segments=tmp{t}.block{b}.r{r}.{s}.segments "
    "--output=tmp{t}.block{b}.r{r}.{s}.maf 2> tmp{t}.block{b}.r{r}.{s}.err"
)

lines = [
    TEMPLATE.format(t=t, b=b, r=r, s=s) for s in ("plus", "minus") for t in (0, 1, 2) for b in (0, 1) for r in (0,)
]


def order_of(shuffled):
    commands = runner.LastzCommands()
    for line in shuffled:
        commands.add(line)
    return [c.output_filename for c in commands.sorted_commands()]


expected = order_of(lines)

# plus must come before minus, matching SegAlign's two loops
first_minus = next(i for i, f in enumerate(expected) if ".minus." in f)
last_plus = max(i for i, f in enumerate(expected) if ".plus." in f)
check("all plus files precede all minus files", last_plus < first_minus, True)
check("first file is tmp0.block0.r0.plus", expected[0], "tmp0.block0.r0.plus.maf")

# and the order must not depend on the order the commands arrived in
rng = random.Random(11)
stable = True
for _ in range(25):
    shuffled = lines[:]
    rng.shuffle(shuffled)
    if order_of(shuffled) != expected:
        stable = False
        break
check("order is identical across 25 shuffled arrival orders", stable, True)

# the unsorted path is what regressed; show it actually differs
unsorted = list(runner.LastzCommands().commands.values())
shuffled = lines[:]
rng.shuffle(shuffled)
c = runner.LastzCommands()
for line in shuffled:
    c.add(line)
insertion = [x.output_filename for x in c.commands.values()]
check("insertion order really is different (so the sort is load-bearing)", insertion != expected, True)

# The checks above exercise sorted_commands() directly, so they stay green if the
# method is correct but main() never calls it -- which is exactly the regression.
# Assert the call site too. (Mutation-checked: reverting either the sort key or
# the call site turns this file red.)
source = RUNNER.read_text()
check(
    "main() concatenates via sorted_commands(), not commands.values()",
    "for lastz_command in lastz_commands.sorted_commands():" in source,
    True,
)
check(
    "the unsorted concatenation is gone",
    "for lastz_command in lastz_commands.commands.values():" in source,
    False,
)

if failures:
    print(f"\n{len(failures)} test(s) failed")
    sys.exit(1)
print("\nall tests passed")
