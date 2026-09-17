#!/usr/bin/env python3
"""A failed lastz command must not leave run_lastz_tarball.py exiting 0.

Run with:  python3 tests/test_lastz_failure_propagation.py

Two independent defects let a failed run report success.

The first was the predicate over future state:

    for f in concurrent.futures.as_completed(futures):
        if not f.done() or f.cancelled() or f.exception() is not None:
            found_falures = True

as_completed only ever yields done futures, so the first clause is dead.  A
worker that ran to completion and *returned* the commands that failed -- the
ordinary case, a lastz exiting nonzero -- is cancelled=False and
exception()=None, indistinguishable from a healthy one.  The failures were
printed to stderr and then discarded, and the tool exited 0 with a short MAF.
No predicate over future *state* alone can see that case; the returned value
has to be read.  Hence collect_failures().

The second was the return code: rc 1 was accepted unconditionally, because
lastz exits 1 with a tolerable warning about truncation.  But that warning is
recognised by reading the stderr file, and when a command has no stderr file
there is nothing to recognise it by -- so rc 1 was accepted with no evidence
at all.  Hence command_succeeded().

These tests call both functions.  A test that greps the source for the fixed
text would pass against any implementation that merely mentions it.
"""

import concurrent.futures
import importlib.util
import pathlib
import sys

SCRIPT = pathlib.Path(__file__).resolve().parent.parent / "scripts" / "run_lastz_tarball.py"
failures: list[str] = []


def load_script():
    # a .pyc left in scripts/__pycache__ by an earlier py_compile shadowed the
    # source here once, and the test silently graded stale bytecode
    importlib.invalidate_caches()
    spec = importlib.util.spec_from_file_location("kegalign_run_lastz_tarball", SCRIPT)
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


def returned(value):
    """A future that ran to completion and returned `value`."""
    future: concurrent.futures.Future = concurrent.futures.Future()
    future.set_result(value)
    return future


def raised(exception):
    future: concurrent.futures.Future = concurrent.futures.Future()
    future.set_exception(exception)
    return future


def cancelled():
    future: concurrent.futures.Future = concurrent.futures.Future()
    assert future.cancel()
    return future


module = load_script()
collect_failures = module.collect_failures
command_succeeded = module.command_succeeded

# --- collect_failures: the three ways a worker can fail ------------------

# This is the case the old state-only predicate could not see, and the only
# one that happens in a normal run.
check(
    "a worker that RETURNS failures is a failure",
    collect_failures([returned(["command failed (rc=1): lastz a b"])]),
    ["command failed (rc=1): lastz a b"],
)
check(
    "a worker that RAISES is a failure",
    collect_failures([raised(RuntimeError("boom"))]),
    ["worker raised: boom"],
)
check(
    "a CANCELLED worker is a failure",
    collect_failures([cancelled()]),
    ["worker was cancelled"],
)
check(
    "a clean run reports nothing",
    collect_failures([returned([]), returned([])]),
    [],
)
check(
    "every failure is collected, not just the first",
    len(collect_failures([returned(["a", "b"]), raised(ValueError("c")), returned([])])),
    3,
)
# cancelled() must be tested before exception(), which raises on a cancelled
# future.  If the order is wrong this call raises instead of returning.
check(
    "a cancelled future among healthy ones does not raise",
    len(collect_failures([returned([]), cancelled(), returned(["x"])])),
    2,
)

# --- command_succeeded: rc 1 needs evidence ------------------------------

check(
    "rc 0 with no stderr file succeeds",
    command_succeeded(0, None, True),
    True,
)
check(
    "rc 1 with NO stderr file fails -- nothing excuses it",
    command_succeeded(1, None, True),
    False,
)
check(
    "rc 1 with a clean stderr file succeeds (tolerated truncation)",
    command_succeeded(1, "lastz.err", True),
    True,
)
check(
    "rc 1 with an unrecognised stderr file fails",
    command_succeeded(1, "lastz.err", False),
    False,
)
check(
    "rc 2 always fails",
    command_succeeded(2, "lastz.err", True),
    False,
)

if failures:
    print(f"\n{len(failures)} test(s) failed")
    sys.exit(1)
print("\nall tests passed")
