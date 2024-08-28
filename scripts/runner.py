#!/usr/bin/env python

import argparse
import collections
import concurrent.futures
import multiprocessing
import os
import queue
import re
import resource
import statistics
import subprocess
import sys
import time
import typing

SENTINEL_VALUE: typing.Final = "SENTINEL"
# CA_SENTEL_VALUE: typing.Final = ChunkAddress(0, 0, 0, 0, 0, SENTINEL_VALUE)
RUSAGE_ATTRS: typing.Final = ["ru_utime", "ru_stime", "ru_maxrss", "ru_minflt", "ru_majflt", "ru_inblock", "ru_oublock", "ru_nvcsw", "ru_nivcsw"]


class LastzCommands:
    def __init__(self) -> None:
        self.commands: dict[str, LastzCommand] = {}
        self.segalign_segments: SegAlignSegments = SegAlignSegments()

    def add(self, line: str) -> None:
        if line not in self.commands:
            self.commands[line] = LastzCommand(line)

        command = self.commands[line]
        self.segalign_segments.add(command.segments_filename)

    def segments(self) -> typing.Iterator["SegAlignSegment"]:
        for segment in self.segalign_segments:
            yield segment


class LastzCommand:
    lastz_command_regex = re.compile(r"lastz (.+?)?ref\.2bit\[nameparse=darkspace\]\[multiple\]\[subset=ref_block(\d+)\.name\] (.+?)?query\.2bit\[nameparse=darkspace\]\[subset=query_block(\d+)\.name] --format=(\S+) --ydrop=(\d+) --gappedthresh=(\d+) --strand=(minus|plus)(?: --ambiguous=(\S+))?(?: --(notrivial))?(?: --scores=(\S+))? --segments=tmp(\d+)\.block(\d+)\.r(\d+)\.(minus|plus)(?:\.split(\d+))?\.segments --output=tmp(\d+)\.block(\d+)\.r(\d+)\.(minus|plus)(?:\.split(\d+))?\.(\S+) 2> tmp(\d+)\.block(\d+)\.r(\d+)\.(minus|plus)(?:\.split(\d+))?\.err")

    def __init__(self, line: str) -> None:
        self.line = line
        self.args: list[str] = []
        self.target_filename: str = ''
        self.query_filename: str = ''
        self.data_folder: str = ''
        self.ref_block: int = 0
        self.query_block: int = 0
        self.output_format: str = ''
        self.ydrop: int = 0
        self.gappedthresh: int = 0
        self.strand: int | None = None
        self.ambiguous: bool = False
        self.nontrivial: bool = False
        self.scoring: str | None = None
        self.segments_filename: str = ''
        self.output_filename: str = ''
        self.error_filename: str = ''

        self._parse_command()

    def _parse_command(self) -> None:
        match = self.lastz_command_regex.match(self.line)
        if not match:
            sys.exit(f"Unkown lastz command format: {self.line}")

        self.data_folder = match.group(1)
        self.ref_block = int(match.group(2))
        self.query_block = int(match.group(4))
        self.output_format = match.group(5)
        self.ydrop = int(match.group(6))
        self.gappedthresh = int(match.group(7))
        strand = match.group(8)

        if strand == 'plus':
            self.strand = 0
        elif strand == 'minus':
            self.strand = 1

        self.target_filename = f"{self.data_folder}ref.2bit[nameparse=darkspace][multiple][subset=ref_block{self.ref_block}.name]"
        self.query_filename = f"{self.data_folder}query.2bit[nameparse=darkspace][subset=query_block{self.query_block}.name]"

        self.args = [
            "lastz",
            self.target_filename,
            self.query_filename,
            f"--format={self.output_format}",
            f"--ydrop={self.ydrop}",
            f"--gappedthresh={self.gappedthresh}",
            f"--strand={strand}"
        ]

        ambiguous = match.group(9)
        if ambiguous is not None:
            self.ambiguous = True
            self.args.append(f"--ambiguous={ambiguous}")

        nontrivial = match.group(10)
        if nontrivial is not None:
            self.nontrivial = True
            self.args.append("--nontrivial")

        scoring = match.group(11)
        if scoring is not None:
            self.scoring = scoring
            self.args.append(f"--scoring={scoring}")

        tmp_no = int(match.group(12))
        block_no = int(match.group(13))
        r_no = int(match.group(14))
        split = match.group(16)

        base_filename = f"tmp{tmp_no}.block{block_no}.r{r_no}.{strand}"

        if split is not None:
            base_filename = f"{base_filename}.split{split}"

        self.segments_filename = f"{base_filename}.segments"
        self.args.append(f"--segments={self.segments_filename}")

        self.output_filename = f"{base_filename}.{self.output_format}"
        self.args.append(f"--output={self.output_filename}")

        self.error_filename = f"{base_filename}.err"


class SegAlignSegments:
    def __init__(self) -> None:
        self.segments: dict[str, SegAlignSegment] = {}

    def add(self, filename: str) -> None:
        if filename not in self.segments:
            self.segments[filename] = SegAlignSegment(filename)

    def __iter__(self) -> "SegAlignSegments":
        return self

    def __next__(self) -> typing.Generator["SegAlignSegment", None, None]:
        for segment in sorted(self.segments.values()):
            yield segment


class SegAlignSegment:
    def __init__(self, filename: str):
        self.filename = filename
        self.tmp: int = 0
        self.block: int = 0
        self.r: int = 0
        self.strand: int = 0
        self.split: int | None = None
        self.fmt: str = ""
        self._parse_filename()

    def _parse_filename(self) -> None:
        match = re.match(r"tmp(\d+)\.block(\d+)\.r(\d+)\.(minus|plus)(?:\.split(\d+))?\.segments$", self.filename)
        if not match:
            sys.exit(f"Unkown segment filename format: {self.filename}")

        self.tmp = int(match.group(1))
        self.block = int(match.group(2))
        self.r = int(match.group(3))

        strand = match.group(4)
        if strand == 'plus':
            self.strand = 0
        if strand == 'minus':
            self.strand = 1

        split = match.group(5)
        if split is None:
            self.split = None
        else:
            self.split = int(split)

    def __lt__(self, other: "SegAlignSegment") -> bool:
        for attr in ['strand', 'tmp', 'block', 'r', 'split']:
            self_value = getattr(self, attr)
            other_value = getattr(other, attr)
            if self_value < other_value:
                return True
            elif self_value > other_value:
                return False

        return False


def main() -> None:
    args, segalign_args = parse_args()
    lastz_commands = LastzCommands()

    if args.diagonal_partition:
        num_diagonal_partitioners = args.num_cpu
    else:
        num_diagonal_partitioners = 0

    with multiprocessing.Manager() as manager:
        segalign_q: queue.Queue[str] = manager.Queue()
        skip_segalign = run_segalign(args, num_diagonal_partitioners, segalign_args, segalign_q, lastz_commands)

        if num_diagonal_partitioners > 0:
            diagonal_partition_q = segalign_q
            segalign_q = manager.Queue()
            run_diagonal_partitioners(args, num_diagonal_partitioners, diagonal_partition_q, segalign_q)

        segalign_q.put(SENTINEL_VALUE)
        output_q = segalign_q
        segalign_q = manager.Queue()

        output_filename = "lastz-commands.txt"
        if args.output_type == "commands":
            output_filename = args.output_file

        with open(output_filename, "w") as f:
            while True:
                line = output_q.get()
                if line == SENTINEL_VALUE:
                    output_q.task_done()
                    break

                # messy, fix this
                if skip_segalign:
                    lastz_commands.add(line)

                if args.output_type != "commands":
                    segalign_q.put(line)

                print(line, file=f)

        if args.output_type == "output":
            run_lastz(args, segalign_q, lastz_commands)

            with open(args.output_file, 'w') as of:
                print("##maf version=1", file=of)
                for lastz_command in lastz_commands.commands.values():
                    with open(lastz_command.output_filename) as f:
                        for line in f:
                            of.write(line)

        if args.output_type == "tarball":
            pass


def run_lastz(args: argparse.Namespace, input_q: queue.Queue[str], lastz_commands: LastzCommands) -> None:
    num_lastz_workers = args.num_cpu

    for _ in range(num_lastz_workers):
        input_q.put(SENTINEL_VALUE)

    if args.debug:
        r_beg = resource.getrusage(resource.RUSAGE_CHILDREN)
        beg: int = time.monotonic_ns()

    try:
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_lastz_workers) as executor:
            for i in range(num_lastz_workers):
                executor.submit(lastz_worker, input_q, i, lastz_commands)
    except Exception as e:
        sys.exit(f"Error: lastz failed: {e}")

    if args.debug:
        ns: int = time.monotonic_ns() - beg
        r_end = resource.getrusage(resource.RUSAGE_CHILDREN)
        print(f"lastz clock time: {ns} ns", file=sys.stderr, flush=True)
        for rusage_attr in RUSAGE_ATTRS:
            value = getattr(r_end, rusage_attr) - getattr(r_beg, rusage_attr)
            print(f"  lastz {rusage_attr}: {value}", file=sys.stderr, flush=True)


def lastz_worker(input_q: queue.Queue[str], instance: int, lastz_commands: LastzCommands) -> None:
    while True:
        line = input_q.get()
        if line == SENTINEL_VALUE:
            input_q.task_done()
            break

        if line not in lastz_commands.commands:
            sys.exit(f"Error: lastz worker {instance} unexpected command format: {line}")

        command = lastz_commands.commands[line]

        if not os.path.exists(command.output_filename):
            process = subprocess.run(command.args, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            for line in process.stdout.splitlines():
                print(line, file=sys.stdout, flush=True)

            if len(process.stderr) != 0:
                for line in process.stderr.splitlines():
                    print(line, file=sys.stderr, flush=True)

            if process.returncode != 0:
                sys.exit(f"Error: lastz {instance} exited with returncode {process.returncode}")


def run_diagonal_partitioners(args: argparse.Namespace, num_workers: int, input_q: queue.Queue[str], output_q: queue.Queue[str]) -> None:
    chunk_size = estimate_chunk_size(args)

    if args.debug:
        print(f"estimated chunk size: {chunk_size}", file=sys.stderr, flush=True)

    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        for i in range(num_workers):
            executor.submit(diagonal_partition_worker(args, input_q, output_q, chunk_size, i))


def diagonal_partition_worker(args: argparse.Namespace, input_q: queue.Queue[str], output_q: queue.Queue[str], chunk_size: int, instance: int) -> None:
    while True:
        line = input_q.get()
        if line == SENTINEL_VALUE:
            input_q.task_done()
            break

        run_args = ["python", f"{args.tool_directory}/diagonal_partition.py", str(chunk_size)]
        for word in line.split():
            run_args.append(word)
        process = subprocess.run(run_args, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=1, text=True)

        for line in process.stdout.splitlines():
            output_q.put(line)

        for line in process.stderr.splitlines():
            print(line, file=sys.stderr, flush=True)

        if process.returncode != 0:
            sys.exit(f"Error: diagonal partitioner {instance} exited with returncode {process.returncode}")


def estimate_chunk_size(args: argparse.Namespace) -> int:
    # only used when segment size is being estimated
    MAX_CHUNK_SIZE = 50000
    chunk_size = -1
    line_size = -1

    if args.debug:
        r_beg = resource.getrusage(resource.RUSAGE_SELF)
        beg: int = time.monotonic_ns()

    # get size of each segment assuming DELETE_AFTER_CHUNKING == True
    # takes into account already split segments
    fdict: typing.DefaultDict[str, int] = collections.defaultdict(int)
    for entry in os.scandir("."):
        if entry.name.endswith(".segments"):
            try:
                file_size = entry.stat().st_size
            except FileNotFoundError:
                continue

            if line_size == -1:
                try:
                    with open(entry.name) as f:
                        line_size = len(f.readline())  # add 1 for newline
                except FileNotFoundError:
                    continue

            fdict[entry.name.split(".split", 1)[0]] += file_size

    if len(fdict) < 7:
        # outliers can heavily skew prediction if <7 data points
        # to be safe, use 50% quantile
        chunk_size = int(statistics.quantiles(fdict.values())[1] // line_size)
    else:
        # otherwise use 75% quantile
        chunk_size = int(statistics.quantiles(fdict.values())[-1] // line_size)

    # if not enough data points, there is a chance of getting unlucky
    # minimize worst case by using MAX_CHUNK_SIZE

    chunk_size = min(chunk_size, MAX_CHUNK_SIZE)

    if args.debug:
        ns: int = time.monotonic_ns() - beg
        r_end = resource.getrusage(resource.RUSAGE_SELF)
        print(f"estimate chunk size clock time: {ns} ns", file=sys.stderr, flush=True)
        for rusage_attr in RUSAGE_ATTRS:
            value = getattr(r_end, rusage_attr) - getattr(r_beg, rusage_attr)
            print(f"  estimate chunk size {rusage_attr}: {value}", file=sys.stderr, flush=True)

    return chunk_size


def run_segalign(args: argparse.Namespace, num_sentinel: int, segalign_args: list[str], segalign_q: queue.Queue[str], commands: LastzCommands) -> bool:
    skip_segalign: bool = False

    # use the currently existing output file if it exists
    if args.debug:
        skip_segalign = load_segalign_output("lastz-commands.txt", segalign_q)

    if not skip_segalign:
        run_args = ["segalign"]
        run_args.extend(segalign_args)
        run_args.append("--num_threads")
        run_args.append(str(args.num_cpu))
        run_args.append("work/")

        if args.debug:
            beg: int = time.monotonic_ns()
            r_beg = resource.getrusage(resource.RUSAGE_CHILDREN)

        process = subprocess.run(run_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=1, text=True)

        for line in process.stdout.splitlines():
            commands.add(line)
            segalign_q.put(line)

        if len(process.stderr) != 0:
            for line in process.stderr.splitlines():
                print(line, file=sys.stderr, flush=True)

        if process.returncode != 0:
            sys.exit(f"Error: segalign exited with returncode {process.returncode}")

        if args.debug:
            ns: int = time.monotonic_ns() - beg
            r_end = resource.getrusage(resource.RUSAGE_CHILDREN)
            print(f"segalign clock time: {ns} ns", file=sys.stderr, flush=True)
            for rusage_attr in RUSAGE_ATTRS:
                value = getattr(r_end, rusage_attr) - getattr(r_beg, rusage_attr)
                print(f"  segalign {rusage_attr}: {value}", flush=True)

    for _ in range(num_sentinel):
        segalign_q.put(SENTINEL_VALUE)

    return skip_segalign


def load_segalign_output(filename: str, segalign_q: queue.Queue[str]) -> bool:
    load_success = False

    r_beg = resource.getrusage(resource.RUSAGE_SELF)
    beg: int = time.monotonic_ns()

    try:
        with open(filename) as f:
            for line in f:
                line = line.rstrip("\n")
                segalign_q.put(line)
        load_success = True
    except FileNotFoundError:
        pass

    if load_success:
        ns: int = time.monotonic_ns() - beg
        r_end = resource.getrusage(resource.RUSAGE_SELF)
        print(f"load output clock time: {ns} ns", file=sys.stderr, flush=True)
        for rusage_attr in RUSAGE_ATTRS:
            value = getattr(r_end, rusage_attr) - getattr(r_beg, rusage_attr)
            print(f"  load output {rusage_attr}: {value}", flush=True)

    return load_success


def parse_args() -> tuple[argparse.Namespace, list[str]]:
    parser = argparse.ArgumentParser(allow_abbrev=False)

    parser.add_argument("--output-type", nargs="?", const="commands", default="commands", type=str, choices=["commands", "output", "tarball"], help="output type (default: %(default)s)")
    parser.add_argument("--output-file", type=str, required=True, help="output pathname")
    parser.add_argument("--diagonal-partition", action="store_true", help="run diagonal partition optimization")
    parser.add_argument("--nogapped", action="store_true", help="don't perform gapped extension stage")
    parser.add_argument("--markend", action="store_true", help="write a marker line just before completion")
    parser.add_argument("--num-gpu", default=-1, type=int, help="number of GPUs to use (default: %(default)s [use all GPUs])")
    parser.add_argument("--num-cpu", default=-1, type=int, help="number of CPUs to use (default: %(default)s [use all CPUs])")
    parser.add_argument("--debug", action="store_true", help="print debug messages")
    parser.add_argument("--tool_directory", type=str, required=True, help="tool directory")

    if not sys.argv[1:]:
        parser.print_help()
        sys.exit(0)

    args, segalign_args = parser.parse_known_args(sys.argv[1:])

    cpus_available = len(os.sched_getaffinity(0))
    if args.num_cpu == -1:
        args.num_cpu = cpus_available
    elif args.num_cpu > cpus_available:
        sys.exit(f"Error: additional {args.num_cpu - cpus_available} CPUs")

    if args.nogapped:
        segalign_args.append("--nogapped")

    if args.markend:
        segalign_args.append("--markend")

    if args.num_gpu != -1:
        segalign_args.extend(["--num_gpu", f"{args.num_gpu}"])

    if args.debug:
        segalign_args.append("--debug")

    return args, segalign_args


if __name__ == "__main__":
    main()
