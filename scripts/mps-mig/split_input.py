#!/usr/bin/env python3

import argparse
import collections
import concurrent.futures
import gzip
import heapq
import math
import os
import re
import resource
import subprocess
import sys
import time
import typing

TEN_MB: typing.Final = 10_000_000
TEN_KB: typing.Final = 10_000

RUSAGE_ATTRS: typing.Final = [
    "ru_utime",
    "ru_stime",
    "ru_maxrss",
    "ru_minflt",
    "ru_majflt",
    "ru_inblock",
    "ru_oublock",
    "ru_nvcsw",
    "ru_nivcsw",
]

FastaSequence = collections.namedtuple(
    "FastaSequence", ["description", "sequence", "length"], defaults=["", "", 0]
)


def debug_start(
    who: int = resource.RUSAGE_SELF, message: str = ""
) -> tuple[resource.struct_rusage, int, int]:
    print(f"DEBUG: {message}", file=sys.stderr, flush=True)
    r_beg = resource.getrusage(who)
    beg = time.monotonic_ns()
    return r_beg, beg, who


def debug_end(
    r_beg: resource.struct_rusage, beg: int, who: int, message: str = ""
) -> None:
    ns = time.monotonic_ns() - beg
    r_end = resource.getrusage(who)
    print(f"DEBUG: {message}: {ns} ns", file=sys.stderr, flush=True)
    for rusage_attr in RUSAGE_ATTRS:
        value = getattr(r_end, rusage_attr) - getattr(r_beg, rusage_attr)
        print(f"DEBUG:   {rusage_attr}: {value}", file=sys.stderr, flush=True)


class FastaFile:
    def __init__(self, pathname: str, debug: bool = False) -> None:
        self.pathname = pathname
        self.sequences: list[FastaSequence] = []
        self._read_fasta()
        self.sequences.sort(key=lambda x: x.length, reverse=True)

    def _read_fasta(self) -> None:
        description = ""
        seqs: list[str] = []

        if args.debug:
            debug_r_beg, debug_beg, debug_who = debug_start(
                resource.RUSAGE_SELF, f"loading fasta {self.pathname}"
            )

        with self._get_open_method() as f:
            for line in f:
                line = line.rstrip()

                if line.startswith(">"):
                    if seqs:
                        sequence = "".join(seqs)
                        self.sequences.append(
                            FastaSequence(description, sequence, len(sequence))
                        )
                        seqs.clear()

                    description = line
                else:
                    seqs.append(line)

            if seqs:
                sequence = "".join(seqs)
                self.sequences.append(
                    FastaSequence(description, sequence, len(sequence))
                )

        if args.debug:
            debug_end(
                debug_r_beg,
                debug_beg,
                debug_who,
                f"loaded fasta {self.pathname}",
            )

    def _get_open_method(self) -> typing.TextIO:
        try:
            with open(self.pathname, "rb") as f:
                if f.read(2) == b"\x1f\x8b":
                    return gzip.open(self.pathname, mode="rt")
        except FileNotFoundError:
            sys.exit(f"ERROR: Unable to read file: {self.pathname}")
        except Exception:
            pass

        return open(self.pathname, mode="rt")

    @property
    def total_bases(self) -> int:
        total = 0
        for fasta_sequence in self.sequences:
            total += fasta_sequence.length

        return total

    def __iter__(self) -> typing.Iterator[FastaSequence]:
        for fasta_sequence in self.sequences:
            yield fasta_sequence

    def to_single_seq(self) -> None:
        description = self.sequences[0].description
        sequence = "".join([seq.sequence for seq in self.sequences])

        self.sequences.clear()
        self.sequences.append(FastaSequence(description, sequence, len(sequence)))

    def discard_sequences_after_and_including(
        self, description: str, debug: bool = False
    ) -> None:
        split_index = -1
        for idx, sequence in enumerate(self.sequences):
            if sequence.description == f">{description}":
                split_index = idx
                break

        if split_index == -1:
            sys.exit(f"ERROR: sequence {description} not found")

        if debug:
            print(
                f"DEBUG: discarding sequences after and including {description}",
                file=sys.stdout,
                flush=True,
            )

        if split_index == 0:
            self.sequences.clear()
        else:
            self.sequences = self.sequences[: split_index - 1]


def chunk_file(
    target_fasta: FastaFile,
    output_dir: str,
    end_chr: str = "",
    save_split: bool = False,
    single_chr: bool = False,
    chunk_size: int = TEN_MB,
    overlap_size: int = TEN_KB,
    debug: bool = False,
) -> None:
    sequence_list = [seq for seq in target_fasta]

    title = sequence_list[0].description
    assert ">" in title

    if single_chr:
        target_fasta.to_single_seq()

    if end_chr:
        target_fasta.discard_sequences_after_and_including(end_chr, debug)

        if save_split:
            save_split_file = f"{target_fasta.pathname}.to_{end_chr}.fa"

            if debug:
                print(
                    f"Saving split file {target_fasta.pathname} to {save_split_file}",
                    file=sys.stdout,
                    flush=True,
                )

            with open(save_split_file, "w") as f:
                for seq in target_fasta:
                    print(f"{seq.description}", file=f)
                    print(f"{seq.sequence}", file=f)


#    # look into this re: overlap
#    for line_index in range(0, len(data), chunk_size - overlap_size):
#        end_line_index = min(len(data), line_index + chunk_size)
#        block_file_name = os.path.join(output_dir, f"block_part{line_index}-{end_line_index}")
#
#        with open(block_file_name, "w") as out_file:
#            to_write = data[line_index:end_line_index]
#
#            if ">" not in to_write[0]:
#                out_file.write(title)
#            out_file.writelines(to_write)
#
#            # find last sequence name for next file
#            for line in reversed(to_write[:-overlap_size]):
#                if ">" in line:
#                    title = line
#                    break


def convert_to_2bit(root_dir: str) -> None:
    commands = []

    for dirpath, dirnames, filenames in os.walk(root_dir):
        for filename in filenames:
            match = re.match(r"chunk_\d+$", filename)
            if match:
                pathname = os.path.join(root_dir, filename)
                new_pathname = f"{pathname}.2bit"
                command = f"faToTwoBit {pathname} {new_pathname}"
                commands.append(command)

    cpus_available = len(os.sched_getaffinity(0))
    num_commands = len(commands)
    if cpus_available > num_commands:
        cpus_available = num_commands

    if args.debug:
        print(
            f"DEBUG: converting to 2bit {cpus_available} CPUs",
            file=sys.stderr,
            flush=True,
        )

    with concurrent.futures.ProcessPoolExecutor(max_workers=cpus_available) as executor:
        for _ in executor.map(twobit_wrapper, commands):
            pass


def twobit_wrapper(command: str) -> None:
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True
    )
    stdout, stderr = process.communicate()
    if stdout:
        print(f"{stdout}", file=sys.stdout)
    if stderr:
        print(f"{stderr}", file=sys.stderr)


def split_chr(
    target_fasta: FastaFile,
    output_dir: str,
    num_chunks: int,
    write_to_output_dir: bool = True,
    debug: bool = False,
) -> list[int]:
    # Longest-processing-time-first Algorithm
    # sequence length is used as a proxy for processing-time
    chunk_size_list = []

    pq: list[tuple[int, int]] = []
    for i in range(num_chunks):
        heapq.heappush(pq, (0, i))  # bin size and bin id

    files: dict[int, list[int]] = {}
    for i in range(num_chunks):
        files[i] = []

    i = 0
    for sequence in target_fasta:
        # get smallest file
        size, bin_no = heapq.heappop(pq)
        size += sequence.length
        heapq.heappush(pq, (size, bin_no))
        files[bin_no].append(i)
        i += 1

    seen_inds = []  # for sanity checking

    for bin_no, chr_indexes in files.items():
        bin_size = 0

        for ind in chr_indexes:
            assert ind not in seen_inds  # sanity check
            seen_inds.append(ind)
            bin_size += target_fasta.sequences[ind].length

        if debug:
            print(
                f"DEBUG: chunk_{bin_no} num bp {bin_size}", file=sys.stderr, flush=True
            )

        chunk_size_list.append(bin_size)
        block_file_name = os.path.join(output_dir, f"chunk_{bin_no}")

        if write_to_output_dir:
            with open(block_file_name, "w") as out_file:
                for ind in chr_indexes:
                    sequence = target_fasta.sequences[ind]
                    print(f"{sequence.description}", file=out_file)
                    print(f"{sequence.sequence}", file=out_file)

    if debug:
        print(
            f"packed {len(target_fasta.sequences)} sequences into {num_chunks} bins",
            file=sys.stderr,
            flush=True,
        )

    assert len(seen_inds) == len(target_fasta.sequences)
    assert len(chunk_size_list) == num_chunks
    return chunk_size_list


def parallel_wrapper(pass_list: tuple[FastaFile, str, int, bool, bool]) -> list[int]:
    return split_chr(*pass_list)


def mse(data: list[int], base: int) -> float:
    mse = [(base - i) ** 2 for i in data]
    return sum(mse) / len(mse)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Input sequence in fasta or fasta.gz format",
    )
    parser.add_argument("--out", type=str, required=True, help="Output directory")

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--end_chr", default="", type=str, help="Use data up until this sequence"
    )
    group.add_argument(
        "--single_chr",
        action="store_true",
        help="Use only first sequence name for all sequences",
    )

    parser.add_argument(
        "--save_split_inputs", action="store_true", help="Save split input files"
    )
    parser.add_argument(
        "--to_2bit", action="store_true", help="Convert data into .2bit format"
    )
    parser.add_argument(
        "--max_chunks",
        default=0,
        type=int,
        help="Calculate chunk size based on input size",
    )
    parser.add_argument(
        "--split_chr",
        action="store_false",
        help="Make each sequence a separate partition",
    )
    parser.add_argument(
        "--goal_bp",
        default=0,
        type=int,
        help="goal basepairs across chunks, calculate best using MSE. Check up to --max_chunks number of bins",
    )
    parser.add_argument("--debug", action="store_true", help="Print debug information")

    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    target_file = args.input
    target_fasta = FastaFile(target_file, debug=True)

    target_block_dir = args.out

    target_chunk_size = TEN_MB
    query_chunk_size = TEN_MB
    overlap_size = TEN_KB

    os.makedirs(target_block_dir, exist_ok=True)

    if args.max_chunks and not args.split_chr:
        total_bases = target_fasta.total_bases
        target_chunk_size = -(total_bases // -args.max_chunks)

        if args.debug:
            print(f"total_bases = {total_bases}", file=sys.stderr, flush=True)
            print(
                f"target chunk size = {target_chunk_size}", file=sys.stderr, flush=True
            )

    chr_end_target = args.end_chr.strip()

    if args.split_chr:
        bin_count = args.max_chunks
        goal_bp = args.goal_bp
        best_bin_count = -1
        best_bin_loss = math.inf
        prev_avg = 0

        if args.goal_bp:
            pass_list = [
                (target_fasta, "", i, False, False) for i in range(1, bin_count + 1)
            ]

            cpus_available = len(os.sched_getaffinity(0))
            if cpus_available > bin_count:
                cpus_available = bin_count

            if args.debug:
                print(
                    f"DEBUG: spliting using {cpus_available} CPUs",
                    file=sys.stderr,
                    flush=True,
                )

            with concurrent.futures.ProcessPoolExecutor(
                max_workers=cpus_available
            ) as executor:
                for bins in executor.map(parallel_wrapper, pass_list):
                    i = len(bins)
                    loss = mse(bins, goal_bp)

                    if args.debug:
                        print(
                            f"DEBUG: * bin count {i}, mse {int(loss)}, bins {bins}",
                            file=sys.stderr,
                            flush=True,
                        )

                    if loss < best_bin_loss:
                        best_bin_count = i
                        best_bin_loss = loss
        else:
            for i in range(1, bin_count + 1):
                bins = split_chr(
                    target_fasta,
                    target_block_dir,
                    num_chunks=i,
                    write_to_output_dir=False,
                    debug=args.debug,
                )
                avg = sum(bins) / len(bins)
                loss = abs(goal_bp - avg)

                if args.debug:
                    print(
                        f"DEBUG: bin count {i}, avg difference {int(loss)}",
                        file=sys.stderr,
                        flush=True,
                    )

                if loss < best_bin_loss:
                    best_bin_count = i
                    best_bin_loss = loss

        bin_count = best_bin_count

        if args.debug:
            print(
                f"DEBUG: bin_count = {bin_count}, loss={best_bin_loss}",
                file=sys.stderr,
                flush=True,
            )

        split_chr(target_fasta, target_block_dir, bin_count, debug=args.debug)
    else:
        chunk_file(
            target_fasta,
            target_block_dir,
            end_chr=chr_end_target,
            save_split=args.save_split_inputs,
            single_chr=args.single_chr,
            chunk_size=target_chunk_size,
            overlap_size=overlap_size,
            debug=args.debug,
        )

    if args.to_2bit:
        # Convert to 2bit format
        convert_to_2bit(target_block_dir)

"""
python split_input.py --input /home/mdl/abg6029/WGA_tests/inputs/hg38.fa --out /home/mdl/abg6029/WGA_tests/inputs/blocked_20_chr_hg38/ --to_2bit true --split_chr true --max_chunks 20

python split_input.py --input hg38.fa --out ./blocked_20_chr_hg38 --to_2bit true --goal_bp 200000000 --max_chunks 20

"""
