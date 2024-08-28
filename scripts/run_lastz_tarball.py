#!/usr/bin/env python

import argparse
import concurrent.futures
import json
import multiprocessing
import os
import queue
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
import typing


lastz_output_format_regex = re.compile(
    r"^(?:axt\+?|blastn|cigar|differences|general-?.+|lav|lav\+text|maf[-+]?|none|paf(?::wfmash)?|rdotplot|sam-?|softsam-?|text)$",
    re.IGNORECASE,
)


# Specifies the output format: lav, lav+text, axt, axt+, maf, maf+, maf-, sam, softsam, sam-, softsam-, cigar, BLASTN, PAF, PAF:wfmash, differences, rdotplot, text, general[:<fields>], or general-[:<fields>].
# ‑‑format=none can be used when no alignment output is desired.


def run_command(
    instance: int,
    input_queue: "queue.Queue[typing.Dict[str, typing.Any]]",
    output_queue: "queue.Queue[float]",
    debug: bool = False,
) -> str | None:
    os.chdir("galaxy/files")

    # These are not considered errors even though
    # we will end up with a segmented alignment
    truncation_regex = re.compile(
        r"truncating alignment (ending|starting) at \(\d+,\d+\);  anchor at \(\d+,\d+\)$"
    )
    truncation_msg = "truncation can be reduced by using --allocate:traceback to increase traceback memory"

    while True:
        command_dict = input_queue.get()

        if not command_dict:
            return None

        args = ["lastz", "--allocate:traceback=1.99G"]
        args.extend(command_dict["args"])

        stdin = command_dict["stdin"]
        if stdin is not None:
            stdin = open(stdin, "r")

        stdout = command_dict["stdout"]
        if stdout is not None:
            stdout = open(stdout, "w")

        stderr = command_dict["stderr"]
        if stderr is not None:
            stderr = open(stderr, "w")

        begin = time.perf_counter()
        p = subprocess.run(args, stdin=stdin, stdout=stdout, stderr=stderr)

        for var in [stdin, stdout, stderr]:
            if var is not None:
                var.close()

        # if there is a stderr_file, make sure it is
        # empty or only contains truncation messages
        stderr_ok = True
        stderr_file = command_dict["stderr"]

        if stderr_file is not None:
            try:
                statinfo = os.stat(stderr_file, follow_symlinks=False)
                if statinfo.st_size != 0:
                    with open(stderr_file) as f:
                        for stderr_line in f:
                            stderr_line = stderr_line.strip()
                            if (not truncation_regex.match(stderr_line) and stderr_line != truncation_msg):
                                stderr_ok = False
            except Exception:
                stderr_ok = False

        if p.returncode in [0, 1] and stderr_ok:
            elapsed = time.perf_counter() - begin
            output_queue.put(elapsed)
        else:
            return f"command failed (rc={p.returncode}): {' '.join(args)}"


class BatchTar:
    def __init__(self, pathname: str, debug: bool = False) -> None:
        self.pathname = pathname
        self.debug = debug
        self.tarfile = None
        self.commands: typing.List[typing.Dict[str, typing.Any]] = []
        self.format_name = "tabular"
        self._extract()
        self._load_commands()
        self._load_format()

    def batch_commands(self) -> typing.Iterator[typing.Dict[str, typing.Any]]:
        for command in self.commands:
            yield command

    def final_output_format(self) -> str:
        return self.format_name

    def _extract(self) -> None:
        try:
            self.tarball = tarfile.open(
                name=self.pathname, mode="r:*", format=tarfile.GNU_FORMAT
            )
        except FileNotFoundError:
            sys.exit(f"ERROR: unable to find input tarball: {self.pathname}")
        except tarfile.ReadError:
            sys.exit(f"ERROR: error reading input tarball: {self.pathname}")

        begin = time.perf_counter()
        self.tarball.extractall(filter="data")
        self.tarball.close()
        elapsed = time.perf_counter() - begin

        if self.debug:
            print(
                f"Extracted tarball in {elapsed} seconds", file=sys.stderr, flush=True
            )

    def _load_commands(self) -> None:
        try:
            f = open("galaxy/commands.json")
        except FileNotFoundError:
            sys.exit(
                f"ERROR: input tarball missing galaxy/commands.json: {self.pathname}"
            )

        begin = time.perf_counter()
        for json_line in f:
            json_line = json_line.rstrip("\n")
            try:
                command_dict = json.loads(json_line)
            except json.JSONDecodeError:
                sys.exit(
                    f"ERROR: bad json line in galaxy/commands.json: {self.pathname}"
                )

            self._load_command(command_dict)

        f.close()
        elapsed = time.perf_counter() - begin

        if self.debug:
            print(
                f"loaded {len(self.commands)} commands in {elapsed} seconds ",
                file=sys.stderr,
                flush=True,
            )

    def _load_command(self, command_dict: typing.Dict[str, typing.Any]) -> None:
        # check command_dict structure
        field_types: typing.Dict[str, typing.List[typing.Any]] = {
            "executable": [str],
            "args": [list],
            "stdin": [str, "None"],
            "stdout": [str, "None"],
            "stderr": [str, "None"],
        }

        bad_format = False
        for field_name in field_types.keys():
            # missing field
            if field_name not in command_dict:
                bad_format = True
                break

            # incorrect field type
            good_type = False
            for field_type in field_types[field_name]:
                if isinstance(field_type, str) and field_type == "None":
                    if command_dict[field_name] is None:
                        good_type = True
                        break
                elif isinstance(command_dict[field_name], field_type):
                    good_type = True
                    break

            if good_type is False:
                bad_format = True

        if not bad_format:
            # all args must be strings
            for arg in command_dict["args"]:
                if not isinstance(arg, str):
                    bad_format = True
                    break

        if bad_format:
            sys.exit(
                f"ERROR: unexpected json format in line in galaxy/commands.json: {self.pathname}"
            )

        self.commands.append(command_dict)

    def _load_format(self) -> None:
        try:
            with open("galaxy/format.txt") as f:
                format_name = f.readline()
                format_name = format_name.rstrip("\n")
        except FileNotFoundError:
            sys.exit(f"ERROR: input tarball missing galaxy/format.txt: {self.pathname}")

        if format_name in ["bam", "maf"]:
            self.format_name = format_name
        elif format_name == "differences":
            self.format_name = "interval"


class TarRunner:
    def __init__(
        self,
        input_pathname: str,
        output_pathname: str,
        parallel: int,
        debug: bool = False,
    ) -> None:
        self.input_pathname = input_pathname
        self.output_pathname = output_pathname
        self.parallel = parallel
        self.debug = debug
        self.batch_tar = BatchTar(self.input_pathname, debug=self.debug)
        self.output_file_format: typing.Dict[str, str] = {}
        self.output_files: typing.Dict[str, typing.List[str]] = {}
        self._set_output()
        self._set_target_query()

    def _set_output(self) -> None:
        for command_dict in self.batch_tar.batch_commands():
            output_file = None
            output_format = None

            for arg in command_dict["args"]:
                if arg.startswith("--format="):
                    output_format = arg[9:]
                elif arg.startswith("--output="):
                    output_file = arg[9:]

            if output_file is None:
                f = tempfile.NamedTemporaryFile(dir="galaxy/files", delete=False)
                output_file = os.path.basename(f.name)
                f.close()
                command_dict["args"].append(f"--output={output_file}")

            if output_format is None:
                output_format = "lav"
                command_dict["args"].append(f"--format={output_format}")

            if not lastz_output_format_regex.match(output_format):
                sys.exit(f"ERROR: invalid output format: {output_format}")

            self.output_file_format[output_file] = output_format

        for output_file, output_format in self.output_file_format.items():
            self.output_files.setdefault(output_format, [])
            self.output_files[output_format].append(output_file)

    def _set_target_query(self) -> None:
        for command_dict in self.batch_tar.batch_commands():
            new_args: typing.List[str] = []

            for arg in command_dict["args"]:
                if arg.startswith("--target="):
                    new_args.insert(0, arg[9:])
                elif arg.startswith("--query="):
                    new_args.insert(1, arg[8:])
                else:
                    new_args.append(arg)

            command_dict["args"] = new_args

    def run(self) -> None:
        run_times = []
        begin = time.perf_counter()

        with multiprocessing.Manager() as manager:
            input_queue: queue.Queue[typing.Dict[str, typing.Any]] = manager.Queue()
            output_queue: queue.Queue[float] = manager.Queue()

            for command_dict in self.batch_tar.batch_commands():
                input_queue.put(command_dict)

            # use the empty dict as a sentinel
            for _ in range(self.parallel):
                input_queue.put({})

            with concurrent.futures.ProcessPoolExecutor(
                max_workers=self.parallel
            ) as executor:
                futures = [
                    executor.submit(
                        run_command,
                        instance,
                        input_queue,
                        output_queue,
                        debug=self.debug,
                    )
                    for instance in range(self.parallel)
                ]

            found_falures = False

            for f in concurrent.futures.as_completed(futures):
                result = f.result()
                if result is not None:
                    print(f"lastz: {result}", file=sys.stderr, flush=True)

                if not f.done() or f.cancelled() or f.exception() is not None:
                    found_falures = True

            while not output_queue.empty():
                run_time = output_queue.get()
                run_times.append(run_time)

            if found_falures:
                sys.exit("lastz command failed")

        elapsed = time.perf_counter() - begin

        if self.debug:
            print(f"elapsed {elapsed}", file=sys.stderr, flush=True)

        self._cleanup()

    def _cleanup(self) -> None:
        num_output_files = len(self.output_files.keys())
        if num_output_files != 1:
            sys.exit(f"ERROR: expecting a single output file, found {num_output_files}")

        final_output_format = self.batch_tar.final_output_format()

        for file_type, file_list in self.output_files.items():
            with open(f"output.{final_output_format}", "w") as ofh:
                if final_output_format == "maf":
                    print("##maf version=1", file=ofh)

                for filename in file_list:
                    with open(f"galaxy/files/{filename}") as ifh:
                        for line in ifh:
                            ofh.write(line)

        src_filename = f"output.{final_output_format}"
        shutil.copy2(src_filename, self.output_pathname)

        output_metadata = {
            "output": {
                "ext": final_output_format,
            }
        }

        with open("galaxy.json", "w") as ofh:
            json.dump(output_metadata, ofh)


def main() -> None:
    if not hasattr(tarfile, "data_filter"):
        sys.exit("ERROR: extracting may be unsafe; consider updating Python")

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, required=True)
    parser.add_argument("--output", type=str, required=True)
    parser.add_argument("--parallel", type=int, default=2, required=False)
    parser.add_argument("--debug", action="store_true", required=False)

    args = parser.parse_args()
    runner = TarRunner(args.input, args.output, args.parallel, args.debug)
    runner.run()


if __name__ == "__main__":
    main()
