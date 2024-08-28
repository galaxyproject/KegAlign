#!/usr/bin/env python

import argparse
import configparser
import json
import os
import resource
import sys
import tarfile
import time
import typing

import bashlex

RUSAGE_ATTRS: typing.Final = ["ru_utime", "ru_stime", "ru_maxrss", "ru_minflt", "ru_majflt", "ru_inblock", "ru_oublock", "ru_nvcsw", "ru_nivcsw"]


class PackageFile:
    def __init__(
        self,
        pathname: str = "data_package.tgz",
        top_dir: str = "galaxy",
        data_dir: str = "files",
        config_file: str = "commands.json",
        format_file: str = "format.txt"
    ) -> None:
        self.pathname: str = os.path.realpath(pathname)
        self.data_root: str = os.path.join(top_dir, data_dir)
        self.config_path: str = os.path.join(top_dir, config_file)
        self.config_file: str = config_file
        self.format_path: str = os.path.join(top_dir, format_file)
        self.format_file: str = format_file
        self.tarfile: typing.Optional[tarfile.TarFile] = None
        self.name_cache: typing.Dict[typing.Any, typing.Any] = {}
        self.working_dir: str = os.path.realpath(os.getcwd())

    def _initialize(self) -> None:
        if self.tarfile is None:
            self.tarfile = tarfile.open(
                name=self.pathname,
                mode="w:gz",
                format=tarfile.GNU_FORMAT,
                compresslevel=6,
            )

    def add_config(self, pathname: str) -> None:
        if self.tarfile is None:
            self._initialize()

        source_path = os.path.realpath(pathname)

        if self.tarfile is not None:
            self.tarfile.add(source_path, arcname=self.config_path, recursive=False)

    def add_file(self, pathname: str, arcname: typing.Optional[str] = None) -> None:
        if self.tarfile is None:
            self._initialize()

        source_path = os.path.realpath(pathname)

        dest_path = None

        if arcname is None:
            dest_path = os.path.join(self.data_root, os.path.basename(source_path))
        else:
            arc_path = os.path.realpath(arcname)
            rel_path = os.path.relpath(arc_path, self.working_dir)
            if not (os.path.isabs(rel_path) or rel_path.startswith("../")):
                dest_path = os.path.join(self.data_root, rel_path)
            else:
                sys.exit("path fail")

        if dest_path is not None:
            if self.tarfile is not None:
                if dest_path not in self.name_cache:
                    try:
                        self.tarfile.add(
                            source_path, arcname=dest_path, recursive=False
                        )
                    except FileNotFoundError:
                        sys.exit(f"missing source file {source_path}")

                    self.name_cache[dest_path] = 1
                    # print(f"added: {dest_path}", flush=True)

    def add_format(self, pathname: str) -> None:
        if self.tarfile is None:
            self._initialize()

        source_path = os.path.realpath(pathname)

        if self.tarfile is not None:
            self.tarfile.add(source_path, arcname=self.format_path, recursive=False)

    def close(self) -> None:
        if self.tarfile is not None:
            self.tarfile.close()
            self.tarfile = None


class bashCommandLineFile:
    def __init__(
        self,
        pathname: str,
        config: configparser.ConfigParser,
        args: argparse.Namespace,
        package_file: PackageFile,
    ) -> None:
        self.pathname: str = pathname
        self.config = config
        self.args = args
        self.package_file = package_file
        self.executable: typing.Optional[str] = None
        self._parse_lines()
        self._write_format()

    def _parse_lines(self) -> None:
        with open("commands.json", "w") as ofh:
            with open(self.pathname) as f:
                line: str
                for line in f:
                    line = line.rstrip("\n")
                    command_dict = self._parse_line(line)
                    # we may want to re-write args here
                    new_args_list = []

                    args_list = command_dict.get("args", [])
                    for arg in args_list:
                        if arg.startswith("--target="):
                            pathname = arg[9:]
                            new_args_list.append(arg)
                            if "[" in pathname:
                                elems = pathname.split("[")
                                sequence_file = elems.pop(0)
                                self.package_file.add_file(sequence_file, sequence_file)
                                for elem in elems:
                                    if elem.endswith("]"):
                                        elem = elem[:-1]
                                        if elem.startswith("subset="):
                                            subset_file = elem[7:]
                                            self.package_file.add_file(subset_file)

                        elif arg.startswith("--query="):
                            pathname = arg[8:]
                            new_args_list.append(arg)
                            if "[" in pathname:
                                elems = pathname.split("[")
                                sequence_file = elems.pop(0)
                                self.package_file.add_file(sequence_file, sequence_file)
                                for elem in elems:
                                    if elem.endswith("]"):
                                        elem = elem[:-1]
                                        if elem.startswith("subset="):
                                            subset_file = elem[7:]
                                            self.package_file.add_file(subset_file)
                        elif arg.startswith("--segments="):
                            pathname = arg[11:]
                            new_args_list.append(arg)
                            self.package_file.add_file(pathname)
                        elif arg.startswith("--scores="):
                            pathname = arg[9:]
                            new_args_list.append("--scores=data/scores.txt")
                            self.package_file.add_file(pathname, "data/scores.txt")
                        else:
                            new_args_list.append(arg)

                    command_dict["args"] = new_args_list
                    print(json.dumps(command_dict), file=ofh)

        self.package_file.add_config("commands.json")

    def _parse_line(self, line: str) -> typing.Dict[str, typing.Any]:
        # resolve shell redirects
        trees: typing.List[typing.Any] = bashlex.parse(line, strictmode=False)
        positions: typing.List[typing.Tuple[int, int]] = []

        for tree in trees:
            visitor = nodevisitor(positions)
            visitor.visit(tree)

        # do replacements from the end so the indicies will be correct
        positions.reverse()

        processed = list(line)
        for start, end in positions:
            processed[start:end] = ""

        processed_line: str = "".join(processed)

        command_dict = self._parse_processed_line(processed_line)
        command_dict["stdin"] = visitor.stdin
        command_dict["stdout"] = visitor.stdout
        command_dict["stderr"] = visitor.stderr

        return command_dict

    def _parse_processed_line(self, line: str) -> typing.Dict[str, typing.Any]:
        argv: typing.List[str] = list(bashlex.split(line))
        self.executable = argv.pop(0)

        parser: argparse.ArgumentParser = argparse.ArgumentParser(add_help=False)
        if "arguments" in self.config:
            arguments_section = self.config["arguments"]

            arg: str
            if "flag_args" in arguments_section:
                for arg in arguments_section["flag_args"].split():
                    parser.add_argument(f"--{arg}", action="store_true")

            if "str_args" in arguments_section:
                for arg in arguments_section["str_args"].split():
                    parser.add_argument(f"--{arg}", type=str)

            if "bool_str_args" in arguments_section:
                for arg in arguments_section["bool_str_args"].split():
                    parser.add_argument(
                        f"--{arg}", nargs="?", const=True, default=False
                    )

            if "int_args" in arguments_section:
                for arg in arguments_section["int_args"].split():
                    parser.add_argument(f"--{arg}", type=int)

            if "bool_int_args" in arguments_section:
                for arg in arguments_section["bool_int_args"].split():
                    parser.add_argument(
                        f"--{arg}", nargs="?", const=True, default=False
                    )

        namespace, rest = parser.parse_known_intermixed_args(argv)
        vars_dict = vars(namespace)

        command_dict: typing.Dict[str, typing.Any] = {
            "executable": self.executable,
            "args": [],
        }

        for var in vars_dict.keys():
            value = vars_dict[var]
            if value is not None:
                if isinstance(value, bool):
                    if value is True:
                        command_dict["args"].append(f"--{var}")
                else:
                    command_dict["args"].append(f"--{var}={value}")

        if len(rest) >= 0:
            value = rest.pop(0)
            command_dict["args"].append(f"--target={value}")

        if len(rest) >= 0:
            value = rest.pop(0)
            command_dict["args"].append(f"--query={value}")

        return command_dict

    def _write_format(self) -> None:
        if self.args.format_selector == "bam":
            format_name = "bam"
        elif self.args.format_selector == "maf":
            format_name = "maf"
        elif self.args.format_selector == "differences":
            format_name = "interval"
        else:
            format_name = "tabular"

        with open("format.txt", "w") as ofh:
            print(f"{format_name}", file=ofh)

        self.package_file.add_format("format.txt")


class nodevisitor(bashlex.ast.nodevisitor):  # type: ignore[misc]
    def __init__(self, positions: typing.List[typing.Tuple[int, int]]) -> None:
        self.positions = positions
        self.stdin = None
        self.stdout = None
        self.stderr = None

    def visitredirect(
        self,
        n: bashlex.ast.node,
        n_input: int,
        n_type: str,
        output: typing.Any,
        heredoc: typing.Any,
    ) -> None:
        if isinstance(n_input, int) and 0 <= n_input <= 2:
            if isinstance(output, bashlex.ast.node) and output.kind == "word":
                self.positions.append(n.pos)
                if n_input == 0:
                    self.stdin = output.word
                elif n_input == 1:
                    self.stdout = output.word
                elif n_input == 2:
                    self.stderr = output.word
            else:
                sys.exit(f"oops 1: {type(n_input)}")
        else:
            sys.exit(f"oops 2: {type(n_input)}")

    def visitheredoc(self, n: bashlex.ast.node, value: typing.Any) -> None:
        pass


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tool_directory", type=str, required=True, help="tool directory")
    parser.add_argument("--format_selector", type=str, required=True, help="format selector")
    parser.add_argument("--debug", action="store_true", help="enable debug messages")
    args = parser.parse_args()

    if args.debug:
        r_beg = resource.getrusage(resource.RUSAGE_SELF)
        beg: int = time.monotonic_ns()

    lastz_command_config_file: str = os.path.join(args.tool_directory, "lastz-cmd.ini")

    config: configparser.ConfigParser = configparser.ConfigParser()
    config.read(lastz_command_config_file)

    package_file = PackageFile()
    lastz_command_file = "lastz-commands.txt"
    bashCommandLineFile(lastz_command_file, config, args, package_file)
    package_file.close()

    if args.debug:
        ns: int = time.monotonic_ns() - beg
        r_end = resource.getrusage(resource.RUSAGE_SELF)
        print(f"package output clock time: {ns} ns", file=sys.stderr, flush=True)
        for rusage_attr in RUSAGE_ATTRS:
            value = getattr(r_end, rusage_attr) - getattr(r_beg, rusage_attr)
            print(f"  package output {rusage_attr}: {value}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
