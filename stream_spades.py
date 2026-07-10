#!/usr/bin/env python3
"""Continuously run SPADES on arriving ONT read files and cumulative alignments."""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from scripts.stream_report import generate_stream_report


READ_SUFFIXES = (
    ".fa",
    ".fasta",
    ".fna",
    ".fq",
    ".fastq",
    ".fa.gz",
    ".fasta.gz",
    ".fna.gz",
    ".fq.gz",
    ".fastq.gz",
)
STATE_VERSION = 1
MANIFEST_FIELDS = (
    "timepoint",
    "observed_at_utc",
    "completed_at_utc",
    "input_file",
    "size_bytes",
    "mtime_utc",
    "raw_reads",
    "filtered_reads",
    "cumulative_raw_reads",
    "cumulative_filtered_reads",
    "qc_json",
    "chunk_bam",
    "cumulative_bam",
    "pathogen_full_tsv",
    "status",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def utc_from_timestamp(timestamp: float) -> str:
    return datetime.fromtimestamp(timestamp, timezone.utc).isoformat(timespec="seconds")


def is_read_file(path: Path) -> bool:
    return path.is_file() and path.name.lower().endswith(READ_SUFFIXES)


def file_signature(path: Path) -> Dict[str, Any]:
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "mtime_utc": utc_from_timestamp(stat.st_mtime),
    }


def database_base(path: Path) -> Path:
    value = str(path)
    for suffix in (".tax.tsv", ".syldb", ".stats", ".zip"):
        if value.endswith(suffix):
            value = value[: -len(suffix)]
            break
    return Path(value)


def atomic_json_write(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def command_text(command: Sequence[str]) -> str:
    return shlex.join(str(part) for part in command)


class StreamSpades:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.input_dir = args.input_dir.resolve()
        self.output_dir = args.outdir.resolve()
        self.run_spades = args.run_spades.resolve()
        self.database = database_base(args.db_path).resolve()
        self.state_path = self.output_dir / "stream_state.json"
        self.manifest_path = self.output_dir / "timepoints.tsv"
        self.report_path = self.output_dir / f"{args.prefix}.stream.html"
        self.observations: Dict[str, Tuple[int, int, float]] = {}
        self.state = self._load_state()
        self._recover_pending_timepoint()
        if self._backfill_qc_metrics():
            self._save_state()
        else:
            self._write_live_report()

    def _initial_state(self) -> Dict[str, Any]:
        return {
            "version": STATE_VERSION,
            "configuration": {
                "input_dir": str(self.input_dir),
                "database": str(self.database),
                "prefix": self.args.prefix,
            },
            "next_timepoint": 1,
            "files": {},
            "timepoints": [],
            "cumulative_bam": "",
            "db_level": "",
            "monitor_status": "stopped",
            "pending": None,
        }

    def _load_state(self) -> Dict[str, Any]:
        if not self.state_path.exists():
            return self._initial_state()
        with self.state_path.open(encoding="utf-8") as handle:
            state = json.load(handle)
        if state.get("version") != STATE_VERSION:
            raise RuntimeError(
                f"Unsupported state version in {self.state_path}: {state.get('version')}"
            )
        expected = self._initial_state()["configuration"]
        if state.get("configuration") != expected:
            raise RuntimeError(
                "The existing stream state belongs to a different input directory, "
                "database, or prefix. Choose another --outdir or restore the original options."
            )
        return state

    def _save_state(self) -> None:
        atomic_json_write(self.state_path, self.state)
        self._write_manifest()
        self._write_live_report()

    def _write_live_report(self) -> None:
        generate_stream_report(self.state, self.report_path)

    def _write_manifest(self) -> None:
        self.output_dir.mkdir(parents=True, exist_ok=True)
        temporary = self.manifest_path.with_name(f".{self.manifest_path.name}.tmp")
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t")
            writer.writeheader()
            for record in sorted(
                self.state["timepoints"], key=lambda item: int(item["timepoint"])
            ):
                writer.writerow({field: record.get(field, "") for field in MANIFEST_FIELDS})
        os.replace(temporary, self.manifest_path)

    @staticmethod
    def _qc_metrics_from_json(path: Path) -> Dict[str, Any]:
        try:
            with path.open(encoding="utf-8") as handle:
                payload = json.load(handle)
            summary = payload["summary"]
            raw_reads = int(summary["before_filtering"]["total_reads"])
            filtered_reads = int(summary["after_filtering"]["total_reads"])
        except (KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
            raise RuntimeError(f"Unable to read QC totals from {path}: {error}") from error
        return {
            "raw_reads": raw_reads,
            "filtered_reads": filtered_reads,
            "qc_json": str(path.resolve()),
        }

    def _find_qc_metrics(
        self, chunk_dir: Path, chunk_prefix: str, required: bool = False
    ) -> Dict[str, Any]:
        candidates = []
        for directory in (chunk_dir / "intermediate", chunk_dir):
            expected = directory / f"{chunk_prefix}.qc.fastq.gz.json"
            if expected.is_file():
                candidates.append(expected)
        candidates = sorted(set(candidates))
        if len(candidates) > 1:
            raise RuntimeError(f"Multiple QC JSON files found for {chunk_prefix}: {candidates}")
        if not candidates:
            if required:
                raise RuntimeError(
                    f"FastP/FastPlong QC JSON was not created for {chunk_prefix}"
                )
            return {"raw_reads": None, "filtered_reads": None, "qc_json": ""}
        return self._qc_metrics_from_json(candidates[0])

    @staticmethod
    def _input_has_qc(input_file: str) -> bool:
        value = input_file.lower()
        return value.endswith((".fq", ".fastq", ".fq.gz", ".fastq.gz"))

    def _recompute_cumulative_qc(self) -> None:
        cumulative_raw = 0
        cumulative_filtered = 0
        for record in sorted(
            self.state.get("timepoints", []), key=lambda item: int(item["timepoint"])
        ):
            if record.get("raw_reads") not in (None, ""):
                cumulative_raw += int(record["raw_reads"])
            if record.get("filtered_reads") not in (None, ""):
                cumulative_filtered += int(record["filtered_reads"])
            record["cumulative_raw_reads"] = cumulative_raw
            record["cumulative_filtered_reads"] = cumulative_filtered

    def _backfill_qc_metrics(self) -> bool:
        changed = False
        for record in self.state.get("timepoints", []):
            if "raw_reads" in record and "filtered_reads" in record:
                continue
            sequence = int(record["timepoint"])
            chunk_dir = self.output_dir / "timepoints" / f"timepoint_{sequence:06d}" / "chunk"
            chunk_prefix = f"{self.args.prefix}.t{sequence:06d}.chunk"
            metrics = self._find_qc_metrics(
                chunk_dir,
                chunk_prefix,
                required=self._input_has_qc(str(record.get("input_file", ""))),
            )
            record.update(metrics)
            changed = True
        if changed:
            self._recompute_cumulative_qc()
        return changed

    def _run(self, command: Sequence[str]) -> None:
        logging.info("Running: %s", command_text(command))
        subprocess.run([str(part) for part in command], check=True)

    def _bam_is_valid(self, path: Optional[Path]) -> bool:
        if path is None or not path.is_file() or path.stat().st_size == 0:
            return False
        return subprocess.run(
            [self.args.samtools, "quickcheck", str(path)],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        ).returncode == 0

    def _finalize_pending(self) -> None:
        pending = self.state["pending"]
        if not pending or pending.get("stage") != "ready":
            raise RuntimeError("Cannot finalize an incomplete pending timepoint")

        candidate_text = pending.get("candidate_bam", "")
        cumulative_text = pending.get("cumulative_bam", "")
        if candidate_text:
            candidate = Path(candidate_text)
            candidate_index = Path(f"{candidate}.bai")
            cumulative = Path(cumulative_text)
            cumulative_index = Path(f"{cumulative}.bai")
            cumulative.parent.mkdir(parents=True, exist_ok=True)

            if candidate.exists():
                os.replace(candidate, cumulative)
            if candidate_index.exists():
                os.replace(candidate_index, cumulative_index)

            expected_size = int(pending["cumulative_size"])
            if not cumulative.is_file() or cumulative.stat().st_size != expected_size:
                raise RuntimeError(
                    f"Unable to recover cumulative BAM for timepoint {pending['timepoint']}"
                )
            if not cumulative_index.is_file() or not self._bam_is_valid(cumulative):
                raise RuntimeError(f"Invalid cumulative BAM after promotion: {cumulative}")

        signature = pending["signature"]
        record = pending["record"]
        self.state["files"][signature["path"]] = {
            "size": signature["size"],
            "mtime_ns": signature["mtime_ns"],
            "timepoint": pending["timepoint"],
        }
        self.state["timepoints"].append(record)
        self.state["next_timepoint"] = int(pending["timepoint"]) + 1
        self.state["cumulative_bam"] = cumulative_text
        if pending.get("db_level"):
            self.state["db_level"] = pending["db_level"]
        self.state["pending"] = None
        self._save_state()

    def _recover_pending_timepoint(self) -> None:
        pending = self.state.get("pending")
        if not pending:
            self._write_manifest()
            self._write_live_report()
            return
        if pending.get("stage") == "ready":
            logging.warning("Recovering completed timepoint %s", pending["timepoint"])
            self._finalize_pending()
            return

        timepoint_dir = Path(pending["timepoint_dir"])
        logging.warning("Removing incomplete timepoint directory: %s", timepoint_dir)
        if timepoint_dir.exists():
            shutil.rmtree(timepoint_dir)
        self.state["pending"] = None
        self._save_state()

    def _already_processed(self, signature: Dict[str, Any]) -> bool:
        previous = self.state["files"].get(signature["path"])
        return bool(
            previous
            and int(previous["size"]) == signature["size"]
            and int(previous["mtime_ns"]) == signature["mtime_ns"]
        )

    def scan(self) -> List[Path]:
        iterator: Iterable[Path]
        iterator = self.input_dir.rglob("*") if self.args.recursive else self.input_dir.iterdir()
        files = []
        for path in iterator:
            try:
                resolved = path.resolve()
                if self.output_dir == resolved or self.output_dir in resolved.parents:
                    continue
                if is_read_file(path):
                    files.append(path)
            except FileNotFoundError:
                continue
        return sorted(files, key=lambda item: (item.stat().st_mtime_ns, str(item)))

    def ready_files(self) -> List[Tuple[Path, Dict[str, Any]]]:
        now_monotonic = time.monotonic()
        now_wall = time.time()
        ready = []
        present = set()
        for path in self.scan():
            try:
                signature = file_signature(path)
            except FileNotFoundError:
                continue
            key = signature["path"]
            present.add(key)
            if self._already_processed(signature):
                self.observations.pop(key, None)
                continue

            previous = self.observations.get(key)
            values = (signature["size"], signature["mtime_ns"])
            if previous is None or previous[:2] != values:
                age = max(0.0, now_wall - path.stat().st_mtime)
                stable_since = now_monotonic - min(age, self.args.settle_seconds)
                self.observations[key] = (*values, stable_since)
                previous = self.observations[key]

            stable_for = now_monotonic - previous[2]
            mtime_age = now_wall - path.stat().st_mtime
            if stable_for >= self.args.settle_seconds and mtime_age >= self.args.settle_seconds:
                ready.append((path, signature))

        for missing in set(self.observations) - present:
            self.observations.pop(missing, None)
        return ready

    def _spades_base_command(self, outdir: Path, prefix: str) -> List[str]:
        command = [
            str(self.run_spades),
            "-o",
            str(outdir),
            "-p",
            prefix,
            "-d",
            str(self.database),
            "-t",
            str(self.args.cpu),
            "--ont",
            "--ont-error-rate",
            str(self.args.ont_error_rate),
            "--min-depth",
            str(self.args.min_depth),
        ]
        if self.args.spades_data:
            command.extend(["--spades-data", str(self.args.spades_data.resolve())])
        if self.args.js_external:
            command.append("--js-external")
        return command

    @staticmethod
    def _find_chunk_bam(chunk_dir: Path, chunk_prefix: str) -> Optional[Path]:
        candidates = list((chunk_dir / "intermediate").glob(f"{chunk_prefix}.gottcha_*.bam"))
        candidates.extend(chunk_dir.glob(f"{chunk_prefix}.gottcha_*.bam"))
        candidates = sorted(set(candidates))
        if len(candidates) > 1:
            raise RuntimeError(
                f"Multiple chunk BAM files found for {chunk_prefix}: {candidates}"
            )
        return candidates[0] if candidates else None

    @staticmethod
    def _bam_level(path: Path) -> str:
        match = re.search(r"\.gottcha_([^.]+)\.bam$", path.name)
        if not match:
            raise RuntimeError(f"Cannot detect database level from BAM name: {path}")
        return match.group(1)

    def process_file(self, path: Path, signature: Dict[str, Any]) -> None:
        sequence = int(self.state["next_timepoint"])
        label = f"timepoint_{sequence:06d}"
        timepoint_dir = self.output_dir / "timepoints" / label
        chunk_dir = timepoint_dir / "chunk"
        profile_dir = timepoint_dir / "profile"
        chunk_prefix = f"{self.args.prefix}.t{sequence:06d}.chunk"
        profile_prefix = f"{self.args.prefix}.t{sequence:06d}"
        observed_at = utc_now()

        if timepoint_dir.exists():
            shutil.rmtree(timepoint_dir)
        chunk_dir.mkdir(parents=True)

        self.state["pending"] = {
            "stage": "started",
            "timepoint": sequence,
            "timepoint_dir": str(timepoint_dir),
            "signature": signature,
        }
        self._save_state()

        logging.info("Processing %s as %s", path, label)
        chunk_command = self._spades_base_command(chunk_dir, chunk_prefix)
        chunk_command[1:1] = ["-i", str(path.resolve())]
        self._run(chunk_command)

        qc_metrics = self._find_qc_metrics(
            chunk_dir,
            chunk_prefix,
            required=self._input_has_qc(str(path)),
        )

        chunk_bam = self._find_chunk_bam(chunk_dir, chunk_prefix)
        valid_chunk_bam = self._bam_is_valid(chunk_bam)
        previous_text = self.state.get("cumulative_bam", "")
        previous_bam = Path(previous_text) if previous_text else None
        if previous_bam is not None and not self._bam_is_valid(previous_bam):
            raise RuntimeError(f"Saved cumulative BAM is missing or invalid: {previous_bam}")

        candidate: Optional[Path] = None
        cumulative: Optional[Path] = None
        db_level = self.state.get("db_level", "")
        status = "profiled"

        if valid_chunk_bam or previous_bam is not None:
            if valid_chunk_bam and chunk_bam is not None:
                chunk_level = self._bam_level(chunk_bam)
                if db_level and db_level != chunk_level:
                    raise RuntimeError(
                        f"BAM database level changed from {db_level} to {chunk_level}"
                    )
                db_level = chunk_level
            if not db_level:
                raise RuntimeError("Cannot determine the database level for cumulative BAM")

            candidate = timepoint_dir / f"candidate.gottcha_{db_level}.bam"
            merged_unsorted = timepoint_dir / "merged.unsorted.bam"
            cumulative = self.output_dir / "cumulative" / (
                f"{self.args.prefix}.gottcha_{db_level}.bam"
            )
            merge_inputs = []
            if previous_bam is not None:
                merge_inputs.append(previous_bam)
            if valid_chunk_bam and chunk_bam is not None:
                merge_inputs.append(chunk_bam)
            merge_command = [
                self.args.samtools,
                "merge",
                "-f",
                "-@",
                str(self.args.cpu),
                str(merged_unsorted),
                *[str(item) for item in merge_inputs],
            ]
            self._run(merge_command)
            self._run(
                [
                    self.args.samtools,
                    "sort",
                    "-@",
                    str(self.args.cpu),
                    "-o",
                    str(candidate),
                    str(merged_unsorted),
                ]
            )
            merged_unsorted.unlink()
            self._run(
                [self.args.samtools, "index", "-@", str(self.args.cpu), str(candidate)]
            )
            if not self._bam_is_valid(candidate):
                raise RuntimeError(f"Merged BAM failed validation: {candidate}")

            profile_dir.mkdir(parents=True)
            profile_command = self._spades_base_command(profile_dir, profile_prefix)
            profile_command[1:1] = ["-b", str(candidate)]
            self._run(profile_command)
        else:
            status = "no_alignments"
            profile_dir.mkdir(parents=True)
            source = chunk_dir / f"{chunk_prefix}.pathogen.full.tsv"
            destination = profile_dir / f"{profile_prefix}.pathogen.full.tsv"
            if source.exists():
                shutil.copy2(source, destination)
            else:
                destination.touch()

        result = profile_dir / f"{profile_prefix}.pathogen.full.tsv"
        if not result.is_file():
            raise RuntimeError(f"Final pathogen result was not created: {result}")

        record = {
            "timepoint": sequence,
            "observed_at_utc": observed_at,
            "completed_at_utc": utc_now(),
            "input_file": signature["path"],
            "size_bytes": signature["size"],
            "mtime_utc": signature["mtime_utc"],
            **qc_metrics,
            "chunk_bam": str(chunk_bam.resolve()) if valid_chunk_bam and chunk_bam else "",
            "cumulative_bam": str(cumulative) if cumulative else previous_text,
            "pathogen_full_tsv": str(result.resolve()),
            "status": status,
        }
        prior_raw_reads = sum(
            int(item["raw_reads"])
            for item in self.state.get("timepoints", [])
            if item.get("raw_reads") not in (None, "")
        )
        prior_filtered_reads = sum(
            int(item["filtered_reads"])
            for item in self.state.get("timepoints", [])
            if item.get("filtered_reads") not in (None, "")
        )
        record["cumulative_raw_reads"] = prior_raw_reads + int(qc_metrics["raw_reads"] or 0)
        record["cumulative_filtered_reads"] = prior_filtered_reads + int(
            qc_metrics["filtered_reads"] or 0
        )
        self.state["pending"] = {
            "stage": "ready",
            "timepoint": sequence,
            "timepoint_dir": str(timepoint_dir),
            "signature": signature,
            "candidate_bam": str(candidate) if candidate else "",
            "cumulative_bam": str(cumulative) if cumulative else previous_text,
            "cumulative_size": candidate.stat().st_size if candidate else 0,
            "db_level": db_level,
            "record": record,
        }
        self._save_state()
        self._finalize_pending()
        self.observations.pop(signature["path"], None)
        logging.info("Saved final timepoint result: %s", result)

    def run(self) -> None:
        processed_this_run = 0
        self.state["monitor_status"] = "running"
        self._save_state()
        try:
            while True:
                ready = self.ready_files()
                if ready:
                    # Process exactly one timepoint at a time. Rescan only after its
                    # chunk run, cumulative merge/profile, and state save complete.
                    path, signature = ready[0]
                    if len(ready) > 1:
                        logging.info(
                            "%d stable files are queued; processing one at a time",
                            len(ready),
                        )
                    self.process_file(path, signature)
                    processed_this_run += 1
                    if self.args.max_files and processed_this_run >= self.args.max_files:
                        return
                    continue
                if self.args.once:
                    return
                time.sleep(self.args.poll_interval)
        finally:
            self.state["monitor_status"] = "stopped"
            self._save_state()


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def nonnegative_float(value: str) -> float:
    parsed = float(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be non-negative")
    return parsed


def build_parser() -> argparse.ArgumentParser:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description=(
            "Monitor an Oxford Nanopore read directory, run SPADES for every stable "
            "FASTA/FASTQ file, merge each BAM into cumulative alignments, and re-run "
            "GOTTCHA2 profiling at every timepoint."
        )
    )
    parser.add_argument("-i", "--input-dir", type=Path, required=True)
    parser.add_argument("-o", "--outdir", type=Path, required=True)
    parser.add_argument("-p", "--prefix", required=True)
    parser.add_argument("-d", "--db-path", type=Path, required=True)
    parser.add_argument("-t", "--cpu", type=positive_int, default=1)
    parser.add_argument(
        "--run-spades", type=Path, default=script_dir / "run_SPADES.sh"
    )
    parser.add_argument("--spades-data", type=Path)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--poll-interval", type=nonnegative_float, default=5.0)
    parser.add_argument("--settle-seconds", type=nonnegative_float, default=30.0)
    parser.add_argument("--ont-error-rate", type=nonnegative_float, default=0.03)
    parser.add_argument("--min-depth", type=positive_int, default=10)
    parser.add_argument("--recursive", action="store_true")
    parser.add_argument("--js-external", action="store_true")
    parser.add_argument(
        "--once", action="store_true", help="Process currently stable files and exit"
    )
    parser.add_argument(
        "--max-files", type=positive_int, help="Exit after this many files"
    )
    parser.add_argument("--verbose", action="store_true")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not args.input_dir.is_dir():
        raise ValueError(f"Input directory not found: {args.input_dir}")
    if not args.run_spades.is_file():
        raise ValueError(f"run_SPADES.sh not found: {args.run_spades}")
    database = str(database_base(args.db_path))
    for suffix in (".syldb", ".zip", ".stats", ".tax.tsv"):
        if not Path(f"{database}{suffix}").is_file():
            raise ValueError(f"Database file not found: {database}{suffix}")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )
    try:
        validate_args(args)
        StreamSpades(args).run()
    except KeyboardInterrupt:
        logging.info("Monitoring stopped by user")
        return 130
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        logging.error("%s", error)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
