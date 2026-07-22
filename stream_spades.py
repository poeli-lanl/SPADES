#!/usr/bin/env python3
"""Continuously run SPADES on arriving ONT read files and cumulative alignments."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from datetime import datetime
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
    "observed_at",
    "completed_at",
    "input_file",
    "size_bytes",
    "mtime",
    "raw_reads",
    "filtered_reads",
    "cumulative_raw_reads",
    "cumulative_filtered_reads",
    "qc_json",
    "chunk_bam",
    "cumulative_bam",
    "cumulative_fasta",
    "pathogen_full_tsv",
    "pathogen_full_html",
    "coverage_html",
    "coverage_tsv",
    "run_log",
    "status",
)


def local_now() -> str:
    """Return the current time in the system timezone, including its offset."""
    return datetime.now().astimezone().isoformat(timespec="seconds")


def local_from_timestamp(timestamp: float) -> str:
    """Render a POSIX timestamp in the system timezone."""
    return datetime.fromtimestamp(timestamp).astimezone().isoformat(timespec="seconds")


def timestamp_in_local_timezone(value: Any) -> str:
    """Normalize old UTC state values and new offset timestamps to local time."""
    if value in (None, ""):
        return ""
    text = str(value)
    try:
        parsed = datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError:
        return text
    if parsed.tzinfo is None:
        parsed = parsed.astimezone()
    return parsed.astimezone().isoformat(timespec="seconds")


def is_read_file(path: Path) -> bool:
    return path.is_file() and path.name.lower().endswith(READ_SUFFIXES)


def count_reads(path: Path) -> int:
    """Count FASTA or FASTQ records without loading the input into memory."""
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
        opener = gzip.open
    else:
        opener = open

    if name.endswith((".fa", ".fasta", ".fna")):
        with opener(path, "rb") as handle:
            return sum(1 for line in handle if line.startswith(b">"))

    if name.endswith((".fq", ".fastq")):
        count = 0
        with opener(path, "rb") as handle:
            while True:
                header = handle.readline()
                if not header:
                    return count
                if not header.startswith(b"@"):
                    raise RuntimeError(
                        f"Invalid FASTQ record {count + 1} in {path}: missing '@' header"
                    )

                sequence_length = 0
                while True:
                    line = handle.readline()
                    if not line:
                        raise RuntimeError(
                            f"Truncated FASTQ sequence for record {count + 1} in {path}"
                        )
                    if line.startswith(b"+"):
                        break
                    sequence_length += len(line.rstrip(b"\r\n"))
                if sequence_length == 0:
                    raise RuntimeError(
                        f"Empty FASTQ sequence for record {count + 1} in {path}"
                    )

                quality_length = 0
                while quality_length < sequence_length:
                    line = handle.readline()
                    if not line:
                        raise RuntimeError(
                            f"Truncated FASTQ quality for record {count + 1} in {path}"
                        )
                    quality_length += len(line.rstrip(b"\r\n"))
                if quality_length != sequence_length:
                    raise RuntimeError(
                        f"FASTQ sequence/quality length mismatch for record "
                        f"{count + 1} in {path}"
                    )
                count += 1

    raise RuntimeError(f"Unsupported read file extension: {path}")


def file_signature(path: Path) -> Dict[str, Any]:
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "mtime": local_from_timestamp(stat.st_mtime),
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
        self.active_run_log: Optional[Path] = None
        self.observations: Dict[str, Tuple[int, int, float]] = {}
        self.state = self._load_state()
        self._recover_pending_timepoint()
        if self._backfill_qc_metrics():
            self._save_state()
        else:
            self._write_live_report()
        logging.info(
            "Open the live report in your browser: %s",
            self._relative_path(self.report_path),
        )

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
            "cumulative_fasta": "",
            "latest_outputs": {},
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
        self._normalize_legacy_timestamps(state)
        return state

    @staticmethod
    def _normalize_legacy_timestamps(state: Dict[str, Any]) -> None:
        """Migrate timestamp field names from early UTC-only stream states."""
        records = list(state.get("timepoints", []))
        pending_record = (state.get("pending") or {}).get("record")
        if pending_record:
            records.append(pending_record)
        for record in records:
            for current, legacy in (
                ("observed_at", "observed_at_utc"),
                ("completed_at", "completed_at_utc"),
                ("mtime", "mtime_utc"),
            ):
                value = record.get(current, record.get(legacy, ""))
                if value not in (None, ""):
                    record[current] = timestamp_in_local_timezone(value)
                record.pop(legacy, None)
        pending = state.get("pending") or {}
        signature = pending.get("signature") or {}
        if "mtime_utc" in signature and "mtime" not in signature:
            signature["mtime"] = timestamp_in_local_timezone(signature["mtime_utc"])
        signature.pop("mtime_utc", None)

    @staticmethod
    def _relative_path(path: Any, base: Optional[Path] = None) -> str:
        if path in (None, ""):
            return ""
        base = (base or Path.cwd()).resolve()
        candidate = Path(str(path)).expanduser()
        try:
            resolved = (
                candidate.resolve()
                if candidate.is_absolute()
                else (Path.cwd() / candidate).resolve()
            )
            return os.path.relpath(resolved, base)
        except (OSError, ValueError):
            return str(candidate)

    def _display_command(self, command: Sequence[str]) -> str:
        display = []
        for part in command:
            text = str(part)
            display.append(self._relative_path(text) if Path(text).is_absolute() else text)
        return command_text(display)

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
                manifest_record = dict(record)
                for field in (
                    "input_file",
                    "qc_json",
                    "chunk_bam",
                    "cumulative_bam",
                    "cumulative_fasta",
                    "pathogen_full_tsv",
                    "pathogen_full_html",
                    "coverage_html",
                    "coverage_tsv",
                    "run_log",
                ):
                    manifest_record[field] = self._relative_path(
                        manifest_record.get(field, ""), self.output_dir
                    )
                writer.writerow(
                    {field: manifest_record.get(field, "") for field in MANIFEST_FIELDS}
                )
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
    def _profile_rows(path: Path) -> List[Dict[str, str]]:
        if not path.is_file() or path.stat().st_size == 0:
            return []
        with path.open(encoding="utf-8", newline="") as handle:
            rows = []
            fields = (
                "LEVEL",
                "NAME",
                "TAXID",
                "READ_COUNT",
                "BEST_SIG_COV",
                "SIG_COV",
                "SNI_SCORE",
                "HUMAN_PATHOGEN",
            )
            for row in csv.DictReader(handle, delimiter="\t"):
                snapshot = {field: str(row.get(field, "") or "") for field in fields}
                pathogenic_info = str(row.get("PATHOGENIC_INFO", "") or "").strip()
                snapshot["PATHOGENIC_INFO"] = (
                    "present"
                    if pathogenic_info and pathogenic_info.lower() != "nan"
                    else ""
                )
                rows.append(snapshot)
            return rows

    @staticmethod
    def _input_has_qc(input_file: str) -> bool:
        value = input_file.lower()
        return value.endswith((".fq", ".fastq", ".fq.gz", ".fastq.gz"))

    def _recompute_cumulative_qc(self) -> None:
        cumulative_raw = 0
        cumulative_filtered = 0
        cumulative_filtered_available = True
        for record in sorted(
            self.state.get("timepoints", []), key=lambda item: int(item["timepoint"])
        ):
            if record.get("raw_reads") not in (None, ""):
                cumulative_raw += int(record["raw_reads"])
            if record.get("filtered_reads") in (None, ""):
                cumulative_filtered_available = False
            elif cumulative_filtered_available:
                cumulative_filtered += int(record["filtered_reads"])
            record["cumulative_raw_reads"] = cumulative_raw
            record["cumulative_filtered_reads"] = (
                cumulative_filtered if cumulative_filtered_available else None
            )

    def _backfill_qc_metrics(self) -> bool:
        changed = False
        for record in self.state.get("timepoints", []):
            has_metrics = "raw_reads" in record and "filtered_reads" in record
            skip_qc = getattr(self.args, "skip_qc", False)
            if has_metrics and (
                not skip_qc
                or (
                    record.get("raw_reads") is not None
                    and record.get("filtered_reads") is None
                )
            ):
                continue
            sequence = int(record["timepoint"])
            chunk_dir = self.output_dir / "timepoints" / f"timepoint_{sequence:06d}" / "chunk"
            chunk_prefix = f"{self.args.prefix}.t{sequence:06d}.chunk"
            input_path = Path(str(record.get("input_file", "")))
            if skip_qc and input_path.is_file():
                read_count = count_reads(input_path)
                metrics = {
                    "raw_reads": read_count,
                    "filtered_reads": None,
                    "qc_json": "",
                }
            elif skip_qc:
                if has_metrics:
                    continue
                metrics = {
                    "raw_reads": None,
                    "filtered_reads": None,
                    "qc_json": "",
                }
            else:
                metrics = self._find_qc_metrics(
                    chunk_dir,
                    chunk_prefix,
                    required=self._input_has_qc(str(input_path)),
                )
            record.update(metrics)
            changed = True
        if changed:
            self._recompute_cumulative_qc()
        return changed

    def _run(self, command: Sequence[str]) -> None:
        command_args = [str(part) for part in command]
        display_command = self._display_command(command_args)
        logging.debug("Running: %s", display_command)
        process_log = (
            self.active_run_log if command_args else None
        )
        if process_log is None:
            subprocess.run(command_args, check=True)
            return

        process_log.parent.mkdir(parents=True, exist_ok=True)
        started_at = local_now()
        with process_log.open("a", encoding="utf-8") as log_handle:
            log_handle.write(
                f"\n[{started_at}] COMMAND {display_command}\n"
                "[combined stdout/stderr]\n"
            )
            log_handle.flush()
            process = subprocess.Popen(
                command_args,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                encoding="utf-8",
                errors="replace",
                bufsize=1,
            )
            try:
                assert process.stdout is not None
                with process.stdout:
                    for line in process.stdout:
                        log_handle.write(line)
                        log_handle.flush()
                returncode = process.wait()
            except BaseException:
                process.terminate()
                process.wait()
                raise
            log_handle.write(
                f"\n[{local_now()}] EXIT {returncode}\n"
            )
            log_handle.flush()
        if returncode:
            raise subprocess.CalledProcessError(returncode, command_args)

    def _bam_is_valid(self, path: Optional[Path]) -> bool:
        if path is None or not path.is_file() or path.stat().st_size == 0:
            return False
        return subprocess.run(
            [self.args.samtools, "quickcheck", str(path)],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        ).returncode == 0

    @staticmethod
    def _fasta_is_valid(path: Optional[Path]) -> bool:
        if path is None or not path.is_file() or path.stat().st_size == 0:
            return False
        try:
            with gzip.open(path, "rb") as handle:
                has_content = False
                while chunk := handle.read(1024 * 1024):
                    has_content = True
                return has_content
        except (OSError, EOFError):
            return False

    def _merge_fastas(self, inputs: Sequence[Path], output: Path) -> None:
        merge_script = Path(__file__).resolve().parent / "scripts" / "merge_fasta.py"
        command = [sys.executable, str(merge_script), *[str(path) for path in inputs]]
        plain = output.with_name(f".{output.name}.merge.tmp")
        compressed = output.with_name(f".{output.name}.tmp")
        output.parent.mkdir(parents=True, exist_ok=True)
        try:
            with plain.open("wb") as handle:
                result = subprocess.run(
                    command,
                    stdout=handle,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=False,
                )
            if result.stderr and self.active_run_log is not None:
                self.active_run_log.parent.mkdir(parents=True, exist_ok=True)
                with self.active_run_log.open("a", encoding="utf-8") as log_handle:
                    log_handle.write(
                        f"\n[{local_now()}] COMMAND {self._display_command(command)}\n"
                        f"[stderr]\n{result.stderr}"
                    )
            if result.returncode:
                raise subprocess.CalledProcessError(
                    result.returncode, command, stderr=result.stderr
                )
            with plain.open("rb") as source, gzip.open(compressed, "wb") as destination:
                shutil.copyfileobj(source, destination)
            os.replace(compressed, output)
        finally:
            plain.unlink(missing_ok=True)
            compressed.unlink(missing_ok=True)

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

        candidate_fasta_text = pending.get("candidate_fasta", "")
        cumulative_fasta_text = pending.get("cumulative_fasta", "")
        if candidate_fasta_text:
            candidate_fasta = Path(candidate_fasta_text)
            cumulative_fasta = Path(cumulative_fasta_text)
            cumulative_fasta.parent.mkdir(parents=True, exist_ok=True)
            if candidate_fasta.exists():
                os.replace(candidate_fasta, cumulative_fasta)
            expected_size = int(pending["cumulative_fasta_size"])
            if (
                not cumulative_fasta.is_file()
                or cumulative_fasta.stat().st_size != expected_size
                or not self._fasta_is_valid(cumulative_fasta)
            ):
                raise RuntimeError(
                    "Unable to recover cumulative FASTA for timepoint "
                    f"{pending['timepoint']}"
                )

        has_profile_plan = "profile_outputs" in pending
        self._promote_profile_outputs(pending.get("profile_outputs", []))

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
        self.state["cumulative_fasta"] = cumulative_fasta_text
        latest_outputs = dict(pending.get("latest_outputs", {}))
        if latest_outputs:
            self.state["latest_outputs"] = latest_outputs
        if pending.get("db_level"):
            self.state["db_level"] = pending["db_level"]
        self.state["pending"] = None

        if has_profile_plan and not getattr(self.args, "keep_timepoints", False):
            timepoint_dir = Path(pending["timepoint_dir"])
            try:
                if timepoint_dir.exists():
                    shutil.rmtree(timepoint_dir)
                if timepoint_dir.parent.exists():
                    timepoint_dir.parent.rmdir()
            except OSError as error:
                logging.warning(
                    "Unable to remove completed timepoint files from %s: %s",
                    timepoint_dir,
                    error,
                )

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
            "-o", str(outdir),
            "-p", prefix,
            "-d", str(self.database),
            "-t", str(self.args.cpu),
            "--ont",
            "--ont-error-rate", str(self.args.ont_error_rate),
            "--min-depth", str(self.args.min_depth),
        ]
        if self.args.spades_data:
            command.extend(["--spades-data", str(self.args.spades_data.resolve())])
        if getattr(self.args, "skip_qc", False):
            command.append("--skip-qc")
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
    def _find_chunk_fasta(chunk_dir: Path, chunk_prefix: str) -> Optional[Path]:
        filename = f"{chunk_prefix}.sylph_extracted.fa.gz"
        candidates = [
            path
            for path in (chunk_dir / filename, chunk_dir / "intermediate" / filename)
            if path.is_file()
        ]
        if len(candidates) > 1:
            raise RuntimeError(
                f"Multiple chunk Sylph-extracted FASTA files found for {chunk_prefix}: "
                f"{candidates}"
            )
        return candidates[0] if candidates else None

    @staticmethod
    def _bam_level(path: Path) -> str:
        match = re.search(r"\.gottcha_([^.]+)\.bam$", path.name)
        if not match:
            raise RuntimeError(f"Cannot detect database level from BAM name: {path}")
        return match.group(1)

    @staticmethod
    def _copy_chunk_profile_outputs(
        chunk_dir: Path,
        chunk_prefix: str,
        profile_dir: Path,
        profile_prefix: str,
    ) -> None:
        """Stage report artifacts from an alignment-free chunk run."""
        exact_suffixes = {
            ".tsv",
            ".full.tsv",
            ".krona.html",
            ".pathogen.tsv",
            ".pathogen.summary.txt",
            ".pathogen.full.tsv",
            ".pathogen.summary.tsv",
            ".pathogen.full.html",
            ".coverage.html",
            ".info",
        }
        for source in sorted(chunk_dir.rglob(f"{chunk_prefix}.*")):
            if not source.is_file():
                continue
            suffix = source.name[len(chunk_prefix) :]
            is_coverage = bool(re.fullmatch(r"\.gottcha_[^.]+\.coverage\.tsv", suffix))
            if suffix not in exact_suffixes and not is_coverage:
                continue
            relative_parent = source.parent.relative_to(chunk_dir)
            destination = profile_dir / relative_parent / f"{profile_prefix}{suffix}"
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)

    def _profile_output_plan(
        self, profile_dir: Path, profile_prefix: str
    ) -> Tuple[List[Dict[str, str]], Dict[str, str]]:
        """Map staged profiling files to their stable cumulative locations."""
        cumulative_dir = self.output_dir / "cumulative"
        outputs: List[Dict[str, str]] = []
        latest_outputs: Dict[str, str] = {}
        destinations = set()
        coverage_pattern = re.compile(
            rf"{re.escape(profile_prefix)}\.gottcha_[^.]+\.coverage\.tsv$"
        )

        for source in sorted(profile_dir.rglob("*")):
            if not source.is_file():
                continue
            relative = source.relative_to(profile_dir)
            destination = cumulative_dir / relative
            if coverage_pattern.fullmatch(source.name):
                destination = cumulative_dir / source.name
            destination_key = str(destination.resolve())
            if destination_key in destinations:
                raise RuntimeError(
                    f"Multiple profiling outputs map to the same destination: {destination}"
                )
            destinations.add(destination_key)
            outputs.append(
                {"source": str(source.resolve()), "destination": destination_key}
            )

            if source.name == f"{profile_prefix}.pathogen.full.tsv":
                latest_outputs["pathogen_full_tsv"] = destination_key
            elif source.name == f"{profile_prefix}.pathogen.full.html":
                latest_outputs["pathogen_full_html"] = destination_key
            elif source.name == f"{profile_prefix}.coverage.html":
                latest_outputs["coverage_html"] = destination_key
            elif coverage_pattern.fullmatch(source.name):
                latest_outputs["coverage_tsv"] = destination_key

        if "pathogen_full_tsv" not in latest_outputs:
            raise RuntimeError("The staged profile is missing its pathogen full TSV")
        return outputs, latest_outputs

    def _promote_profile_outputs(self, outputs: Sequence[Dict[str, str]]) -> None:
        """Publish a staged profile, retaining its staged copy only on request."""
        keep_timepoints = getattr(self.args, "keep_timepoints", False)
        for item in outputs:
            source = Path(item["source"])
            destination = Path(item["destination"])
            destination.parent.mkdir(parents=True, exist_ok=True)
            if source.is_file():
                if keep_timepoints:
                    temporary = destination.with_name(f".{destination.name}.tmp")
                    shutil.copy2(source, temporary)
                    os.replace(temporary, destination)
                else:
                    os.replace(source, destination)
            elif not destination.is_file():
                raise RuntimeError(f"Unable to recover profiling output: {destination}")

    def process_file(self, path: Path, signature: Dict[str, Any]) -> None:
        sequence = int(self.state["next_timepoint"])
        label = f"timepoint_{sequence:06d}"
        timepoint_dir = self.output_dir / "timepoints" / label
        chunk_dir = timepoint_dir / "chunk"
        profile_dir = timepoint_dir / "profile"
        chunk_prefix = f"{self.args.prefix}.t{sequence:06d}.chunk"
        profile_prefix = self.args.prefix
        log_prefix = f"{self.args.prefix}.t{sequence:06d}"
        observed_at = local_now()
        # Keep process output outside the disposable in-progress timepoint tree so
        # failed/interrupted attempts remain available after recovery.
        run_log = self.output_dir / "logs" / f"{log_prefix}.run_SPADES.log"
        self.active_run_log = run_log

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

        logging.info("Processing %s as %s", self._relative_path(path), label)
        chunk_command = self._spades_base_command(chunk_dir, chunk_prefix)
        chunk_command[1:1] = ["-i", str(path.resolve())]
        self._run(chunk_command)

        if getattr(self.args, "skip_qc", False):
            read_count = count_reads(path)
            qc_metrics = {
                "raw_reads": read_count,
                "filtered_reads": None,
                "qc_json": "",
            }
        else:
            qc_metrics = self._find_qc_metrics(
                chunk_dir,
                chunk_prefix,
                required=self._input_has_qc(str(path)),
            )

        chunk_bam = self._find_chunk_bam(chunk_dir, chunk_prefix)
        valid_chunk_bam = self._bam_is_valid(chunk_bam)
        chunk_fasta = self._find_chunk_fasta(chunk_dir, chunk_prefix)
        if valid_chunk_bam and not self._fasta_is_valid(chunk_fasta):
            raise RuntimeError(
                f"Valid chunk BAM has no valid Sylph-extracted FASTA: {chunk_fasta}"
            )
        previous_text = self.state.get("cumulative_bam", "")
        previous_bam = Path(previous_text) if previous_text else None
        if previous_bam is not None and not self._bam_is_valid(previous_bam):
            raise RuntimeError(f"Saved cumulative BAM is missing or invalid: {previous_bam}")
        previous_fasta_text = self.state.get("cumulative_fasta", "")
        if previous_bam is not None and not previous_fasta_text:
            inferred_fasta = (
                self.output_dir / "cumulative" / f"{self.args.prefix}.sylph_extracted.fa.gz"
            )
            if inferred_fasta.is_file():
                previous_fasta_text = str(inferred_fasta)
        previous_fasta = Path(previous_fasta_text) if previous_fasta_text else None
        if previous_bam is not None and not self._fasta_is_valid(previous_fasta):
            raise RuntimeError(
                "Saved cumulative BAM has no valid cumulative Sylph-extracted FASTA; "
                "rebuild the streaming output to call variants cumulatively"
            )

        candidate: Optional[Path] = None
        cumulative: Optional[Path] = None
        candidate_fasta: Optional[Path] = None
        cumulative_fasta: Optional[Path] = None
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
                "-@", str(self.args.cpu),
                str(merged_unsorted),
                *[str(item) for item in merge_inputs],
            ]
            self._run(merge_command)
            self._run(
                [
                    self.args.samtools,
                    "sort",
                    "-@", str(self.args.cpu),
                    "-o", str(candidate),
                    str(merged_unsorted),
                ]
            )
            merged_unsorted.unlink()
            self._run(
                [self.args.samtools, "index", "-@", str(self.args.cpu), str(candidate)]
            )
            if not self._bam_is_valid(candidate):
                raise RuntimeError(f"Merged BAM failed validation: {candidate}")
            else:
                logging.info(
                    "Merged cumulative BAM for timepoint %s: %s",
                    sequence,
                    self._relative_path(candidate),
                )

            candidate_fasta = timepoint_dir / "candidate.sylph_extracted.fa.gz"
            cumulative_fasta = (
                self.output_dir
                / "cumulative"
                / f"{self.args.prefix}.sylph_extracted.fa.gz"
            )
            fasta_inputs = []
            if previous_fasta is not None:
                fasta_inputs.append(previous_fasta)
            if valid_chunk_bam and chunk_fasta is not None:
                fasta_inputs.append(chunk_fasta)
            self._merge_fastas(fasta_inputs, candidate_fasta)
            if not self._fasta_is_valid(candidate_fasta):
                raise RuntimeError(
                    f"Merged Sylph-extracted FASTA failed validation: {candidate_fasta}"
                )
            else:
                logging.info(
                    "Merged cumulative Sylph-extracted FASTA for timepoint %s: %s",
                    sequence,
                    self._relative_path(candidate_fasta),
                )

            profile_dir.mkdir(parents=True)
            profile_reference = (
                profile_dir / f"{profile_prefix}.sylph_extracted.fa.gz"
            )
            # BAM mode discovers its variant-calling reference by this output name.
            # A hard link exposes the merged FASTA without copying it.
            os.link(candidate_fasta, profile_reference)
            profile_command = self._spades_base_command(profile_dir, profile_prefix)
            profile_command[1:1] = ["-b", str(candidate)]
            try:
                self._run(profile_command)
            finally:
                # run_SPADES.sh moves prefix-matched intermediates after profiling.
                # Remove only the temporary reference link; candidate_fasta is
                # promoted atomically once the full cumulative profile is ready.
                profile_reference.unlink(missing_ok=True)
                (profile_dir / "intermediate" / profile_reference.name).unlink(
                    missing_ok=True
                )
        else:
            status = "no_alignments"
            profile_dir.mkdir(parents=True)
            self._copy_chunk_profile_outputs(
                chunk_dir, chunk_prefix, profile_dir, profile_prefix
            )
            destination = profile_dir / f"{profile_prefix}.pathogen.full.tsv"
            if not destination.exists():
                destination.touch()

        result = profile_dir / f"{profile_prefix}.pathogen.full.tsv"
        if not result.is_file():
            raise RuntimeError(f"Final pathogen result was not created: {result}")
        profile_outputs, latest_outputs = self._profile_output_plan(
            profile_dir, profile_prefix
        )
        final_result = Path(latest_outputs["pathogen_full_tsv"])
        keep_timepoints = getattr(self.args, "keep_timepoints", False)

        record = {
            "timepoint": sequence,
            "observed_at": observed_at,
            "completed_at": local_now(),
            "input_file": signature["path"],
            "size_bytes": signature["size"],
            "mtime": signature["mtime"],
            **qc_metrics,
            "chunk_bam": (
                str(chunk_bam.resolve())
                if keep_timepoints and valid_chunk_bam and chunk_bam
                else ""
            ),
            "cumulative_bam": str(cumulative) if cumulative else previous_text,
            "cumulative_fasta": (
                str(cumulative_fasta) if cumulative_fasta else previous_fasta_text
            ),
            "pathogen_full_tsv": str(final_result),
            "pathogen_full_html": latest_outputs.get("pathogen_full_html", ""),
            "coverage_html": latest_outputs.get("coverage_html", ""),
            "coverage_tsv": latest_outputs.get("coverage_tsv", ""),
            "profile_rows": self._profile_rows(result),
            "run_log": str(run_log.resolve()),
            "status": status,
        }
        if not keep_timepoints:
            record["qc_json"] = ""
        prior_raw_reads = sum(
            int(item["raw_reads"])
            for item in self.state.get("timepoints", [])
            if item.get("raw_reads") not in (None, "")
        )
        previous_timepoints = self.state.get("timepoints", [])
        prior_filtered_reads = sum(
            int(item["filtered_reads"])
            for item in previous_timepoints
            if item.get("filtered_reads") not in (None, "")
        )
        record["cumulative_raw_reads"] = prior_raw_reads + int(qc_metrics["raw_reads"] or 0)
        filtered_reads_available = (
            qc_metrics["filtered_reads"] not in (None, "")
            and all(
                item.get("filtered_reads") not in (None, "")
                for item in previous_timepoints
            )
        )
        record["cumulative_filtered_reads"] = (
            prior_filtered_reads + int(qc_metrics["filtered_reads"])
            if filtered_reads_available
            else None
        )
        self.state["pending"] = {
            "stage": "ready",
            "timepoint": sequence,
            "timepoint_dir": str(timepoint_dir),
            "signature": signature,
            "candidate_bam": str(candidate) if candidate else "",
            "cumulative_bam": str(cumulative) if cumulative else previous_text,
            "cumulative_size": candidate.stat().st_size if candidate else 0,
            "candidate_fasta": str(candidate_fasta) if candidate_fasta else "",
            "cumulative_fasta": (
                str(cumulative_fasta) if cumulative_fasta else previous_fasta_text
            ),
            "cumulative_fasta_size": (
                candidate_fasta.stat().st_size if candidate_fasta else 0
            ),
            "db_level": db_level,
            "profile_outputs": profile_outputs,
            "latest_outputs": latest_outputs,
            "record": record,
        }
        self._save_state()
        self._finalize_pending()
        self.observations.pop(signature["path"], None)
        logging.info(
            "Saved latest cumulative profile: %s", self._relative_path(final_result)
        )

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
        ),
        epilog="""
How it works:
  New or changed FASTA/FASTQ files become timepoints after their size and
  modification time remain unchanged for --settle-seconds. Results and state
  are retained in OUTDIR, so restarting with the same options resumes safely.

Example (continuous monitoring):
  %(prog)s --input-dir /data/ont/fastq_pass \\
    --outdir /results/run_01_stream --prefix run_01 \\
    --db-path /db/gottcha_db.species.fna --cpu 8

Example (process the current files and exit):
  %(prog)s -i /data/ont/fastq_pass -o /results/run_01_stream \\
    -p run_01 -d /db/gottcha_db.species.fna --once --settle-seconds 0

Open OUTDIR/PREFIX.stream.html to follow the live report. Pipeline output is
saved under OUTDIR/logs/. Press Ctrl-C to stop continuous monitoring; run the
same command again to resume.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input-dir", type=Path, required=True, metavar="DIR",
        help="Directory containing arriving ONT FASTA/FASTQ files",
    )
    parser.add_argument(
        "-o", "--outdir", type=Path, required=True, metavar="DIR",
        help="Output directory for cumulative results, state, logs, and the live report",
    )
    parser.add_argument(
        "-p", "--prefix", required=True, metavar="NAME",
        help="Filename prefix for reports and analysis results",
    )
    parser.add_argument(
        "-d", "--db-path", type=Path, required=True, metavar="PATH",
        help="GOTTCHA2 database base path, such as gottcha_db.species.fna",
    )
    parser.add_argument(
        "-t", "--cpu", type=positive_int, default=1, metavar="N",
        help="Worker threads passed to SPADES and samtools (default: %(default)s)",
    )
    parser.add_argument(
        "--run-spades", type=Path, default=script_dir / "run_SPADES.sh", metavar="PATH",
        help="Path to run_SPADES.sh (default: next to this script)",
    )
    parser.add_argument(
        "--spades-data", type=Path, metavar="DIR",
        help="Data directory passed to run_SPADES.sh",
    )
    parser.add_argument(
        "--samtools", default="samtools", metavar="COMMAND",
        help="samtools executable or path (default: %(default)s)",
    )
    parser.add_argument(
        "--poll-interval", type=nonnegative_float, default=5.0, metavar="SECONDS",
        help="Delay between input-directory scans (default: %(default)s)",
    )
    parser.add_argument(
        "--settle-seconds", type=nonnegative_float, default=30.0, metavar="SECONDS",
        help="Time a file must remain unchanged before processing (default: %(default)s)",
    )
    parser.add_argument(
        "--ont-error-rate", type=nonnegative_float, default=0.03, metavar="RATE",
        help="ONT error rate passed to run_SPADES.sh (default: %(default)s)",
    )
    parser.add_argument(
        "--min-depth", type=positive_int, default=10, metavar="N",
        help="Minimum variant-calling depth passed to run_SPADES.sh (default: %(default)s)",
    )
    parser.add_argument(
        "--skip-qc", action="store_true",
        help=(
            "Skip fastp/fastplong quality control and preprocessing; count reads "
            "directly from each input"
        ),
    )
    recursive_group = parser.add_mutually_exclusive_group()
    recursive_group.add_argument(
        "--recursive", dest="recursive", action="store_true", help=argparse.SUPPRESS,
    )
    recursive_group.add_argument(
        "--no-recursive", dest="recursive", action="store_false",
        help="Watch only the input directory, not its subdirectories",
    )
    parser.set_defaults(recursive=True)
    parser.add_argument(
        "--keep-timepoints", action="store_true",
        help="Keep per-timepoint chunk and staged profiling files (default: discard them)",
    )
    parser.add_argument(
        "--once", action="store_true", help="Process currently stable files and exit"
    )
    parser.add_argument(
        "--max-files", type=positive_int, metavar="N",
        help="Exit after processing this many files"
    )
    parser.add_argument(
        "--verbose", action="store_true", help="Show detailed diagnostic messages"
    )
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
        format="[%(asctime)s] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    logging.info("Starting SPADES-GOTTCHA2 in stream mode...")
    logging.info("Monitoring input directory: %s", args.input_dir)
    logging.info("Writing results, state, and logs to: %s", args.outdir)
    if args.once:
        logging.info("One-shot mode: process currently stable files, then exit")
    else:
        logging.info("Continuous mode: press Ctrl-C to stop; rerun to resume")

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
