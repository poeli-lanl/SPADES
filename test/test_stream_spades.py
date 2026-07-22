import gzip
import json
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path

import stream_spades
from scripts.stream_report import build_payload


class FakeStreamSpades(stream_spades.StreamSpades):
    def _bam_is_valid(self, path):
        return bool(path and path.is_file() and path.stat().st_size > 0)

    def _run(self, command):
        command = [str(part) for part in command]
        if command[0] == str(self.run_spades):
            outdir = Path(command[command.index("-o") + 1])
            prefix = command[command.index("-p") + 1]
            outdir.mkdir(parents=True, exist_ok=True)
            if "-i" in command:
                self.processing_report = (
                    self.output_dir / "stream.stream.html"
                ).read_text(encoding="utf-8")
                source = Path(command[command.index("-i") + 1])
                intermediate = outdir / "intermediate"
                intermediate.mkdir()
                bam = intermediate / f"{prefix}.gottcha_species.bam"
                bam.write_bytes(source.read_bytes())
                Path(f"{bam}.bai").touch()
                with gzip.open(
                    outdir / f"{prefix}.sylph_extracted.fa.gz", "wt"
                ) as handle:
                    handle.write(f">reference_{source.name}\nACGT\n")
                if "--skip-qc" not in command:
                    (intermediate / f"{prefix}.qc.fastq.gz.json").write_text(
                        json.dumps(
                            {
                                "summary": {
                                    "before_filtering": {"total_reads": 10},
                                    "after_filtering": {"total_reads": 8},
                                }
                            }
                        ),
                        encoding="utf-8",
                    )
                (outdir / f"{prefix}.pathogen.full.tsv").write_text(
                    "LEVEL\tNAME\tTAXID\tREAD_COUNT\tSIG_COV\tSNI_SCORE\tHUMAN_PATHOGEN\n",
                    encoding="utf-8",
                )
                (outdir / f"{prefix}.pathogen.full.html").write_text(
                    "<html><body>No alignments</body></html>\n", encoding="utf-8"
                )
                (outdir / f"{prefix}.coverage.html").write_text(
                    "<html><body>No coverage</body></html>\n", encoding="utf-8"
                )
                (outdir / f"{prefix}.gottcha_species.coverage.tsv").write_text(
                    "#rname\tstartpos\tendpos\n", encoding="utf-8"
                )
            else:
                sequence = int(self.state["pending"]["timepoint"])
                (outdir / f"{prefix}.pathogen.full.tsv").write_text(
                    "LEVEL\tNAME\tTAXID\tREAD_COUNT\tSIG_COV\tSNI_SCORE\tHUMAN_PATHOGEN\n"
                    f"species\tTest pathogen\t1234\t{sequence * 10}\t"
                    f"{sequence / 10:.2f}\t{0.95 + sequence / 100:.3f}\tYes\n",
                    encoding="utf-8",
                )
                (outdir / f"{prefix}.pathogen.full.html").write_text(
                    f"<html><body>Batch {sequence} report</body></html>\n",
                    encoding="utf-8",
                )
                (outdir / f"{prefix}.coverage.html").write_text(
                    f"<html><body>Batch {sequence} coverage</body></html>\n",
                    encoding="utf-8",
                )
                intermediate = outdir / "intermediate"
                intermediate.mkdir()
                (intermediate / f"{prefix}.gottcha_species.coverage.tsv").write_text(
                    f"#rname\tstartpos\tendpos\nref\t1\t{sequence}\n",
                    encoding="utf-8",
                )
            return

        if command[1] == "merge":
            candidate = Path(command[5])
            inputs = [Path(item) for item in command[6:]]
            candidate.write_bytes(b"".join(item.read_bytes() for item in inputs))
        elif command[1] == "sort":
            source = Path(command[-1])
            destination = Path(command[command.index("-o") + 1])
            destination.write_bytes(source.read_bytes())
        elif command[1] == "index":
            bam = Path(command[-1])
            Path(f"{bam}.bai").touch()
        else:
            raise AssertionError(f"Unexpected command: {command}")


class StreamSpadesTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.input_dir = self.root / "incoming"
        self.output_dir = self.root / "output"
        self.input_dir.mkdir()
        self.run_spades = self.root / "run_SPADES.sh"
        self.run_spades.touch()
        self.database = self.root / "gottcha_db.species.fna"
        self.args = Namespace(
            input_dir=self.input_dir,
            outdir=self.output_dir,
            prefix="stream",
            db_path=self.database,
            cpu=2,
            run_spades=self.run_spades,
            run_spades_version="test-1.2.3",
            spades_data=None,
            samtools="samtools",
            poll_interval=0.0,
            settle_seconds=0.0,
            ont_error_rate=0.03,
            min_depth=10,
            skip_qc=False,
            recursive=False,
            keep_timepoints=False,
            js_external=False,
            once=True,
            max_files=None,
            verbose=False,
        )

    def tearDown(self):
        self.temporary.cleanup()

    def test_read_file_extensions(self):
        accepted = ["a.fa", "b.fasta.gz", "c.fna", "d.fq.gz", "e.fastq"]
        for name in accepted:
            path = self.input_dir / name
            path.touch()
            self.assertTrue(stream_spades.is_read_file(path), name)
        rejected = self.input_dir / "notes.txt"
        rejected.touch()
        self.assertFalse(stream_spades.is_read_file(rejected))

    def test_count_reads_supports_multiline_fasta_and_gzipped_fastq(self):
        fasta = self.input_dir / "records.fasta"
        fasta.write_bytes(b">first\nAC\nGT\n>second\nNN\n")
        fastq = self.input_dir / "records.fastq.gz"
        with gzip.open(fastq, "wb") as handle:
            handle.write(
                b"@first\nAC\nGT\n+\n!!!!\n"
                b"@second\nNN\n+second\n##\n"
            )

        self.assertEqual(stream_spades.count_reads(fasta), 2)
        self.assertEqual(stream_spades.count_reads(fastq), 2)

    def test_recursive_watch_is_enabled_by_default_and_can_be_disabled(self):
        parser = stream_spades.build_parser()
        required = [
            "--input-dir", str(self.input_dir),
            "--outdir", str(self.output_dir),
            "--prefix", "stream",
            "--db-path", str(self.database),
        ]

        self.assertTrue(parser.parse_args(required).recursive)
        self.assertTrue(parser.parse_args([*required, "--recursive"]).recursive)
        self.assertFalse(parser.parse_args([*required, "--no-recursive"]).recursive)

    def test_skip_qc_is_disabled_by_default_and_can_be_enabled(self):
        parser = stream_spades.build_parser()
        required = [
            "--input-dir", str(self.input_dir),
            "--outdir", str(self.output_dir),
            "--prefix", "stream",
            "--db-path", str(self.database),
        ]

        self.assertFalse(parser.parse_args(required).skip_qc)
        self.assertTrue(parser.parse_args([*required, "--skip-qc"]).skip_qc)

    def test_skip_qc_is_forwarded_and_read_counts_are_calculated(self):
        self.args.skip_qc = True
        input_path = self.input_dir / "chunk.fastq.gz"
        with gzip.open(input_path, "wb") as handle:
            handle.write(b"@first\nACGT\n+\n!!!!\n@second\nNN\n+\n##\n")

        runner = FakeStreamSpades(self.args)
        runner.run()

        record = runner.state["timepoints"][0]
        self.assertEqual(record["raw_reads"], 2)
        self.assertIsNone(record["filtered_reads"])
        self.assertEqual(record["cumulative_raw_reads"], 2)
        self.assertIsNone(record["cumulative_filtered_reads"])
        self.assertEqual(record["qc_json"], "")
        payload = build_payload(runner.state)
        self.assertEqual(payload["summary"]["raw_reads"], 2)
        self.assertIsNone(payload["summary"]["filtered_reads"])
        report = (self.output_dir / "stream.stream.html").read_text(encoding="utf-8")
        self.assertIn("formatNullable(payload.summary.filtered_reads)", report)
        self.assertIn("run_SPADES v{{ payload.run_spades_version", report)
        self.assertIn('"run_spades_version":"test-1.2.3"', report)

    def test_recursive_scan_finds_reads_in_subdirectories(self):
        self.args.recursive = True
        nested = self.input_dir / "run" / "pass"
        nested.mkdir(parents=True)
        read = nested / "chunk.fastq.gz"
        read.touch()

        runner = FakeStreamSpades(self.args)

        self.assertEqual(runner.scan(), [read.resolve()])

    def test_two_timepoints_merge_and_resume(self):
        first = self.input_dir / "chunk_1.fastq.gz"
        second = self.input_dir / "chunk_2.fasta"
        first.write_bytes(b"first\n")

        runner = FakeStreamSpades(self.args)
        runner.run()
        self.assertIn("Processing batch 000001", runner.processing_report)
        self.assertIn("chunk_1.fastq.gz", runner.processing_report)
        cumulative = self.output_dir / "cumulative" / "stream.gottcha_species.bam"
        self.assertEqual(cumulative.read_bytes(), b"first\n")
        cumulative_fasta = (
            self.output_dir / "cumulative" / "stream.sylph_extracted.fa.gz"
        )
        with gzip.open(cumulative_fasta, "rt") as handle:
            self.assertEqual(
                handle.read(), ">reference_chunk_1.fastq.gz\nACGT\n"
            )

        second.write_bytes(b"second\n")
        resumed = FakeStreamSpades(self.args)
        resumed.run()
        self.assertEqual(cumulative.read_bytes(), b"first\nsecond\n")
        with gzip.open(cumulative_fasta, "rt") as handle:
            merged_fasta = handle.read()
        self.assertIn(">reference_chunk_1.fastq.gz\nACGT\n", merged_fasta)
        self.assertIn(">reference_chunk_2.fasta\nACGT\n", merged_fasta)
        self.assertEqual(len(resumed.state["timepoints"]), 2)
        self.assertEqual(resumed.state["next_timepoint"], 3)

        manifest = (self.output_dir / "timepoints.tsv").read_text(encoding="utf-8")
        self.assertFalse((self.output_dir / "timepoints").exists())
        self.assertIn("chunk_1.fastq.gz", manifest)
        self.assertIn("chunk_2.fasta", manifest)
        self.assertIn("cumulative/stream.pathogen.full.tsv", manifest)
        self.assertIn("cumulative/stream.sylph_extracted.fa.gz", manifest)
        self.assertIn("raw_reads", manifest.splitlines()[0])
        self.assertEqual(resumed.state["timepoints"][0]["raw_reads"], 10)
        self.assertEqual(resumed.state["timepoints"][1]["filtered_reads"], 8)
        self.assertEqual(resumed.state["timepoints"][1]["cumulative_raw_reads"], 20)

        report = (self.output_dir / "stream.stream.html").read_text(encoding="utf-8")
        self.assertIn("Latest profiling results", report)
        self.assertIn("Analysis batches", report)
        self.assertIn("Post-QC reads", report)
        self.assertIn("Test pathogen", report)
        self.assertIn('"read_count":20', report)
        self.assertIn('"best_sig_cov":0.2', report)
        self.assertIn('"sni_score":0.97', report)
        self.assertIn("Taxa report", report)
        self.assertIn("Coverage", report)
        self.assertIn(
            '"pathogen_full_html":"cumulative/stream.pathogen.full.html"', report
        )
        self.assertIn(
            '"coverage_tsv":"cumulative/stream.gottcha_species.coverage.tsv"',
            report,
        )
        self.assertIn(
            '"coverage_html":"cumulative/stream.coverage.html"', report
        )
        self.assertTrue(
            (self.output_dir / "cumulative" / "stream.pathogen.full.html").is_file()
        )
        coverage = (
            self.output_dir / "cumulative" / "stream.gottcha_species.coverage.tsv"
        )
        self.assertIn("ref\t1\t2", coverage.read_text(encoding="utf-8"))

        payload = build_payload(resumed.state)
        self.assertEqual(payload["filter_defaults"]["levels"], ["species"])
        self.assertEqual(payload["summary"]["raw_reads"], 20)
        self.assertEqual(payload["summary"]["filtered_reads"], 16)
        self.assertEqual(len(payload["pathogens"]), 1)
        self.assertEqual(len(payload["pathogens"][0]["history"]), 2)
        self.assertEqual(
            [point["read_count"] for point in payload["pathogens"][0]["history"]],
            [10, 20],
        )

        no_duplicates = FakeStreamSpades(self.args)
        no_duplicates.run()
        self.assertEqual(len(no_duplicates.state["timepoints"]), 2)

    def test_timepoint_outputs_can_be_kept_explicitly(self):
        self.args.keep_timepoints = True
        (self.input_dir / "chunk.fastq.gz").write_bytes(b"reads\n")

        runner = FakeStreamSpades(self.args)
        runner.run()

        timepoint = self.output_dir / "timepoints" / "timepoint_000001"
        self.assertTrue((timepoint / "chunk").is_dir())
        self.assertTrue((timepoint / "profile" / "stream.pathogen.full.tsv").is_file())
        self.assertTrue(
            (
                timepoint
                / "profile"
                / "intermediate"
                / "stream.gottcha_species.coverage.tsv"
            ).is_file()
        )

    def test_alignment_free_chunk_publishes_placeholder_links(self):
        (self.input_dir / "empty.fastq.gz").write_bytes(b"reads\n")
        runner = FakeStreamSpades(self.args)
        runner._bam_is_valid = lambda path: False

        runner.run()

        record = runner.state["timepoints"][0]
        self.assertEqual(record["status"], "no_alignments")
        self.assertTrue(Path(record["pathogen_full_html"]).is_file())
        self.assertTrue(Path(record["coverage_tsv"]).is_file())
        report = (self.output_dir / "stream.stream.html").read_text(encoding="utf-8")
        self.assertIn("cumulative/stream.pathogen.full.html", report)
        self.assertIn("cumulative/stream.gottcha_species.coverage.tsv", report)

    def test_multiple_ready_files_are_processed_sequentially(self):
        first = self.input_dir / "01.fastq.gz"
        second = self.input_dir / "02.fastq.gz"
        third = self.input_dir / "03.fastq.gz"
        first.write_bytes(b"one\n")
        second.write_bytes(b"two\n")
        third.write_bytes(b"three\n")

        runner = FakeStreamSpades(self.args)
        runner.run()

        records = runner.state["timepoints"]
        self.assertEqual([record["timepoint"] for record in records], [1, 2, 3])
        self.assertEqual(
            [Path(record["input_file"]).name for record in records],
            ["01.fastq.gz", "02.fastq.gz", "03.fastq.gz"],
        )
        cumulative = self.output_dir / "cumulative" / "stream.gottcha_species.bam"
        self.assertEqual(cumulative.read_bytes(), b"one\ntwo\nthree\n")

    def test_running_pending_state_is_reported_as_current_activity(self):
        state = {
            "configuration": {"prefix": "stream", "database": "db"},
            "monitor_status": "running",
            "timepoints": [],
            "pending": {
                "timepoint": 1,
                "signature": {"path": "/incoming/chunk.fastq.gz"},
            },
        }
        payload = build_payload(state)
        self.assertEqual(payload["current"]["tone"], "active")
        self.assertIn("Processing batch 000001", payload["current"]["label"])
        self.assertEqual(payload["current"]["detail"], "chunk.fastq.gz")

    def test_latest_findings_are_distinguished_from_historical_findings(self):
        first_result = self.root / "first.pathogen.full.tsv"
        second_result = self.root / "second.pathogen.full.tsv"
        header = (
            "LEVEL\tNAME\tTAXID\tREAD_COUNT\tSIG_COV\tSNI_SCORE\tHUMAN_PATHOGEN\n"
        )
        first_result.write_text(
            header
            + "species\tHistorical organism\t111\t12\t0.2\t0.96\tYes\n",
            encoding="utf-8",
        )
        second_result.write_text(
            header + "species\tCurrent organism\t222\t24\t0.4\t0.98\tYes\n",
            encoding="utf-8",
        )
        state = {
            "configuration": {"prefix": "stream", "database": "db"},
            "monitor_status": "stopped",
            "pending": None,
            "timepoints": [
                {
                    "timepoint": 1,
                    "pathogen_full_tsv": str(first_result),
                    "status": "profiled",
                    "completed_at_utc": "2026-07-09T10:00:00+00:00",
                },
                {
                    "timepoint": 2,
                    "pathogen_full_tsv": str(second_result),
                    "status": "profiled",
                    "completed_at_utc": "2026-07-09T10:10:00+00:00",
                },
            ],
        }
        payload = build_payload(state)
        self.assertEqual(payload["summary"]["taxa"], 2)
        self.assertEqual(payload["summary"]["taxa_latest"], 1)
        self.assertEqual(payload["summary"]["human_pathogens"], 2)
        self.assertEqual(payload["summary"]["human_pathogens_latest"], 1)
        self.assertEqual(payload["pathogens"][0]["name"], "Current organism")
        self.assertTrue(payload["pathogens"][0]["present_latest"])
        self.assertFalse(payload["pathogens"][1]["present_latest"])


if __name__ == "__main__":
    unittest.main()
