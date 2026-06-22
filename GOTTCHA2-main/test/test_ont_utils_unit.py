import gzip
import io
import os
import tempfile
import unittest
from types import SimpleNamespace

from gottcha.utils import ont_utils


class TestOntUtils(unittest.TestCase):
    def test_detect_format(self):
        self.assertEqual(ont_utils.detect_format(b">seq1\n"), "fasta")
        self.assertEqual(ont_utils.detect_format(b"@seq1\n"), "fastq")
        with self.assertRaises(ValueError):
            ont_utils.detect_format(b"")
        with self.assertRaises(ValueError):
            ont_utils.detect_format(b"seq1\n")

    def test_write_chunks_fasta_with_tail(self):
        out = io.BytesIO()
        written = ont_utils.write_chunks_fasta(
            out=out,
            rid=b"read1",
            seq=b"ACGTACGTAA",
            L=4,
            step=4,
            drop_tail=False,
            min_tail=2,
            prefix=b"np_",
        )
        self.assertEqual(written, 3)
        content = out.getvalue().decode()
        self.assertIn(">np_read1|chunk=0", content)
        self.assertIn(">np_read1|chunk=2", content)
        self.assertIn("ACGT", content)
        self.assertIn("AA", content)

    def test_split_to_fasta_fastq_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            inp = os.path.join(tmp, "reads.fastq")
            outp = os.path.join(tmp, "chunks.fasta")
            with open(inp, "wb") as f:
                f.write(b"@r1\nACGTAC\n+\nIIIIII\n")

            count = ont_utils.split_to_fasta(
                input_path=inp,
                output_path=outp,
                split_length=4,
                step_length=4,
                drop_tail=False,
                min_tail=2,
                prefix="x_",
            )
            self.assertEqual(count, 2)
            with open(outp, "rb") as f:
                data = f.read().decode()
            self.assertIn(">x_r1|chunk=0", data)
            self.assertIn(">x_r1|chunk=1", data)

    def test_preprocess_nanopore_reads(self):
        with tempfile.TemporaryDirectory() as tmp:
            inp = os.path.join(tmp, "ont.fasta")
            with open(inp, "wb") as f:
                f.write(b">r1\n" + b"A" * 300 + b"\n")

            out_reads = ont_utils.preprocess_nanopore_reads(
                reads=[inp],
                outdir=tmp,
                prefix="sample1",
                silent=True,
            )
            self.assertEqual(len(out_reads), 1)
            out_path = out_reads[0]
            self.assertTrue(out_path.endswith(".split_reads.fasta.gz"))
            self.assertTrue(os.path.exists(out_path))

            with gzip.open(out_path, "rt") as f:
                text = f.read()
            self.assertIn(">r1|chunk=0", text)
            self.assertIn(">r1|chunk=1", text)


if __name__ == "__main__":
    unittest.main()
