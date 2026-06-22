import os
import sys
import tempfile
import types
import unittest


if "pysam" not in sys.modules:
    sys.modules["pysam"] = types.ModuleType("pysam")

from gottcha.utils import profile


class TestProfileUtils(unittest.TestCase):
    def test_parse_args_defaults_for_short_reads(self):
        with tempfile.TemporaryDirectory() as tmp:
            db_prefix = os.path.join(tmp, "gottcha_db.species")
            read_path = os.path.join(tmp, "reads.fastq")
            open(db_prefix + ".mmi", "w").close()
            open(db_prefix + ".tax.tsv", "w").close()
            open(db_prefix + ".stats", "w").close()
            open(read_path, "w").close()

            args = profile.parse_args("test", ["profile", "-i", read_path, "-d", db_prefix])

            self.assertEqual(args.dbLevel, "species")
            self.assertEqual(args.matchIdentity, 0.95)
            self.assertEqual(args.errorRate, 0.005)
            self.assertEqual(args.prefix, "reads")
            self.assertEqual(args.input[0], os.path.abspath(read_path))

    def test_parse_args_nanopore_defaults(self):
        with tempfile.TemporaryDirectory() as tmp:
            db_prefix = os.path.join(tmp, "gottcha_db.species")
            read_path = os.path.join(tmp, "ont.fastq")
            open(db_prefix + ".mmi", "w").close()
            open(db_prefix + ".tax.tsv", "w").close()
            open(db_prefix + ".stats", "w").close()
            open(read_path, "w").close()

            args = profile.parse_args("test", ["profile", "-i", read_path, "-d", db_prefix, "-np"])

            self.assertEqual(args.matchIdentity, 0.85)
            self.assertEqual(args.matchFraction, 0.85)
            self.assertEqual(args.errorRate, 0.03)

    def test_parse_args_extractfullref_and_nocutoff(self):
        with tempfile.TemporaryDirectory() as tmp:
            db_prefix = os.path.join(tmp, "gottcha_db.species")
            read_path = os.path.join(tmp, "reads.fastq")
            open(db_prefix + ".mmi", "w").close()
            open(db_prefix + ".tax.tsv", "w").close()
            open(db_prefix + ".stats", "w").close()
            open(read_path, "w").close()

            args = profile.parse_args(
                "test",
                ["profile", "-i", read_path, "-d", db_prefix, "-ef", "-nc"],
            )
            self.assertEqual(args.extract, "all:20:fasta")
            self.assertEqual(args.sniScore, "0,0,0")

    def test_load_database_stats_with_header_row(self):
        with tempfile.TemporaryDirectory() as tmp:
            stats_path = os.path.join(tmp, "db.stats")
            with open(stats_path, "w") as f:
                f.write("Rank\tName\tTaxid\tSK\tNum\tMax\tMin\tTotalLength\tGenomeSize\tNote\n")
                f.write("species\tX\t123\tB\t1\t0\t0\t1000\t1200\tok\n")

            df = profile.load_database_stats(stats_path)
            self.assertIn("123", df.index)
            self.assertEqual(int(df.loc["123", "TotalLength"]), 1000)
            self.assertEqual(int(df.loc["123", "GenomeSize"]), 1200)

    def test_load_acc_list_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "acc.txt")
            open(p, "w").close()
            self.assertEqual(profile.load_acc_list(p), set())


if __name__ == "__main__":
    unittest.main()
