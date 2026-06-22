import os
import sys
import tempfile
import types
import unittest
from collections import defaultdict

import pandas as pd


if "pysam" not in sys.modules:
    sys.modules["pysam"] = types.ModuleType("pysam")

from gottcha.utils import extract_reads


class TestExtractReadsUtils(unittest.TestCase):
    def setUp(self):
        extract_reads.ref_to_extract_taxid = defaultdict(list)

    def _res_df(self):
        return pd.DataFrame(
            [
                {
                    "LEVEL": "strain",
                    "NAME": "Alpha strain",
                    "SNI_SCORE": 0.995,
                    "TAXID": "111",
                    "PARENT_TAXID": "11",
                },
                {
                    "LEVEL": "species",
                    "NAME": "Alpha taxa",
                    "SNI_SCORE": 0.0,
                    "TAXID": "11",
                    "PARENT_TAXID": "1",
                },
                {
                    "LEVEL": "strain",
                    "NAME": "Beta strain",
                    "SNI_SCORE": 0.995,
                    "TAXID": "222",
                    "PARENT_TAXID": "22",
                },
                {
                    "LEVEL": "species",
                    "NAME": "Beta taxa",
                    "SNI_SCORE": 0.0,
                    "TAXID": "22",
                    "PARENT_TAXID": "1",
                },
                {
                    "LEVEL": "genus",
                    "NAME": "Genus 1",
                    "SNI_SCORE": 0.0,
                    "TAXID": "1",
                    "PARENT_TAXID": "2",
                },
                {
                    "LEVEL": "family",
                    "NAME": "Family 2",
                    "SNI_SCORE": 0.0,
                    "TAXID": "2",
                    "PARENT_TAXID": "3",
                },
                {
                    "LEVEL": "order",
                    "NAME": "Order 3",
                    "SNI_SCORE": 0.0,
                    "TAXID": "3",
                    "PARENT_TAXID": "4",
                },
                {
                    "LEVEL": "class",
                    "NAME": "Class 4",
                    "SNI_SCORE": 0.0,
                    "TAXID": "4",
                    "PARENT_TAXID": "5",
                },
                {
                    "LEVEL": "phylum",
                    "NAME": "Phylum 5",
                    "SNI_SCORE": 0.0,
                    "TAXID": "5",
                    "PARENT_TAXID": "6",
                },
                {
                    "LEVEL": "superkingdom",
                    "NAME": "Superkingdom 6",
                    "SNI_SCORE": 0.0,
                    "TAXID": "6",
                    "PARENT_TAXID": "0",
                },
            ]
        )

    def test_parse_taxids_filters_by_selected_taxon(self):
        res_df = self._res_df()
        taxa_dict, ref_to_extract = extract_reads.parse_taxids(
            "11",
            res_df,
            "unused.tsv",
            sni_score_cutoff=0.9,
            sni_score_species=0.95,
            sni_score_strain=0.99,
        )

        self.assertIn("11", taxa_dict)
        self.assertEqual(taxa_dict["11"]["LEVEL"], "species")
        self.assertEqual(taxa_dict["11"]["NAME"], "Alpha taxa")
        self.assertEqual(dict(ref_to_extract), {"111": ["11"]})

    def test_parse_taxids_from_file(self):
        res_df = self._res_df()
        with tempfile.TemporaryDirectory() as tmp:
            taxid_file = os.path.join(tmp, "taxids.txt")
            with open(taxid_file, "w") as f:
                f.write("# comment\n22\n")

            taxa_dict, ref_to_extract = extract_reads.parse_taxids(
                f"@{taxid_file}",
                res_df,
                "unused.tsv",
                sni_score_cutoff=0.9,
                sni_score_species=0.95,
                sni_score_strain=0.99,
            )

        self.assertIn("22", taxa_dict)
        self.assertEqual(dict(ref_to_extract), {"222": ["22"]})

    def test_iter_tasks_aoi_filters_and_lineage(self):
        extract_reads.ref_to_extract_taxid = defaultdict(list, {"111": ["11"], "222": []})
        refs = [
            "ACC1|1|100|111|A",
            "ACC2|1|100|111|A",
            "ACC3|1|100|222|A",
            "BADREF",
        ]

        tasks_out = list(extract_reads._iter_tasks(refs, 5, {"ACC2"}, "filter_out"))
        tasks_in = list(extract_reads._iter_tasks(refs, 5, {"ACC2"}, "filter_in"))

        self.assertEqual(tasks_out, [("ACC1|1|100|111|A", 5, False)])
        self.assertEqual(tasks_in, [("ACC2|1|100|111|A", 5, True)])


if __name__ == "__main__":
    unittest.main()
