import os
import sys
import unittest

import pandas as pd

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.data_cleaning import (
    check_reference_match,
    convert_three_letter_change,
    standardize_hector_rows,
    standardize_kgp_rows,
)


class TestDataCleaning(unittest.TestCase):
    def test_convert_three_letter_change(self):
        self.assertEqual(convert_three_letter_change("Ile3Phe"), "I3F")
        self.assertEqual(convert_three_letter_change("Tyr24Cys"), "Y24C")
        self.assertEqual(convert_three_letter_change("bad-value"), None)

    def test_standardize_kgp_rows(self):
        df = pd.DataFrame(
            {
                "HGVSp": ["ENSP00000369325.3:p.Ile3Phe"],
                "Protein change": ["I3F"],
                "CHROM": ["chrX"],
                "POS": [18507103],
            }
        )
        standardized = standardize_kgp_rows(df)
        self.assertEqual(standardized.loc[0, "Gene(s)"], "CDKL5")
        self.assertEqual(standardized.loc[0, "Condition(s)"], "Healthy")
        self.assertEqual(standardized.loc[0, "Source"], "The 1000 Genomes Project")
        self.assertEqual(standardized.loc[0, "GRCh38Chromosome"], "X")
        self.assertEqual(standardized.loc[0, "Germline classification"], "Benign")

    def test_standardize_hector_rows(self):
        df = pd.DataFrame(
            {
                "Mutation": ["H36R"],
                "Consequence": ["Benign/Likely benign"],
                "Source": ["PMID_29264392"],
            }
        )
        standardized = standardize_hector_rows(df)
        self.assertEqual(standardized.loc[0, "Name"], "PMID_29264392")
        self.assertEqual(standardized.loc[0, "Protein change"], "H36R")
        self.assertEqual(standardized.loc[0, "Germline classification"], "Benign/Likely benign")
        self.assertEqual(standardized.loc[0, "Source"], "PMID_29264392")

    def test_check_reference_match(self):
        sequence = "MIN"
        self.assertEqual(check_reference_match("M", 1, sequence), "Match")
        self.assertEqual(check_reference_match("A", 1, sequence), "Mismatch (Seq: M)")
        self.assertEqual(check_reference_match("M", 10, sequence), "Invalid")


if __name__ == "__main__":
    unittest.main()
