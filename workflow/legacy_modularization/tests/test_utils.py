import unittest
import pandas as pd
import os
import sys

# Add the project root to sys.path to allow imports from src
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.data_utils import (
    add_mutation_components,
    get_absolute_path,
    get_workflow_root,
    load_config,
    parse_mutation,
)

class TestUtils(unittest.TestCase):
    def test_parse_mutation(self):
        self.assertEqual(parse_mutation("F13S"), ("F", 13, "S"))
        self.assertEqual(parse_mutation("I3F"), ("I", 3, "F"))
        self.assertEqual(parse_mutation("A100T"), ("A", 100, "T"))
        self.assertEqual(parse_mutation("invalid"), (None, None, None))
        self.assertEqual(parse_mutation(None), (None, None, None))

    def test_add_mutation_components(self):
        df = pd.DataFrame({"Mutation": ["F13S", "I3F"]})
        df = add_mutation_components(df)
        self.assertIn("wild", df.columns)
        self.assertIn("position", df.columns)
        self.assertIn("mutant", df.columns)
        self.assertEqual(df.loc[0, "wild"], "F")
        self.assertEqual(df.loc[0, "position"], 13)
        self.assertEqual(df.loc[1, "mutant"], "F")

    def test_load_config(self):
        config = load_config()
        self.assertIn("project", config)
        self.assertIn("paths", config)
        self.assertEqual(config["project"]["name"], "cdkl5-variants")

    def test_workflow_relative_paths_resolve_inside_workflow(self):
        workflow_root = get_workflow_root()
        resolved = get_absolute_path("outputs/processed/cleaned_variants.xlsx")
        self.assertTrue(resolved.startswith(workflow_root))
        self.assertTrue(resolved.endswith("outputs/processed/cleaned_variants.xlsx"))

if __name__ == "__main__":
    unittest.main()
