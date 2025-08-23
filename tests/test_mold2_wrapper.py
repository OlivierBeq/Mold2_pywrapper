# -*- coding: utf-8 -*-


import unittest
from unittest.mock import patch

from rdkit import Chem
from Mold2_pywrapper import Mold2


class Mold2TestCase(unittest.TestCase):
    """
    Test suite for the Mold2 wrapper.
    These tests are designed to run without needing to download the executables.
    """

    @patch("src.Mold2_pywrapper.mold2_wrapper.Mold2._download_executables")
    def setUp(self, mock_download):
        """
        Set up a Mold2 instance for testing.
        The patch prevents the actual download during test runs.
        """
        # We mock the download, so the test will fail if the code can't run
        # without the executables. More advanced tests could mock the executable call itself.
        self.mold2 = Mold2()

    def test_instantiation(self):
        """
        Tests that the Mold2 object can be instantiated.
        """
        self.assertIsNotNone(self.mold2)

    @patch("src.Mold2_pywrapper.mold2_wrapper.Mold2._run_command")
    @patch("src.Mold2_pywrapper.mold2_wrapper.Mold2._prepare_input")
    def test_calculate_returns_dataframe(self, mock_prepare_input, mock_run_command):
        """
        Test that the calculate method returns a pandas DataFrame with the correct shape.
        This test mocks the actual calculation to run quickly and without the executables.
        """
        # Create a dummy RDKit molecule
        smiles = "CCO"  # Ethanol
        mol = Chem.MolFromSmiles(smiles)

        # We can't easily test the real calculation output without the executables,
        # so for this example, we will focus on mocking the internal calls.
        # A more advanced integration test would be needed to verify the actual results.

        # Let's assume a successful calculation would produce a result
        # We will mock the parser to return a dummy result.
        with patch(
            "src.Mold2_pywrapper.mold2_wrapper.Mold2._parse_result"
        ) as mock_parse:
            import pandas as pd
            import numpy as np

            # Create a fake DataFrame that _parse_result will return
            mock_df = pd.DataFrame(
                np.random.rand(1, 777), columns=[f"D{i:03d}" for i in range(1, 778)]
            )
            mock_parse.return_value = mock_df

            # Now, call the function we are testing
            result_df = self.mold2.calculate([mol], show_banner=False)

            # --- Assertions ---
            # 1. Check that the result is a pandas DataFrame
            self.assertIsInstance(result_df, pd.DataFrame)

            # 2. Check that the DataFrame has the correct shape (1 molecule, 777 descriptors)
            self.assertEqual(result_df.shape, (1, 777))

            # 3. Check that the column names are correct
            self.assertEqual(result_df.columns[0], "D001")
            self.assertEqual(result_df.columns[-1], "D777")
