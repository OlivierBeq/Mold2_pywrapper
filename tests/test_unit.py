# -*- coding: utf-8 -*-

import unittest
from unittest.mock import patch
import pandas as pd
import numpy as np

# Import the module and the specific functions/classes needed to be tested
import Mold2_pywrapper
from Mold2_pywrapper.mold2_wrapper import convert_to_numeric, fix_bond_blocks


class TestHelperFunctions(unittest.TestCase):
    """Unit tests for standalone helper functions."""

    def test_convert_to_numeric_options(self):
        series = pd.Series(["1.0", "2.0", "NA", "3.0", 3.40282e38], dtype=str)
        result_mean = convert_to_numeric(series.copy(), fill_value="mean")
        self.assertAlmostEqual(result_mean[2], 2.0)
        result_median = convert_to_numeric(series.copy(), fill_value="median")
        self.assertAlmostEqual(result_median[2], 2.0)
        result_float = convert_to_numeric(series.copy(), fill_value=-999.0)
        self.assertEqual(result_float[2], -999.0)
        result_nan = convert_to_numeric(series.copy(), fill_value="nan")
        self.assertTrue(np.isnan(result_nan[2]))

    def test_fix_bond_blocks_error(self):
        fake_sdf_content = "Mol\n\n\n  0  0  0  0  0  0  0  0  0  0 V3000\n"
        with patch(
            "builtins.open", unittest.mock.mock_open(read_data=fake_sdf_content)
        ):
            with self.assertRaisesRegex(ValueError, "not a V2000 SD file"):
                fix_bond_blocks("dummy_path.sdf")


class TestMold2Unit(unittest.TestCase):
    """Unit tests for the Mold2 class methods, with downloads mocked."""

    # Patch the 'platform' object *where it is used* inside the mold2_wrapper module.
    @patch("Mold2_pywrapper.mold2_wrapper.platform", "darwin")
    @patch("Mold2_pywrapper.mold2_wrapper.Mold2._download_executables")
    def test_init_error_on_bad_platform(self, mock_download):
        """Test that instantiation fails on an unsupported platform like macOS."""
        with self.assertRaisesRegex(
            RuntimeError, "can only be calculated on Windows and Linux"
        ):
            Mold2_pywrapper.Mold2()

    @patch("Mold2_pywrapper.mold2_wrapper.Mold2._download_executables")
    def test_init_error_on_bad_fill_na(self, mock_download):
        """Test ValueError for invalid fill_na arguments."""
        with self.assertRaisesRegex(ValueError, "Can only fill undefined values"):
            Mold2_pywrapper.Mold2(fill_na="invalid_string")

    @patch("Mold2_pywrapper.mold2_wrapper.Mold2._download_executables")
    def test_descriptor_detail_errors(self, mock_download):
        """Test that descriptor_detail raises errors for out-of-bounds indices."""
        mold2 = Mold2_pywrapper.Mold2()
        with self.assertRaises(ValueError):
            mold2.descriptor_detail(0)
        with self.assertRaises(ValueError):
            mold2.descriptor_detail(778)

    @patch("Mold2_pywrapper.mold2_wrapper.Mold2._download_executables")
    def test_parse_result_no_header(self, mock_download):
        """Test parsing a result file where the header is missing."""
        mold2 = Mold2_pywrapper.Mold2()

        # The DataFrame the mock returns must have 777 columns, because that's what
        # the real _parse_result function would receive after `usecols=range(1, 778)`.
        mock_data = np.random.rand(2, 777)
        mock_df = pd.DataFrame(mock_data)
        # Give it a non-numeric first column to ensure the 'no header' logic triggers.
        mock_df.iloc[0, 0] = "This is not D001"

        with patch("pandas.read_table", return_value=mock_df) as mock_read_table:
            result = mold2._parse_result("dummy_path.txt")

            mock_read_table.assert_called_with(
                "dummy_path.txt", header=None, usecols=range(1, 778), low_memory=False
            )

            self.assertEqual(len(result.columns), 777)
            self.assertEqual(result.columns[0], "D001")
            self.assertEqual(result.columns[-1], "D777")
