# tests/test_unit.py

import unittest
from unittest.mock import patch, MagicMock
import pandas as pd
import numpy as np
import io

import Mold2_pywrapper
from Mold2_pywrapper.mold2_wrapper import convert_to_numeric, fix_bond_blocks


class TestHelperFunctions(unittest.TestCase):
    """Unit tests for standalone helper functions."""

    def test_convert_to_numeric_options(self):
        series = pd.Series(
            ["1.0", "2.0", None, "3.0", 3.4028199999999998e38], dtype=object
        )
        result_mean = convert_to_numeric(series.copy(), fill_value="mean")
        self.assertAlmostEqual(result_mean.iloc[2], 2.0)
        result_median = convert_to_numeric(series.copy(), fill_value="median")
        self.assertAlmostEqual(result_median.iloc[2], 2.0)
        result_float = convert_to_numeric(series.copy(), fill_value=-999.0)
        self.assertEqual(result_float.iloc[2], -999.0)
        result_nan = convert_to_numeric(series.copy(), fill_value="nan")
        self.assertTrue(np.isnan(result_nan.iloc[2]))

    def test_fix_bond_blocks_v3000_error(self):
        fake_sdf_content = "Mol\n\n\n  0  0  0  0  0  0  0  0  0  0 V3000\n"
        with patch(
            "builtins.open", unittest.mock.mock_open(read_data=fake_sdf_content)
        ):
            with self.assertRaisesRegex(ValueError, "not a V2000 SD file"):
                fix_bond_blocks("dummy_path.sdf")


class TestMold2Unit(unittest.TestCase):
    """Unit tests for the Mold2 class methods, with downloads completely mocked."""

    def setUp(self):
        # This setup ensures that _download_executables is mocked for every test in this class.
        patcher = patch(
            "Mold2_pywrapper.mold2_wrapper.Mold2._download_executables", MagicMock()
        )
        self.addCleanup(patcher.stop)
        patcher.start()

    def test_error_handling(self):
        """Test that the class raises errors on bad input."""
        mold2 = Mold2_pywrapper.Mold2(verbose=False)
        with self.assertRaises(ValueError):
            mold2.descriptor_detail(999)
        with self.assertRaises(ValueError):
            Mold2_pywrapper.Mold2(fill_na="bad_value")

    def test_descriptor_details(self):
        """Test the helper methods for retrieving descriptor information."""
        with patch("os.path.isfile", return_value=True):  # Pretend the zip file exists
            mold2 = Mold2_pywrapper.Mold2(verbose=False)
            detail = mold2.descriptor_detail(15)
            self.assertEqual(detail, "rotatable bond fraction")
            all_details = mold2.descriptor_details()
            self.assertIsInstance(all_details, dict)
            self.assertEqual(len(all_details), 777)

    @patch("requests.session")
    def test_download_path(self, mock_session):
        """Test that the download logic is called when the zip is missing."""
        mock_get = mock_session.return_value.get
        mock_get.return_value.raise_for_status.return_value = None
        mock_get.return_value.iter_content.return_value = [b"fake zip content"]

        with patch(
            "os.path.isfile", return_value=False
        ):  # Pretend zip file does NOT exist
            mold2 = Mold2_pywrapper.Mold2(verbose=False)
            mock_get.assert_called_once()

    @patch("shutil.rmtree")
    def test_del_method_cleans_up(self, mock_rmtree):
        mold2 = Mold2_pywrapper.Mold2()
        mold2._dir = "/fake/temp/dir"
        del mold2
        mock_rmtree.assert_called_once_with("/fake/temp/dir")

    def test_parse_result_no_header(self):
        mold2 = Mold2_pywrapper.Mold2()
        mock_df = pd.DataFrame(np.random.rand(2, 777))
        mock_df.iloc[0, 0] = "This is not D001"
        with patch("pandas.read_table", return_value=mock_df):
            result = mold2._parse_result("dummy_path.txt")
            self.assertEqual(result.columns[0], "D001")

    @patch("Mold2_pywrapper.mold2_wrapper.platform", "linux")
    @patch("Mold2_pywrapper.mold2_wrapper.architecture", return_value=("32bit",))
    def test_prepare_command_on_linux32(self, mock_arch, mock_platform):
        mold2 = Mold2_pywrapper.Mold2()
        mold2._dir = "/fake/dir"
        with patch("os.path.isdir", return_value=True), patch(
            "os.path.isfile", return_value=True
        ):
            command = mold2._prepare_command("input.sdf", "output.txt")
            self.assertIn("Linux_x86-32", command)
