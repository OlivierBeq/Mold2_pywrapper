# -*- coding: utf-8 -*-

import gc
import os
import unittest
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd

import Mold2_pywrapper
from Mold2_pywrapper.mold2_wrapper import convert_to_numeric, fix_bond_blocks


class TestHelperFunctions(unittest.TestCase):
    """Unit tests for standalone helper functions."""

    def test_convert_to_numeric_options(self):
        """Test all the fill_na options for the convert_to_numeric function."""
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
        """Test that fix_bond_blocks correctly raises an error for non-V2000 files."""
        fake_sdf_content = "Mol\n\n\n  0  0  0  0  0  0  0  0  0  0 V3000\n"
        with patch(
            "builtins.open", unittest.mock.mock_open(read_data=fake_sdf_content)
        ):
            with self.assertRaisesRegex(ValueError, "not a V2000 SD file"):
                fix_bond_blocks("dummy_path.sdf")


@patch("Mold2_pywrapper.mold2_wrapper.Mold2._download_executables", MagicMock())
class TestMold2Unit(unittest.TestCase):
    """
    Unit tests for the Mold2 class methods.
    All external dependencies (downloads, executables) are mocked.
    """

    @patch("Mold2_pywrapper.mold2_wrapper.Mold2._download_executables", MagicMock())
    def test_init_error_on_bad_fill_na(self):
        """Test ValueError for invalid fill_na arguments."""
        with self.assertRaisesRegex(ValueError, "Can only fill undefined values"):
            Mold2_pywrapper.Mold2(fill_na="invalid_string")

    @patch("Mold2_pywrapper.mold2_wrapper.Mold2._download_executables", MagicMock())
    def test_descriptor_detail_errors(self):
        """Test that descriptor_detail raises errors for out-of-bounds indices."""
        mold2 = Mold2_pywrapper.Mold2()
        with self.assertRaises(ValueError):
            mold2.descriptor_detail(0)
        with self.assertRaises(ValueError):
            mold2.descriptor_detail(778)

    @unittest.skipIf(os.name == "nt", "Parallel processing test is skipped on Windows")
    @patch("Mold2_pywrapper.mold2_wrapper.BoundedProcessPoolExecutor")
    @patch("Mold2_pywrapper.mold2_wrapper.Mold2._download_executables", MagicMock())
    def test_parallel_calculation_path(self, mock_executor):
        """Cover the parallel processing block (njobs > 1)."""
        # Configure the mock to simulate the entire parallel execution flow correctly
        mock_future = MagicMock()
        mock_future.result.return_value = (
            pd.DataFrame()
        )  # Ensure result() returns a DataFrame
        mock_executor.return_value.__enter__.return_value.submit.return_value = (
            mock_future
        )

        mold2 = Mold2_pywrapper.Mold2()
        # The show_banner=False is important to prevent stdout issues in some runners
        mold2.calculate(mols=[MagicMock()], njobs=2, show_banner=False)

        # Assert that our code *attempted* to use the process pool.
        mock_executor.assert_called_once()

    @patch("Mold2_pywrapper.mold2_wrapper.Mold2._download_executables", MagicMock())
    def test_lazy_loading_of_descriptor_details(self):
        """Cover the lazy-loading logic for descriptor details."""
        mold2 = Mold2_pywrapper.Mold2()
        self.assertFalse(hasattr(mold2, "_details"))
        mold2.descriptor_details()
        self.assertTrue(hasattr(mold2, "_details"))

    @patch("os.path.isdir", return_value=True)
    @patch("shutil.rmtree")
    @patch("Mold2_pywrapper.mold2_wrapper.Mold2._download_executables", MagicMock())
    def test_del_method_cleans_up(self, mock_isdir, mock_rmtree):
        """Cover the __del__ method for directory cleanup."""
        mold2 = Mold2_pywrapper.Mold2()
        mold2._dir = "/fake/temp/dir"
        del mold2
        # Force garbage collection to run, which will trigger __del__.
        gc.collect()
        mock_rmtree.assert_called_once_with("/fake/temp/dir")

    @patch("Mold2_pywrapper.mold2_wrapper.Mold2._download_executables", MagicMock())
    def test_parse_result_no_header(self):
        """Cover the logic for parsing a result file with no header."""
        mold2 = Mold2_pywrapper.Mold2()
        # Create mock DataFrame with object dtype to avoid FutureWarning
        mock_df = pd.DataFrame(np.random.rand(2, 777), dtype=object)
        mock_df.iloc[0, 0] = "This is not D001"
        with patch("pandas.read_table", return_value=mock_df):
            result = mold2._parse_result("dummy_path.txt")
            self.assertEqual(result.columns[0], "D001")

    @patch("Mold2_pywrapper.mold2_wrapper.platform", "linux")
    @patch("Mold2_pywrapper.mold2_wrapper.architecture", return_value=("32bit",))
    @patch("Mold2_pywrapper.mold2_wrapper.Mold2._download_executables", MagicMock())
    def test_prepare_command_on_linux32(self, mock_arch):
        """Cover the command preparation logic for 32-bit Linux."""
        mold2 = Mold2_pywrapper.Mold2()
        mold2._dir = "/fake/dir"
        with patch("os.path.isdir", return_value=True), patch(
            "os.path.isfile", return_value=True
        ):
            command = mold2._prepare_command("input.sdf", "output.txt")
            self.assertIn("Linux_x86-32", command)

    @patch("Mold2_pywrapper.mold2_wrapper.platform", "linux")
    @patch("os.chmod")
    def test_extract_executables_on_linux(self, mock_chmod):
        """Cover the os.chmod calls made on Linux during extraction."""
        mold2 = Mold2_pywrapper.Mold2()
        mold2._dir = "/fake/dir"
        with patch("zipfile.ZipFile"), patch("os.rename"):
            mold2._extract_executables("dummy.zip")
            self.assertGreater(mock_chmod.call_count, 0)

    def test_show_banner_path(self):
        """Cover the `show_banner=True` path in calculate()."""
        mold2 = Mold2_pywrapper.Mold2()
        # Mock the internal methods to isolate the banner logic
        with patch.object(mold2, "_show_banner") as mock_banner, patch.object(
            mold2, "_calculate", return_value=pd.DataFrame()
        ):
            mold2.calculate(mols=[], show_banner=True)
            mock_banner.assert_called_once()

    def test_del_method_no_dir(self):
        """Cover the __del__ path where _dir was never created."""
        mold2 = Mold2_pywrapper.Mold2()
        # The _dir attribute is not set, so __del__ should do nothing.
        # This test passes if no error is raised during garbage collection.
        del mold2
        gc.collect()

    def test_descriptor_detail_lazy_loading(self):
        """Cover the lazy-loading path specifically for descriptor_detail."""
        mold2 = Mold2_pywrapper.Mold2()
        self.assertFalse(hasattr(mold2, "_details"))
        # Calling descriptor_detail should trigger the lazy load.
        detail = mold2.descriptor_detail(1)
        self.assertTrue(hasattr(mold2, "_details"))
        self.assertIsNotNone(detail)

    @patch("shutil.copy")
    def test_from_executable_return_value(self, mock_copy):
        """Cover the return statement of the from_executable static method."""
        # We patch shutil.copy because we only care about the return value here.
        result = Mold2_pywrapper.Mold2.from_executable("dummy_path.zip")
        self.assertIsInstance(result, Mold2_pywrapper.Mold2)

    @patch("Mold2_pywrapper.mold2_wrapper.platform", "win32")
    def test_prepare_command_on_windows(self):
        """Cover the command preparation logic for Windows."""
        mold2 = Mold2_pywrapper.Mold2()
        mold2._dir = "/fake/dir"
        with patch("os.path.isdir", return_value=True), patch(
            "os.path.isfile", return_value=True
        ):
            command = mold2._prepare_command("input.sdf", "output.txt")
            self.assertIn("Mold2.exe", command)
            self.assertIn("echo.", command)

    def test_prepare_command_bad_paths(self):
        """Cover the error handling for non-existent paths."""
        mold2 = Mold2_pywrapper.Mold2()
        mold2._dir = "/fake/dir"
        # Test case where mold2_folder does not exist
        with patch("os.path.isdir", return_value=False):
            with self.assertRaises(RuntimeError):
                mold2._prepare_command("input.sdf", "output.txt")
        # Test case where input_path does not exist
        with patch("os.path.isdir", return_value=True), patch(
            "os.path.isfile", return_value=False
        ):
            with self.assertRaises(ValueError):
                mold2._prepare_command("input.sdf", "output.txt")


class TestMold2Download(unittest.TestCase):
    """A dedicated test case for the download logic."""

    @patch("requests.session")
    def test_download_path(self, mock_session):
        """Test that the download logic is called when the zip is missing."""
        mock_get = mock_session.return_value.get
        mock_get.return_value.raise_for_status.return_value = None
        mock_get.return_value.iter_content.return_value = [b"fake zip content"]

        with patch("os.path.isfile", return_value=False):
            # We must also patch the final `open` call to prevent a FileNotFoundError
            with patch("builtins.open", unittest.mock.mock_open()):
                Mold2_pywrapper.Mold2(verbose=True)
                mock_get.assert_called_once()
