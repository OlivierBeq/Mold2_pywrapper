# -*- coding: utf-8 -*-

import os
import shutil
import unittest
from unittest.mock import patch
import tempfile
import pathlib

import pandas as pd
from rdkit import Chem

import Mold2_pywrapper


# This environment variable will be set by the CI workflow
DOWNLOADED_MOLD2_ZIP = os.environ.get("DOWNLOADED_MOLD2_ZIP_PATH")


# If the pre-downloaded asset is not available (e.g., local run), skip these tests.
@unittest.skipIf(
    not DOWNLOADED_MOLD2_ZIP,
    "Mold2 executable not pre-downloaded; skipping integration tests.",
)
class TestMold2Integration(unittest.TestCase):
    """
    An integration test suite for the Mold2 wrapper.
    Each test uses a fully isolated environment by patching the external
    dependency `pystow.join`.
    """

    def setUp(self):
        """
        This method runs BEFORE every single test.
        It creates a new temporary directory, sets it as the cache home,
        and forces all relevant modules to re-read this new configuration.
        """
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.original_zipfile_path = Mold2_pywrapper.Mold2._zipfile
        self.addCleanup(setattr, Mold2_pywrapper.Mold2, '_zipfile', self.original_zipfile_path)

        self.patched_path = pathlib.Path(self.temp_dir.name) / "Mold2" / "Mold2-Executable-File.zip"
        Mold2_pywrapper.Mold2._zipfile = str(self.patched_path)
        os.makedirs(self.patched_path.parent, exist_ok=True)

    def test_end_to_end_calculation(self):
        """
        Test the full workflow: instantiation (with download), calculation, and result parsing.
        """
        shutil.copy(DOWNLOADED_MOLD2_ZIP, self.patched_path)
        # 3. Instantiate the calculator.
        mold2 = Mold2_pywrapper.Mold2(verbose=False)
        # 4. Prepare input molecules
        smiles_list = [
            "CCO",  # Ethanol
            "c1ccccc1",  # Benzene
        ]
        mols = [Chem.MolFromSmiles(s) for s in smiles_list]
        # 5. Perform the calculation
        result_df = mold2.calculate(mols, show_banner=False, njobs=1)
        # 4. Assert the results are correct
        self.assertIsInstance(result_df, pd.DataFrame)
        self.assertEqual(result_df.shape, (2, 777))
        self.assertTrue(pd.api.types.is_numeric_dtype(result_df["D001"]))
        # 6. Assert specific, known descriptor values for high confidence
        # For Ethanol (index 0):
        self.assertEqual(result_df.loc[0, "D024"], 2)  # Number of Carbons
        self.assertEqual(result_df.loc[0, "D026"], 1)  # Number of Oxygens
        self.assertEqual(result_df.loc[0, "D014"], 2)  # Number of rotatable bonds
        # For Benzene (index 1):
        self.assertEqual(
            result_df.loc[1, "D001"], 1
        )  # Number of 6-membered aromatic rings
        self.assertEqual(result_df.loc[1, "D024"], 6)  # Number of Carbons
        self.assertEqual(result_df.loc[1, "D017"], 6)  # Number of aromatic bonds
        self.assertEqual(result_df.loc[1, "D014"], 0)  # Number of rotatable bonds

    def test_descriptor_details(self):
        """
        Test the helper methods for retrieving descriptor information.
        """
        self.patched_path.touch() # Create dummy file to prevent download
        # Test retrieving a single descriptor
        mold2 = Mold2_pywrapper.Mold2(verbose=False)
        detail = mold2.descriptor_detail(15)
        self.assertEqual(detail, "rotatable bond fraction")
        # Test retrieving all descriptors
        all_details = mold2.descriptor_details()
        self.assertIsInstance(all_details, dict)
        self.assertEqual(len(all_details), 777)
        self.assertEqual(
            all_details["D001"],
            "number of 6-membered aromatic rings (only carbon atoms)",
        )

    def test_parallel_calculation(self):
        """Test the `njobs > 1` parallel processing path."""
        # Pre-seed the cache by copying our known-good zip file into place.
        shutil.copy(DOWNLOADED_MOLD2_ZIP, self.patched_path)
        mold2 = Mold2_pywrapper.Mold2(verbose=False)
        mols = [Chem.MolFromSmiles(s) for s in ["CCO", "c1ccccc1"]]
        # Call with njobs > 1 to engage the parallel code
        result_df = mold2.calculate(mols, show_banner=False, njobs=2)
        # Assert the results are correct
        self.assertIsInstance(result_df, pd.DataFrame)
        self.assertEqual(result_df.shape, (2, 777))
        self.assertTrue(pd.api.types.is_numeric_dtype(result_df["D001"]))
        # Assert specific, known descriptor values for high confidence
        # For Ethanol (index 0):
        self.assertEqual(result_df.loc[0, "D024"], 2)  # Number of Carbons
        self.assertEqual(result_df.loc[0, "D026"], 1)  # Number of Oxygens
        self.assertEqual(result_df.loc[0, "D014"], 2)  # Number of rotatable bonds
        # For Benzene (index 1):
        self.assertEqual(
            result_df.loc[1, "D001"], 1
        )  # Number of 6-membered aromatic rings
        self.assertEqual(result_df.loc[1, "D024"], 6)  # Number of Carbons
        self.assertEqual(result_df.loc[1, "D017"], 6)  # Number of aromatic bonds
        self.assertEqual(result_df.loc[1, "D014"], 0)  # Number of rotatable bonds

    def test_from_executable_installation(self):
        """
        Test that `from_executable` correctly copies a user-provided ZIP file
        into its own clean environment.
        """
        self.assertFalse(self.patched_path.exists())
        # 2. The source file is our known-good, pre-downloaded zip
        Mold2_pywrapper.Mold2.from_executable(DOWNLOADED_MOLD2_ZIP)
        self.assertTrue(self.patched_path.exists())

    @patch('requests.session')
    def test_download_path(self, mock_session):
        """Test that the download logic is called when the zip is missing."""
        # We don't pre-seed the cache, so os.path.isfile will be false.
        self.assertFalse(self.patched_path.exists())
        # Configure the mock to simulate a successful download
        mock_get = mock_session.return_value.get
        mock_get.return_value.raise_for_status.return_value = None
        # Provide a real zip file's content to the mock
        with open(DOWNLOADED_MOLD2_ZIP, 'rb') as f:
            mock_get.return_value.iter_content.return_value = [f.read()]
        # This will now trigger the download block inside __init__
        mold2 = Mold2_pywrapper.Mold2(verbose=False)
        # Assert that the download was attempted and the file was created
        mock_get.assert_called_once()
        self.assertTrue(self.patched_path.exists())

    def test_error_handling(self):
        """Test that the class raises errors on bad input."""
        self.patched_path.touch()  # Create dummy file to prevent download
        mold2 = Mold2_pywrapper.Mold2(verbose=False)
        # Test bad descriptor index
        with self.assertRaises(ValueError):
            mold2.descriptor_detail(999)
        # Test bad fill_na value
        with self.assertRaises(ValueError):
            Mold2_pywrapper.Mold2(fill_na="bad_value")
