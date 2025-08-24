# tests/test_integration.py

import os
import pathlib
import shutil
import tempfile
import unittest

from rdkit import Chem

import Mold2_pywrapper

# This environment variable is set by the CI workflow
DOWNLOADED_MOLD2_ZIP = os.environ.get("DOWNLOADED_MOLD2_ZIP_PATH")


@unittest.skipIf(
    not DOWNLOADED_MOLD2_ZIP,
    "Mold2 executable not pre-downloaded; skipping integration tests.",
)
class TestMold2Integration(unittest.TestCase):
    """
    An integration test suite that uses a real, pre-downloaded executable
    to test the end-to-end workflow.
    """

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.original_zipfile_path = Mold2_pywrapper.Mold2._zipfile
        self.addCleanup(
            setattr, Mold2_pywrapper.Mold2, "_zipfile", self.original_zipfile_path
        )
        self.patched_path = (
            pathlib.Path(self.temp_dir.name) / "Mold2" / "Mold2-Executable-File.zip"
        )
        Mold2_pywrapper.Mold2._zipfile = str(self.patched_path)
        os.makedirs(self.patched_path.parent, exist_ok=True)

    def test_end_to_end_calculation(self):
        """Test the main calculation workflow with njobs=1."""
        shutil.copy(DOWNLOADED_MOLD2_ZIP, self.patched_path)
        mold2 = Mold2_pywrapper.Mold2(verbose=False)
        mols = [Chem.MolFromSmiles(s) for s in ["CCO", "c1ccccc1"]]
        result_df = mold2.calculate(mols, show_banner=False, njobs=1)
        self.assertEqual(result_df.shape, (2, 777))
        self.assertEqual(result_df.loc[0, "D024"], 2)

    @unittest.skipIf(os.name == "nt", "Parallel processing test is skipped on Windows")
    def test_parallel_calculation(self):
        """Test the parallel processing path with njobs > 1."""
        shutil.copy(DOWNLOADED_MOLD2_ZIP, self.patched_path)
        mold2 = Mold2_pywrapper.Mold2(verbose=False)
        mols = [Chem.MolFromSmiles(s) for s in ["CCO", "c1ccccc1"]]
        result_df = mold2.calculate(mols, show_banner=False, njobs=2)
        self.assertEqual(result_df.shape, (2, 777))
        self.assertEqual(result_df.loc[0, "D024"], 2)

    def test_from_executable_installation(self):
        """Test installation from a user-provided zip."""
        self.assertFalse(self.patched_path.exists())
        Mold2_pywrapper.Mold2.from_executable(DOWNLOADED_MOLD2_ZIP)
        self.assertTrue(self.patched_path.exists())
