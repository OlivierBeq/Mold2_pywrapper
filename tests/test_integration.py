import gc
import os
import shutil
import tempfile

import pytest
from rdkit import Chem

from Mold2_pywrapper.mold2_wrapper import Mold2

pytestmark = pytest.mark.integration


class TestEndToEnd:
    def test_known_molecules(self, mold2_real):
        mols = [Chem.MolFromSmiles(s) for s in ["CCO", "c1ccccc1"]]
        result = mold2_real.calculate(mols, show_banner=False, njobs=1)
        assert result.shape == (2, 777)
        assert result.loc[0, "D024"] == 2

    def test_small_batch(self, mold2_real, small_mols):
        result = mold2_real.calculate(small_mols, show_banner=False, njobs=1)
        assert result.shape == (len(small_mols), 777)
        assert list(result.columns[:2]) == ["D001", "D002"]
        # Values should be numeric (not left as raw strings)
        assert result.select_dtypes(include="number").shape[1] == 777

    def test_fill_na_options(self, isolated_zipfile_path, real_mold2_zip, small_mols):
        os.makedirs(os.path.dirname(isolated_zipfile_path), exist_ok=True)
        shutil.copy(real_mold2_zip, isolated_zipfile_path)
        mold2 = Mold2(fill_na="mean", verbose=False)
        result = mold2.calculate(small_mols, show_banner=False, njobs=1)
        assert not result.isna().any().any()

    def test_parallel_matches_sequential(self, mold2_real, small_mols):
        sequential = mold2_real.calculate(small_mols, show_banner=False, njobs=1)
        parallel = mold2_real.calculate(small_mols, show_banner=False, njobs=2, chunksize=5)
        assert sequential.shape == parallel.shape
        assert list(sequential.columns) == list(parallel.columns)

    def test_larger_sample(self, mold2_real, sample_mols):
        subset = sample_mols[:100]
        result = mold2_real.calculate(subset, show_banner=False, njobs=2, chunksize=25)
        assert result.shape == (len(subset), 777)

    def test_del_cleans_up_real_temp_dir(self, isolated_zipfile_path, real_mold2_zip, small_mols):
        # Built locally rather than via the `mold2_real` fixture: a fixture's return
        # value is kept alive by pytest's internal bookkeeping until teardown, so
        # `del` inside the test body would never drop the last reference.
        os.makedirs(os.path.dirname(isolated_zipfile_path), exist_ok=True)
        shutil.copy(real_mold2_zip, isolated_zipfile_path)
        mold2 = Mold2(verbose=False)
        mold2.calculate(small_mols, show_banner=False, njobs=1)
        real_dir = mold2._dir
        assert os.path.isdir(real_dir)
        del mold2
        gc.collect()
        assert not os.path.isdir(real_dir)


class TestFromExecutable:
    def test_install_and_use(self, isolated_zipfile_path, real_mold2_zip, small_mols):
        assert not os.path.isfile(isolated_zipfile_path)
        mold2 = Mold2.from_executable(real_mold2_zip)
        assert os.path.isfile(isolated_zipfile_path)
        result = mold2.calculate(small_mols[:5], show_banner=False, njobs=1)
        assert result.shape == (5, 777)


class TestExtractExecutables:
    def test_real_extraction_renames_and_permissions(self, mold2_real):
        mold2_real._dir = tempfile.mkdtemp(prefix="Mold2_test_")
        mold2_real._extract_executables(mold2_real._zipfile_path())
        base = os.path.join(mold2_real._dir, "extras", "Mold2")
        assert os.path.isfile(os.path.join(base, "Windows", "Mold2.exe"))
        assert os.path.isfile(os.path.join(base, "Windows", "Mold2.bat"))
        assert os.path.isfile(os.path.join(base, "Linux_x86-64", "Mold2"))
        mode = os.stat(os.path.join(base, "Linux_x86-64", "Mold2")).st_mode
        assert mode & 0o100  # executable by owner
