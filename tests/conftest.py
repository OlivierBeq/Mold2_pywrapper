import os
import shutil
from pathlib import Path

import pytest
from rdkit import Chem

from Mold2_pywrapper.mold2_wrapper import Mold2

DATA_DIR = Path(__file__).parent / "data"

# A real Mold2 executable ZIP, pointed to via this environment variable, either by a
# developer running tests locally or by CI (see .github/workflows/ci.yml). Downloaded
# from https://www.fda.gov/science-research/bioinformatics-tools/mold2
_ZIP_PATH_ENV_VAR = "DOWNLOADED_MOLD2_ZIP_PATH"


@pytest.fixture(scope="session")
def real_mold2_zip() -> str:
    """Path to a real Mold2 executables ZIP, or skip if none is available."""
    path = os.environ.get(_ZIP_PATH_ENV_VAR)
    if not path or not Path(path).is_file():
        pytest.skip(f"No real Mold2 executable ZIP available (set {_ZIP_PATH_ENV_VAR}).")
    return path


@pytest.fixture
def isolated_zipfile_path(tmp_path, monkeypatch):
    """Redirect Mold2's ZIP cache to an isolated temp dir for the duration of a test.

    Uses PYSTOW_HOME rather than monkeypatching Mold2._zipfile_path directly: worker
    processes spawned by ProcessPoolExecutor (njobs > 1) don't inherit in-memory
    monkeypatches on Windows (multiprocessing start method "spawn" re-imports the
    module fresh), but they do inherit the parent's environment variables.
    """
    monkeypatch.setenv("PYSTOW_HOME", str(tmp_path / "pystow-home"))
    return Mold2._zipfile_path()


@pytest.fixture
def mold2_real(isolated_zipfile_path, real_mold2_zip):
    """A real, working Mold2 instance backed by an isolated copy of the real executables."""
    os.makedirs(os.path.dirname(isolated_zipfile_path), exist_ok=True)
    shutil.copy(real_mold2_zip, isolated_zipfile_path)
    instance = Mold2(verbose=False)
    yield instance
    del instance


@pytest.fixture(scope="session")
def sample_smiles() -> list[str]:
    """~2000 diverse molecules sampled from a large real-world SMILES library."""
    return (DATA_DIR / "sample_molecules.smi").read_text().splitlines()


@pytest.fixture(scope="session")
def sample_mols(sample_smiles) -> list[Chem.Mol]:
    mols = [Chem.MolFromSmiles(smi) for smi in sample_smiles]
    return [mol for mol in mols if mol is not None]


@pytest.fixture(scope="session")
def small_mols(sample_mols) -> list[Chem.Mol]:
    """A small, fast-to-run subset, plus two molecules with known reference descriptor values."""
    known = [Chem.MolFromSmiles("CCO"), Chem.MolFromSmiles("c1ccccc1")]
    return known + sample_mols[:18]
