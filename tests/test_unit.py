import gc
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem

import Mold2_pywrapper
from Mold2_pywrapper.mold2_wrapper import (
    Mold2,
    _get_user_agent,
    convert_to_numeric,
    fix_bond_blocks,
    mktempfile,
)

DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture
def mold2_mocked(monkeypatch):
    """A Mold2 instance with the download step disabled, for logic that doesn't need a real exe."""
    monkeypatch.setattr(Mold2, "_download_executables", lambda self, verbose=True: None)
    return Mold2()


class TestConvertToNumeric:
    SERIES = ["1.0", "2.0", None, "3.0", 3.4028199999999998e38]

    @pytest.mark.parametrize(
        "fill_value, expected",
        [
            ("mean", 2.0),
            ("MEAN", 2.0),
            ("median", 2.0),
            (-999.0, -999.0),
            ("nan", np.nan),
        ],
    )
    def test_fill_options(self, fill_value, expected):
        series = pd.Series(self.SERIES, dtype=object)
        result = convert_to_numeric(series, fill_value=fill_value)
        if np.isnan(expected):
            assert np.isnan(result.iloc[2])
            assert np.isnan(result.iloc[4])  # the Mold2 "undefined" sentinel is also replaced
        else:
            assert result.iloc[2] == pytest.approx(expected)
            assert result.iloc[4] == pytest.approx(expected)

    def test_no_fill_leaves_nans(self):
        series = pd.Series(self.SERIES, dtype=object)
        result = convert_to_numeric(series, fill_value=None)
        assert np.isnan(result.iloc[2])
        assert result.iloc[4] == pytest.approx(3.4028199999999998e38)

    def test_numeric_values_preserved(self):
        series = pd.Series(self.SERIES, dtype=object)
        result = convert_to_numeric(series, fill_value=None)
        assert result.iloc[0] == pytest.approx(1.0)
        assert result.iloc[1] == pytest.approx(2.0)
        assert result.iloc[3] == pytest.approx(3.0)


class TestFixBondBlocks:
    def test_raises_on_v3000(self, tmp_path):
        path = tmp_path / "mol.sdf"
        path.write_text("Mol\n\n\n  0  0  0  0  0  0  0  0  0  0 V3000\n")
        with pytest.raises(ValueError, match="not a V2000 SD file"):
            fix_bond_blocks(str(path))

    def test_fixes_bond_lines_on_v2000(self, tmp_path):
        path = tmp_path / "mol.sdf"
        # A minimal, RDKit-style V2000 molblock: header (4 lines), 1 atom, 1 bond line
        # missing the trailing 3 fields RDKit does not write.
        content = (
            "\n     RDKit          2D\n\n"
            "  1  0  0  0  0  0  0  0  0  0999 V2000\n"
            "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "  1  2  1  0\n"
            "M  END\n"
        )
        path.write_text(content)
        fix_bond_blocks(str(path))
        lines = path.read_text().splitlines()
        assert lines[5] == "  1  2  1  0  0  0  0"

    def test_non_bond_lines_untouched(self, tmp_path):
        path = tmp_path / "mol.sdf"
        content = "\n     RDKit          2D\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n"
        path.write_text(content)
        fix_bond_blocks(str(path))
        assert path.read_text() == content


class TestMktempfile:
    def test_creates_writeable_file_with_suffix(self):
        path = mktempfile(".sdf")
        try:
            assert os.path.isfile(path)
            assert path.endswith(".sdf")
            with open(path, "wt") as fh:
                fh.write("hello")
        finally:
            os.remove(path)


class TestUserAgent:
    def test_returns_a_string(self):
        ua = _get_user_agent()
        assert isinstance(ua, str)
        assert "Mozilla" in ua

    def test_falls_back_on_error(self, monkeypatch):
        class BoomUserAgent:
            def __init__(self, *args, **kwargs):
                raise RuntimeError("no data available")

        monkeypatch.setattr("Mold2_pywrapper.mold2_wrapper.UserAgent", BoomUserAgent)
        ua = _get_user_agent()
        assert ua == Mold2_pywrapper.mold2_wrapper._FALLBACK_USER_AGENT


class TestInit:
    def test_rejects_bad_fill_na(self, monkeypatch):
        monkeypatch.setattr(Mold2, "_download_executables", lambda self, verbose=True: None)
        with pytest.raises(ValueError, match="Can only fill undefined values"):
            Mold2(fill_na="invalid_string")

    @pytest.mark.parametrize("fill_na", ["mean", "MEDIAN", "NaN", 0.0, -1])
    def test_accepts_valid_fill_na(self, monkeypatch, fill_na):
        monkeypatch.setattr(Mold2, "_download_executables", lambda self, verbose=True: None)
        mold2 = Mold2(fill_na=fill_na)
        assert mold2.fill_na == fill_na

    def test_rejects_unsupported_platform(self, monkeypatch):
        monkeypatch.setattr("Mold2_pywrapper.mold2_wrapper.sys.platform", "darwin")
        with pytest.raises(RuntimeError, match="Windows and Linux platforms"):
            Mold2()


class TestDescriptorDetails:
    def test_out_of_bounds_indices_raise(self, mold2_mocked):
        with pytest.raises(ValueError):
            mold2_mocked.descriptor_detail(0)
        with pytest.raises(ValueError):
            mold2_mocked.descriptor_detail(778)

    def test_lazy_loading(self, mold2_mocked):
        assert not hasattr(mold2_mocked, "_details")
        mold2_mocked.descriptor_details()
        assert hasattr(mold2_mocked, "_details")

    def test_lazy_loading_via_descriptor_detail(self, mold2_mocked):
        assert not hasattr(mold2_mocked, "_details")
        mold2_mocked.descriptor_detail(1)
        assert hasattr(mold2_mocked, "_details")

    def test_content_matches_ground_truth(self, mold2_mocked):
        descriptors_path = DATA_DIR.parent.parent / "src" / "Mold2_pywrapper" / "descriptors.json"
        ground_truth = json.loads(descriptors_path.read_text())
        details = mold2_mocked.descriptor_details()
        assert details == ground_truth
        assert mold2_mocked.descriptor_detail(1) == ground_truth["D001"]
        assert mold2_mocked.descriptor_detail(777) == ground_truth["D777"]


class TestDel:
    # NOTE: these tests build their own Mold2 instance instead of using the
    # `mold2_mocked` fixture: a fixture's return value is kept alive by pytest's
    # internal bookkeeping until teardown, so `del` inside the test body would
    # never drop the last reference and __del__ would never fire.

    def test_cleans_up_real_directory(self, monkeypatch, tmp_path):
        monkeypatch.setattr(Mold2, "_download_executables", lambda self, verbose=True: None)
        mold2 = Mold2()
        real_dir = tmp_path / "Mold2_instance"
        real_dir.mkdir()
        (real_dir / "marker.txt").write_text("x")
        mold2._dir = str(real_dir)
        assert real_dir.is_dir()
        del mold2
        gc.collect()
        assert not real_dir.is_dir()

    def test_noop_without_dir(self, monkeypatch):
        monkeypatch.setattr(Mold2, "_download_executables", lambda self, verbose=True: None)
        mold2 = Mold2()
        # _dir was never set; __del__ should do nothing and raise nothing.
        del mold2
        gc.collect()

    def test_survives_removal_error(self, monkeypatch, tmp_path):
        monkeypatch.setattr(Mold2, "_download_executables", lambda self, verbose=True: None)
        mold2 = Mold2()
        real_dir = tmp_path / "Mold2_instance"
        real_dir.mkdir()
        mold2._dir = str(real_dir)

        def boom(path):
            raise OSError("locked")

        monkeypatch.setattr("shutil.rmtree", boom)
        del mold2
        gc.collect()  # must not raise


class TestPrepareCommand:
    def test_missing_mold2_folder_raises(self, mold2_mocked, tmp_path):
        mold2_mocked._dir = str(tmp_path)
        with pytest.raises(RuntimeError, match="Could not locate Mold2 executables"):
            mold2_mocked._prepare_command("input.sdf", "output.txt")

    def test_missing_input_file_raises(self, mold2_mocked, tmp_path):
        mold2_mocked._dir = str(tmp_path)
        os.makedirs(tmp_path / "extras")
        with pytest.raises(ValueError, match="Input file does not exist"):
            mold2_mocked._prepare_command("nope.sdf", "output.txt")

    def test_missing_output_dir_raises(self, mold2_mocked, tmp_path):
        mold2_mocked._dir = str(tmp_path)
        os.makedirs(tmp_path / "extras")
        input_file = tmp_path / "input.sdf"
        input_file.write_text("x")
        with pytest.raises(ValueError, match="Path to output file does not exist"):
            mold2_mocked._prepare_command(str(input_file), "/no/such/dir/output.txt")

    def test_windows_command(self, mold2_mocked, tmp_path, monkeypatch):
        mold2_mocked._dir = str(tmp_path)
        os.makedirs(tmp_path / "extras")
        input_file = tmp_path / "input.sdf"
        input_file.write_text("x")
        monkeypatch.setattr("Mold2_pywrapper.mold2_wrapper.sys.platform", "win32")
        command = mold2_mocked._prepare_command(str(input_file), str(tmp_path / "output.txt"))
        assert command[0].endswith("Mold2.exe")
        assert command[-1] == "NUL"

    def test_linux64_command(self, mold2_mocked, tmp_path, monkeypatch):
        mold2_mocked._dir = str(tmp_path)
        os.makedirs(tmp_path / "extras")
        input_file = tmp_path / "input.sdf"
        input_file.write_text("x")
        monkeypatch.setattr("Mold2_pywrapper.mold2_wrapper.sys.platform", "linux")
        monkeypatch.setattr("Mold2_pywrapper.mold2_wrapper.architecture", lambda: ("64bit", ""))
        command = mold2_mocked._prepare_command(str(input_file), str(tmp_path / "output.txt"))
        assert "Linux_x86-64" in command[0]
        assert command[-1] == "/dev/null"

    def test_linux32_command(self, mold2_mocked, tmp_path, monkeypatch):
        mold2_mocked._dir = str(tmp_path)
        os.makedirs(tmp_path / "extras")
        input_file = tmp_path / "input.sdf"
        input_file.write_text("x")
        monkeypatch.setattr("Mold2_pywrapper.mold2_wrapper.sys.platform", "linux")
        monkeypatch.setattr("Mold2_pywrapper.mold2_wrapper.architecture", lambda: ("32bit", ""))
        command = mold2_mocked._prepare_command(str(input_file), str(tmp_path / "output.txt"))
        assert "Linux_x86-32" in command[0]

    def test_unsupported_platform_raises(self, mold2_mocked, tmp_path, monkeypatch):
        mold2_mocked._dir = str(tmp_path)
        os.makedirs(tmp_path / "extras")
        input_file = tmp_path / "input.sdf"
        input_file.write_text("x")
        monkeypatch.setattr("Mold2_pywrapper.mold2_wrapper.sys.platform", "darwin")
        with pytest.raises(RuntimeError, match="not supported"):
            mold2_mocked._prepare_command(str(input_file), str(tmp_path / "output.txt"))


class TestRunCommand:
    def test_raises_on_nonzero_exit(self, mold2_mocked):
        with pytest.raises(Exception):
            mold2_mocked._run_command(["python3", "-c", "import sys; sys.exit(1)"])

    def test_succeeds_on_zero_exit(self, mold2_mocked):
        mold2_mocked._run_command(["python3", "-c", "pass"])


class TestParseResult:
    def test_header_present(self, mold2_mocked, tmp_path):
        path = tmp_path / "out.txt"
        header = "\t".join(["Name"] + [f"D{i:03d}" for i in range(1, 778)])
        row = "\t".join(["mol1"] + ["1.0"] * 777)
        path.write_text(header + "\n" + row + "\n")
        result = mold2_mocked._parse_result(str(path))
        assert list(result.columns[:3]) == ["D001", "D002", "D003"]
        assert len(result) == 1

    def test_header_absent(self, mold2_mocked, tmp_path):
        path = tmp_path / "out.txt"
        row = "\t".join(["mol1"] + ["1.0"] * 777)
        path.write_text(row + "\n")
        result = mold2_mocked._parse_result(str(path))
        assert list(result.columns[:3]) == ["D001", "D002", "D003"]
        assert len(result) == 1


class TestPrepareInput:
    def test_writes_v2000_sdf(self, mold2_mocked):
        mols = [Chem.MolFromSmiles("CCO"), Chem.MolFromSmiles("c1ccccc1")]
        path = mold2_mocked._prepare_input(mols)
        try:
            content = Path(path).read_text()
            assert "V2000" in content
            assert content.count("$$$$") == 2
        finally:
            os.remove(path)


class TestDownloadExecutables:
    def test_skips_when_already_downloaded(self, tmp_path, monkeypatch):
        # Bypass __init__ and the mold2_mocked fixture entirely: the fixture's own
        # monkeypatch of _download_executables to a no-op would otherwise shadow the
        # real method this test needs to exercise, and would make the assertion below
        # pass vacuously regardless of whether the skip logic actually works. Building
        # the instance directly also guarantees the real default cache is never touched.
        mold2 = Mold2.__new__(Mold2)
        existing = tmp_path / "Mold2-Executable-File.zip"
        existing.write_bytes(b"already here")
        monkeypatch.setattr(Mold2, "_zipfile_path", staticmethod(lambda: str(existing)))
        called = []
        monkeypatch.setattr("requests.session", lambda: called.append(True))
        mold2._download_executables(verbose=False)
        assert not called

    def test_downloads_when_missing(self, tmp_path, monkeypatch):
        # Bypass __init__ (and its own call to _download_executables) entirely: this
        # test targets the real _download_executables implementation directly.
        mold2 = Mold2.__new__(Mold2)
        target = tmp_path / "sub" / "Mold2-Executable-File.zip"
        monkeypatch.setattr(Mold2, "_zipfile_path", staticmethod(lambda: str(target)))

        class FakeResponse:
            def raise_for_status(self):
                pass

            def iter_content(self, chunk_size):
                return [b"fake", b"zip", b"content"]

        class FakeSession:
            def get(self, *args, **kwargs):
                return FakeResponse()

        monkeypatch.setattr("requests.session", FakeSession)
        mold2._download_executables(verbose=True)
        assert target.is_file()
        assert target.read_bytes() == b"fakezipcontent"


class TestZipfilePath:
    def test_real_default_path(self, tmp_path, monkeypatch):
        # Redirect pystow's home dir so this never touches (or creates) the real
        # default cache directory, and can never see a real cached executable.
        monkeypatch.setenv("PYSTOW_HOME", str(tmp_path))
        path = Mold2._zipfile_path()
        assert path.endswith("Mold2-Executable-File.zip")
        assert "Mold2" in path
        assert path.startswith(str(tmp_path))


class TestShowBanner:
    def test_prints_citation_notice(self, mold2_mocked, monkeypatch, capsys):
        mold2_mocked._show_banner()
        captured = capsys.readouterr()
        assert "Mold2" in captured.out
        assert "DOI: 10.1021/ci800038f" in captured.out

    def test_calculate_shows_banner_by_default(self, mold2_mocked, monkeypatch, capsys):
        monkeypatch.setattr(mold2_mocked, "_calculate", lambda mols: pd.DataFrame())
        mold2_mocked.calculate(mols=[], show_banner=True)
        captured = capsys.readouterr()
        assert "Mold2" in captured.out


class TestFromExecutable:
    def test_copies_zip_to_default_location(self, tmp_path, monkeypatch):
        monkeypatch.setattr(Mold2, "_download_executables", lambda self, verbose=True: None)
        source = tmp_path / "source.zip"
        source.write_bytes(b"zip-bytes")
        destination = tmp_path / "nested" / "Mold2-Executable-File.zip"
        monkeypatch.setattr(Mold2, "_zipfile_path", staticmethod(lambda: str(destination)))

        result = Mold2.from_executable(str(source))

        assert isinstance(result, Mold2)
        assert destination.is_file()
        assert destination.read_bytes() == b"zip-bytes"
