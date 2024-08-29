# -*- coding: utf-8

"""Python wrapper for Mold2 descriptors"""

from __future__ import annotations

import json
import os
import re
import shutil
import stat
import subprocess
import tempfile
import zipfile
from platform import architecture
from sys import platform
from typing import Iterable, List, Optional, Union

import more_itertools
import numpy as np
import pandas as pd
import pystow
import requests
from bounded_pool_executor import BoundedProcessPoolExecutor
from rdkit import Chem
from rdkit.Chem import AllChem


class Mold2:
    """Mold2 wrapper to obtain molecular descriptors."""

    # Default folder for Mold2 executables
    _zipfile = os.path.abspath(os.path.join(pystow.join('Mold2').as_posix(),
                                            'Mold2-Executable-File.zip'))

    def __init__(self, fill_na: Union[str, float] = None, verbose: bool = True):
        """Instantiate a wrapper to calculate Mold2 molecular descriptors.

        :param fill_na: value to fill undefined values by; may be a number or one of {'median', 'mean', 'NaN'}
        :param verbose: Should details about the download of executables be printed out
        """
        if platform not in ['win32', 'linux']:
            raise RuntimeError(f'Mold2 descriptors can only be calculated on Windows and Linux platforms.')
        if fill_na is not None and (not np.isreal(fill_na) and fill_na not in ['median', 'mean', 'NaN']):
            raise ValueError("Can only fill undefined values by a float or one of {'median', 'mean', 'NaN'}")
        self.fill_na = fill_na
        # Ensure executables are available
        self._download_executables(verbose)

    def __del__(self):
        """Remove downloaded executables."""
        if hasattr(self, '_dir') and os.path.isdir(self._dir):
            shutil.rmtree(self._dir)

    def calculate(self, mols: Iterable[Chem.Mol], show_banner: bool = True, njobs: int = 1,
                  chunksize: Optional[int] = 1000) -> pd.DataFrame:
        """Caclulate Mold2 descriptors.

        :param mols: RDkit molecules for which Mold2 descriptors should be calculated
        :param show_banner: If True, show notice on Mold2 descriptors usage
        :param njobs: number of concurrent processes
        :param chunksize: number of molecules to be processed by a process; ignored if njobs is 1
        :return: a pandas DataFrame containing all Mold2 descriptor values
        """
        if show_banner:
            self._show_banner()
        # Parallelize should need be
        if njobs > 1:
            with BoundedProcessPoolExecutor(max_workers=njobs) as worker:
                futures = [worker.submit(self._calculate, list(chunk))
                           for chunk in more_itertools.batched(mols, chunksize)
                           ]
            return pd.concat([future.result()
                              for future in futures]
                             ).reset_index(drop=True)
        # Single process
        return self._calculate(list(mols))

    def descriptor_detail(self, index: int) -> str:
        """Details of a Mold2 descriptor."""
        if not (0 < index < 778):
            raise ValueError('Mold2 descriptor index must be strictly greater than 0 and strictly less than 778.')
        # Lazy loading to avoid too large memory footprint
        if not hasattr(self, '_details'):
            self._parse_details()
        return self._details[f'D{index:03d}']

    def descriptor_details(self):
        # Lazy loading to avoid too large memory footprint
        if not hasattr(self, '_details'):
            self._parse_details()
        return self._details

    @staticmethod
    def from_executable(zipfile_path: str) -> Mold2:
        """Instantiate a Mold2 object from the user-downloaded mold2 zip file containing binaries.

        The provided ZIP file is extracted, so that default instantiation of Mold2 is henceforth possible.
        The ZIP file must be downloaded from "https://www.fda.gov/science-research/bioinformatics-tools/mold2"

        :param zipfile_path: Path to the zip file containing Mold2 binaries
        :return: a Mold2 calculator object
        """
        shutil.copy(zipfile_path, Mold2._zipfile)
        return Mold2()

    def _show_banner(self):
        """Print info message for citing."""
        print("""Mold2 calculates a large and diverse set of molecular descriptors encoding two-
dimensional chemical structure information. Comparative analysis of Mold2 descriptors
with those calculated from commercial software on several published datasets
demonstrated that Mold2 descriptors convey sufficient structural information. In addition,
better models were generated using Mold2 descriptors than the compared commercial
software packages. This publicly available software is developed by the Center for
Bioinformatics, which is led by Dr. Weida Tong, at the National Center for Toxicological
Research (NCTR).
    
Mold2 is a product designed and produced by the National Center for Toxicological
Research (NCTR).  FDA and NCTR retain ownership of this product.

Please address any questions or suggestions to Dr. Huixiao Hong, National Center for Toxicological
Research, at 870-543-7296 or Huixiao.Hong@fda.hhs.gov.

###################################

Should you publish results based on the Mold² descriptors, please cite:

Mold², Molecular Descriptors from 2D Structures for Chemoinformatics and Toxicoinformatics
Huixiao Hong, Qian Xie, Weigong Ge, Feng Qian, Hong Fang, Leming Shi, Zhenqiang Su, Roger Perkins, and Weida Tong
Journal of Chemical Information and Modeling 2008 48 (7), 1337-1344
DOI: 10.1021/ci800038f

###################################

""")

    def _parse_details(self):
        """Parse descriptor details."""
        with open(os.path.abspath(os.path.join(__file__, os.pardir, 'descriptors.json'))) as jfile:
            self._details = json.load(jfile)

    def _download_executables(self, verbose: bool = True) -> None:
        """Download executables from the FDA website.

        :param verbose: Should details about the download of executables be printed out
        """
        if not os.path.isfile(self._zipfile):
            if verbose:
                # Display information
                print('The executables will be installed and are being downloaded from\n'
                      'https://www.fda.gov/science-research/bioinformatics-tools/mold2')
            # Download Mold2
            session = requests.session()
            res = session.get("https://www.fda.gov/files/science%20&%20research/published/Mold2-Executable-File.zip",
                              headers={"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                                                     "AppleWebKit/537.36 (KHTML, like Gecko) "
                                                     "Chrome/128.0.0.0 "
                                                     "Safari/537.36",
                                       "cache-control": "public, max-age=27440",
                                       "content-type": "text/html; charset=UTF-8",
                                       "accept": "text/html,application/xhtml+xml,application/xml;q=0.9,"
                                                 "image/avif,image/webp,image/apng,*/*;q=0.8,"
                                                 "application/signed-exchange;v=b3;q=0.7",
                                       "accept-encoding": "gzip, deflate, br, zstd",
                                       "accept-language": "en-US;q=0.9,en;q=0.8"},
                              stream=True, verify=True)
            # Save ZIP file
            with open(self._zipfile, 'wb') as fh:
                for chunk in res.iter_content(chunk_size=1024):
                    fh.write(chunk)
            if verbose:
                # Display information
                print('Download and installation are now complete.\n')

    def _extract_executables(self, zip_path: str) -> None:
        """Extract executables from the ZIP file.

        :param zip_path: Path to the zip file containing Mold2 executables
        """
        mold2_folder = os.path.abspath(os.path.join(self._dir, 'extras'))
        # Extract Zip file
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(mold2_folder)
        # Rename files
        os.rename(os.path.join(mold2_folder, 'Mold2', 'Windows', 'Mold2.doc'),
                  os.path.join(mold2_folder, 'Mold2', 'Windows', 'Mold2.exe'))
        os.rename(os.path.join(mold2_folder, 'Mold2', 'Windows', 'Mold2.txt'),
                  os.path.join(mold2_folder, 'Mold2', 'Windows', 'Mold2.bat'))
        # Make files executable
        if platform in ('linux', 'darwin'):
            os.chmod(os.path.join(mold2_folder, 'Mold2', 'Windows', 'Mold2.exe'), stat.S_IXUSR)
            os.chmod(os.path.join(mold2_folder, 'Mold2', 'Linux_x86-32', 'Mold2'), stat.S_IXUSR)
            os.chmod(os.path.join(mold2_folder, 'Mold2', 'Linux_x86-64', 'Mold2'), stat.S_IXUSR)

    def _prepare_input(self, mols: List[Chem.Mol]) -> str:
        """Create temporary V2000 SD input file.

        :param mols: Molecules to obtain Mold2 descriptors of.
        :return: Path to the temporary file
        """
        # Write SD file (input of Mold2)
        sdf = mktempfile('.sdf')
        with open(sdf, 'wt') as outfile, AllChem.SDWriter(outfile) as writer:
            writer.SetForceV3000(False)
            writer.SetProps(('',))
            for mol in mols:
                writer.write(mol)
        fix_bond_blocks(sdf)
        return sdf

    def _prepare_command(self, input_path: str, output_path: str, log: Optional[str] = None) -> str:
        """Obtain command to run to calculate molecular descriptors.

        :param input_path: Path to V2000 SD input file containing molecules to calculate Mold2 descriptors from.
        :param output_path: Path to output file
        :param log: If not None, path to the log file tobe created
        :return: Internal command to be run to obtain descriptors
        """
        mold2_folder = os.path.abspath(os.path.join(self._dir, 'extras'))
        # Exceptions to ensure execution security of subprocess
        if not os.path.isdir(mold2_folder):
            raise RuntimeError('Could not locate Mold2 executables. Were they downloaded first?')
        if not os.path.isfile(input_path):
            raise ValueError(f'Input file does not exist: {input_path}')
        if not os.path.isdir(os.path.abspath(os.path.join(output_path, os.pardir))):
            raise ValueError(f'''Path to output file does not exist: {os.path.abspath(os.path.join(output_path,
                                                                                                   os.pardir))}''')
        # Determine components for the command
        if platform.startswith('win32'):
            exec_path = os.path.join(mold2_folder, 'Mold2', 'Windows', 'Mold2.exe')
            log_file = log if log is not None else 'NUL'
            echo_cmd = 'echo.'
        elif platform.startswith('linux'):
            log_file = log if log is not None else '/dev/null'
            echo_cmd = r'echo -e "\n"'
            if architecture()[0].startswith('32'):
                exec_path = os.path.join(mold2_folder, 'Mold2', 'Linux_x86-32', 'Mold2')
            else:
                exec_path = os.path.join(mold2_folder, 'Mold2', 'Linux_x86-64', 'Mold2')
        else:
            raise RuntimeError(f'Platform ({platform}) not supported.')
        # Create command
        command = f'{echo_cmd} | {exec_path} -i {input_path} -o {output_path} -r {log_file}'
        return command

    def _run_command(self, command: str) -> None:
        """Run the internal command to obtain molecular descriptors.

        :param command: Internal command to be run.
        """
        # Run calculation
        with open(os.devnull, 'wb') as devnull:
            _ = subprocess.check_output(command, shell=True, stderr=devnull)  # noqa: S602

    def _parse_result(self, output_path: str) -> pd.DataFrame:
        """Read a Mold2 output file and convert to a DataFrame.

        :param output_path: Path to a Mold2 output file
        :return: a pandas DataFrame containing all Mold2 desciptor values
        """
        # Read results
        data = pd.read_table(output_path, header=None, usecols=range(1, 778), low_memory=False)
        if data.iloc[0, 0] != 'D001':
            # Header is not provided if first molecule failed
            data.columns = [f'D{x:03d}' for x in range(1, 778)]
        else:
            data.columns = data.iloc[0, :]
            data.drop(index=0, inplace=True)
        data.reset_index(drop=True, inplace=True)
        return data

    def _calculate(self, mols: List[Chem.Mol]) -> pd.DataFrame:
        """Caclulate Mold2 descriptors on one process.

        :param mols: RDkit molecules for which Mold2 descriptors should be calculated
        :return: a pandas DataFrame containing all Mold2 desciptor values
        """
        # Copy of executables for this instance
        if not hasattr(self, '_dir'):
            self._dir = tempfile.mkdtemp(prefix='Mold2_')
        if not os.path.isdir(os.path.join(self._dir, 'extras')):
            # Extract executables
            self._extract_executables(self._zipfile)
        input_path = self._prepare_input(mols)
        output = mktempfile('.txt')
        mold2_command = self._prepare_command(input_path=input_path, output_path=output, log=None)
        self._run_command(mold2_command)
        data = self._parse_result(output_path=output)
        data = data.apply(convert_to_numeric, axis=0, args=(self.fill_na,))
        os.remove(input_path)
        os.remove(output)
        return data


def mktempfile(suffix: str = None) -> str:
    """Return the path to a writeable temporary file."""
    file = tempfile.mkstemp(suffix=suffix)
    os.close(file[0])
    return file[1]


def fix_bond_blocks(path: str) -> None:
    """Fix RDKit bond blocks of molecules in the given SD file.

    Ensures the fields 5 to 7 are included as RDKit does not output them:
    - field 5: not used
    - field 6: bond topology
    - field 7: reacting center status

    :param path: Path to the RDKit SD file
    """
    temp_out = mktempfile()
    bond_line = re.compile('^(?:\s+\d+){4}$')
    with open(temp_out, 'wt') as outfile, open(path, 'rt') as infile:
        lines = [infile.readline() for _ in range(4)]
        outfile.writelines(lines)
        if not lines[-1].strip().endswith('V2000'):
            raise ValueError(f'File provided is not a V2000 SD file: {path}')
        for line in infile:
            if re.match(bond_line, line):
                outfile.write(f'{line.rstrip()}  0  0  0\n')
            else:
                outfile.write(line)
    os.remove(path)
    os.rename(temp_out, path)


def convert_to_numeric(series: pd.Series, fill_value: Optional[Union[float, str]] = None) -> pd.Series:
    """Convert a ``pandas.Series`` to numeric dtype.
    ``pd.to_numeric`` does not handle nans, so this function masks the nans, converts and then
    reinserts them.

    :param series: series to convert.
    :param fill_value: value to replace NaN values by
    :return: the converted series
    """
    # Identify NAs, convert to numeric and set NAs back
    nan_mask = pd.isna(series)
    series[nan_mask] = 0
    column = pd.to_numeric(series, errors='coerce')
    column[nan_mask] = np.NaN
    # Replace NaNs by custom value
    if fill_value is not None:
        mask = pd.isna(column) | (column == 3.4028199999999998e+38)  # undefined Mold2 value
        if np.isreal(fill_value):
            fn = lambda x: fill_value
        elif isinstance(fill_value, str) and fill_value.lower() == 'nan':
            fn     = lambda x: np.NaN
        elif fill_value == 'median':
            fn = np.median
        elif fill_value == 'mean':
            fn = np.mean
        column[mask] = fn(column[~mask])
    return column
