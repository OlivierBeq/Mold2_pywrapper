# Mold² Python Wrapper

<!-- Badges -->
<div align="center">

[![PyPI version](https://img.shields.io/pypi/v/mold2-pywrapper.svg)](https://pypi.org/project/mold2-pywrapper/)
[![Supported Python versions](https://img.shields.io/pypi/pyversions/mold2-pywrapper.svg)](https://pypi.org/project/mold2-pywrapper/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![CI](https://github.com/OlivierBeq/Mold2_pywrapper/actions/workflows/ci.yml/badge.svg)](https://github.com/OlivierBeq/Mold2_pywrapper/actions/workflows/ci.yml)
[![Code Quality](https://github.com/OlivierBeq/Mold2_pywrapper/actions/workflows/linting.yml/badge.svg)](https://github.com/OlivierBeq/Mold2_pywrapper/actions/workflows/linting.yml)
[![codecov](https://codecov.io/gh/OlivierBeq/Mold2_pywrapper/graph/badge.svg?token=Q49XHK5FLB)](https://codecov.io/gh/OlivierBeq/Mold2_pywrapper)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![DOI](https://zenodo.org/badge/592541962.svg)](https://doi.org/10.5281/zenodo.16936866)
<br>
![Windows](https://img.shields.io/badge/Windows-0078D6?style=for-the-badge&logo=windows&logoColor=white)
![Linux](https://img.shields.io/badge/Linux-FCC624?style=for-the-badge&logo=linux&logoColor=black)

</div>


A simple and reliable Python wrapper for calculating Mold² molecular descriptors. This library automates the process of downloading, installing, and running the Mold² executables for seamless integration into your cheminformatics workflows.

## ✨ Features

-   **Automated Installation**: Automatically downloads and sets up the correct Mold² executables for your operating system (Windows or Linux).
-   **Simple API**: A clean and easy-to-use interface for calculating descriptors from RDKit molecule objects.
-   **Parallel Processing**: Calculate descriptors for large datasets quickly using built-in support for multiprocessing.
-   **Data Handling**: Returns descriptors in a convenient Pandas DataFrame, ready for analysis.
-   **Descriptor Information**: Easily access details and descriptions for all 777 Mold² descriptors.

## ✍️ Copyright and Citation Notice

#### Mold² Software
Mold² is a product designed and produced by the National Center for Toxicological Research (NCTR). FDA and NCTR retain ownership of this product. This Python wrapper is the work of Olivier J. M. Béquignon and is **not** affiliated with the original authors.

#### Citing
If you use this wrapper in your research, please cite the original Mold² publication in addition to this software package:

1.  **Original Mold² Paper:**
    > Hong, H., Xie, Q., Ge, W., Qian, F., Fang, H., Shi, L., Su, Z., Perkins, R., & Tong, W. (2008). Mold², Molecular Descriptors from 2D Structures for Chemoinformatics and Toxicoinformatics. *Journal of Chemical Information and Modeling*, 48(7), 1337–1344.
    > [DOI: 10.1021/ci800038f](https://doi.org/10.1021/ci800038f)

2.  **This Wrapper:**
    > Please refer to the `CITATION.cff` file for details on how to cite this library.


## 🛠️ Requirements

-   Python 3.9+
-   [RDKit](https://www.rdkit.org/docs/Install.html)

## 📦 Installation

You can install the wrapper from PyPI:

```bash
pip install mold2-pywrapper
```

## 💡  Usage

Here is a basic example of how to calculate Mold² descriptors for a list of molecules.

```python
from rdkit import Chem
from Mold2_pywrapper import Mold2

# A list of SMILES strings
smiles_list = [
    # erlotinib
    "n1cnc(c2cc(c(cc12)OCCOC)OCCOC)Nc1cc(ccc1)C#C",
    # midecamycin
    "CCC(=O)O[C@@H]1CC(=O)O[C@@H](C/C=C/C=C/[C@@H]([C@@H](C[C@@H]([C@@H]([C@H]1OC)O[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)C)O[C@H]3C[C@@]([C@H]([C@@H](O3)C)OC(=O)CC)(C)O)N(C)C)O)CC=O)C)O)C",
    # selenofolate
    "C1=CC(=CC=C1C(=O)NC(CCC(=O)OCC[Se]C#N)C(=O)O)NCC2=CN=C3C(=N2)C(=O)NC(=N3)N",
    # cisplatin
    "N.N.Cl[Pt]Cl"
]

# Convert SMILES to RDKit molecule objects
mols = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

# Instantiate the Mold2 calculator
# This will automatically download the executables on first run.
mold2 = Mold2()
# Calculate descriptors
descriptors_df = mold2.calculate(mols)

print(descriptors_df.head())
```


### ⚙️ Working with Pre-Downloaded Executables

If you have already downloaded the Mold² ZIP file from the [FDA website](https://www.fda.gov/science-research/bioinformatics-tools/mold2), you can install the executables directly from the file:

```python
path_to_zipfile = 'path/to/your/Mold2-Executable-File.zip'
mold2 = Mold2.from_executable(path_to_zipfile)

# The executables are now installed for future use.
# You can now instantiate Mold2() normally.
```

### 🔍 Getting descriptors details

You can retrieve the description for any of the 777 descriptors by its index or get a dictionary of all of them.

```python
# Get details for a single descriptor
print(mold2.descriptor_detail(15))
# Output: rotatable bond fraction

# Get details for all descriptors
all_descriptors = mold2.descriptor_details()
print(all_descriptors['D001'])
# Output: number of 6-membered aromatic rings (only carbon atoms)
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/OlivierBeq/Mold2_pywrapper/blob/master/LICENSE) file for details.

## 📚 API Documentation

```python
def calculate(mols, show_banner=True, njobs=1, chunksize=100):
```

Default method to calculate #### Parameters

- ***mols  : Iterable[Chem.Mol]***  
  RDKit molecule objects for which to obtain Mold2 descriptors.
- ***show_banner  : bool***  
  Displays default notice about Mold2 descriptors.
- ***njobs  : int***  
  Maximum number of simultaneous processes.
- ***chunksize  : int***  
  Maximum number of molecules each process is charged of.
- ***return_type  : pd.DataFramce***  
  Pandas DataFrame containing Mold2 molecular descriptors.Mold2 descriptors.
  If executables have not previously been downloaded, attempts to download and install them.

________________

```python
def descriptor_detail(index):
```

Obtain detils about one descriptor.

#### Parameters

- ***index  : int***  
  Index of the descriptor.
- ***return_type  : str***  
  Description of the descriptor.

________________

```python
def descriptor_details():
```

Obtain details about all descriptors.

#### Parameters

- ***return_type  : dict***  
  Mapping of molecular descriptors with their details.

________________

```python
def from_executable(zipfile_path):
```

Install executables and instantiate a Mold2 calculator.

#### Parameters

- ***zipfile_path  : str***  
  Path to the user-downloaded ZIP file containing Mold2 executables.
