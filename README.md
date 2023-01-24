[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Mold² Python wrapper

Python wrapper to ease the calculation of Mold2 molecular descriptors.

## Copyright notice

Mold2 is a product designed and produced by the National Center for Toxicological
Research (NCTR).<br/>FDA and NCTR retain ownership of this product.

Olivier J. M. Béquignon is **neither** the copyright holder of Mold2 **nor** responsible for it.

Only the Python wrapper is the work of Olivier J. M. Béquignon.

## Installation

From source:

    git clone https://github.com/OlivierBeq/Mold2_pywrapper.git
    pip install ./Mold2_pywrapper

with pip:

```bash
pip install mold2-pywrapper
```

### Get started

```python
from Mold2_pywrapper import Mold2

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
mols = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

mold2 = Mold2()
print(mold2.calculate(mols))
```

Instantiating a Mold2 object ensures that the executables for your platform are accessible.
If this is not the case, an attempt to download them from the
[website of the FDA](https://www.fda.gov/science-research/bioinformatics-tools/mold2) is made.

Should one have downloaded the original ZIP file available from the
[website of the FDA](https://www.fda.gov/science-research/bioinformatics-tools/mold2), the executables can be installed
using the following:

```python
path_to_zipfile = '...'  # Replace by the path to the ZIP file on your machine
mold2 = Mold2.from_executable(path_to_zipfile)
print(mold2.calculate(mols))
```

Executables will be installed for future use. From then on, default instanciation may be carried out:

```python
mold2 = Mold2()
print(mold2.calculate(mols))
```

### Details about descriptors

Any detail about the 777 Mold2 descriptors can be obtained either for a single descriptor, by providing its index:

```python
print(mold2.descriptor_detail(15))
# rotatable bond fraction
```

Or for all at once:

```python
print(mold2.descriptor_details())
# {"D001": "number of 6-membered aromatic rings (only carbon atoms)",
#  "D002": "number of 03-membered rings",
#  ...
#  }
```

## Documentation

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


