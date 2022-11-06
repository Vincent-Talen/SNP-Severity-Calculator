# SNP Severity Calculator
**Hanzehogeschool Groningen: Bioinformatics Year 3, Period 9**

A Python script that calculates the severity of a given Single Nucleotide Polymorphism (SNP) in a given genomic sequence protein based on the conservation in its protein family.

## About the assignment
This assignment is part of the course Bio-Informatics 3 (BFVH19BIN3), constituting to 20 percent of the final grade.  
The goal was to create a script that when given a coding genomic sequence, the position and nucleotide for an SNP can determine if this mutation is deleterious based on the conservation of the coding genomic sequence's protein family.

## Repository File Structure
### Project Tree
```bash
SNP-Severity-Calculator
├── main.py
├── classes.py
├── LICENSE
├── README.md
└── testdata
    └── *
```

### / testdata
This directory contains protein families' protein sequences and a gene from each family as coding DNA sequence.
These protein families can be used to test the functionality of the script and one is used as defaults when no family and sequence is given by the user.

## Installation
This script was created in, and for, Python 3.10, backwards compatibility is not guaranteed!
Python can be installed from the following link: [Python.org Downloads](https://www.python.org/downloads/)

### Clustal Omega for Multiple Sequence Alignment
If it is desired to input multi-fasta files with just the proteins of a family and 
let the script create multiple sequence alignments for the family, then it is required to install `Clustal Omega`.
Clustal Omega can be downloaded from http://www.clustal.org/omega/. 

> For Windows users please follow the instructions found there and make sure to give the `-c, --clustalo_path` argument the path to the executable file.

> Mac users should download the 'Standalone Mac Binary' file and place it into the root directory of this repository. 
Rename the file to `clustalo`, open a terminal in the same directory, type or paste `chmod u+x clustalo` in the terminal and then press enter.
If typing `./clustalo` in the terminal does not work, the following step needs to be taken:   
*In Finder right-click clustalo, which is now an executable, navigate to 'Open With' and click on Terminal. 
The computer will now show a prompt that the author of the file is unrecognized, click on 'Open' to override this. 
After this clustalo can be used normally.*

### Required Packages
The following Python packages are required for this project and should be installed using pip in a terminal with `pip install package_name`:
- biopython

## Usage
To use the SNP Severity Check, the `main.py` script has to be run from the terminal.  
For further help use the `-h` argument with the script like `python3 main.py -h`.  

Examples:  
> $python3 main.py 42 G

> $python3 main.py 42 G -s testdata/ZCCHC17_C-Lupus_DNA.fasta -f testdata/ZCCHC17_Protein_Family.fasta

> $python3 main.py 42 G --sequence testdata/ZCCHC17_C-Lupus_DNA.fasta --family testdata/ZCCHC17_Protein_Family.fasta

## Useful links
* [HomoloGene: ZCCHC17 Protein Family](https://www.ncbi.nlm.nih.gov/homologene/32319)
* [HomoloGene: HIST1H2BA Protein Family](https://www.ncbi.nlm.nih.gov/homologene/69356)
* [Biopython webpage](https://biopython.org/)

## Contact
Vincent Talen  
v.k.talen@st.hanze.nl
