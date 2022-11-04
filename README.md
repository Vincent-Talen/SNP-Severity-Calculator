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
This directory contains two protein family's protein sequences and one of the genes from the families as their coding DNA sequence.
These can be used to test the functionality of the script and are used as defaults when no family and sequence is given by the user.

## Installation
This script was created in, and for, Python 3.10, backwards compatibility is not guaranteed!
Python can be installed from the following link: [Python.org Downloads](https://www.python.org/downloads/)

### Required Packages
The following Python packages are required for this project and should be installed using pip in a terminal with `pip install package_name`:
- biopython

## Usage
To use the SNP Severity Check, the `main.py` script has to be run from the terminal.  
For further help use the `-h` argument with the script like `python3 main.py -h`.  

Examples:  

## Useful links
* [Biopython webpage](https://biopython.org/)

## Contact
Vincent Talen  
v.k.talen@st.hanze.nl
