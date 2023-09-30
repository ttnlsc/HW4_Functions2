# Ultimate Protein Tools

### Overview

This project contains a `ultimate_protein_tools.py` program, which implements the `run_ultimate_protein_tools()` function, the input of which is protein sequences and the action that needs to be applied to them. 

### Installation

TODO import from module!!!!

### Usage

To run the script, just call it from the directory where the program is located:

```
python ultimate_protein_tools.py
```

To exit the program, type `exit` in the console.

While running the program, you can use next arguments for ***run_ultimate_protein_tools*** function (next – *main function*):

- `get_length_of_protein`: Calculates the length of a protein.
- `count_aa`: Counts the number of amino acids in a protein sequence.
- `get_fracture_of_aa`: Returns the fracture or percentage of amino acids in a protein sequence.

- `read_seq_from_fasta`: 
- `find_sites`: 
- `get_protein_rnas_number`: 
- `get_protein_rnas`: 
- `get_frameshift_proteins`: 

- `calculate_protein_mass`: Calculates the molecular mass of a protein based on its amino acid sequence and a dictionary of amino acid masses.
- `get_atomic_mass`: Calculates the molecular mass of a biological molecule, primarily an amino acid, based on a simple chemical formula.
- `convert_aa_name`: Converts a sequence of one-letter amino acid codes to three-letter designations.

For additional information about arguments and output, please read the docstring for the desired function.

### Input of data

During each run of the main function, the user is required to enter a **protein sequence / sequences** that must be processed using the procedures listed above.

The program involves the analysis of protein sequences consisting of <u>**20 canonical amino acids**</u>.

If the data is entered incorrectly, an appropriate error will be displayed.

```python
run_ultimate_protein_tools('AZAZA', get_length_of_protein) -> ValueError #TODO add error message
```

### Examples

```python
run_ultimate_protein_tools('MAGDVLAGTTTSDRAAGALGTLGTAATLRAATDGLLQR', get_length_of_protein) -> 38
run_ultimate_protein_tools('MAGDVLAGTTTSDRAAGALGTLGTAATLRAATDGLLQR', aminoacids='AT', count_aa) -> {'A': 9, 'T': 7}
run_ultimate_protein_tools('MAGDVLAGTTTSDRAAGALGTLGTAATLRAATDGLLQR', aminoacids='L', get_fracture_of_aa) -> {'L': 0.1579}
run_ultimate_protein_tools('MAGDVLAGTTTSDRAAGALGTLGTAATLRAATDGLLQR', aminoacids='DRG', get_fracture_of_aa, show_as_percentage=True) -> {'D': 7.89, 'R': 7.89, 'G': 15.79}

#TODO examples for other functions
```

### Troubleshooting

TODO change this

If the program doesn't work – try to scream like opossum.

### Contacts
![Wonderful Team](https://github.com/ArtemVaska/HW4_Functions2/tree/HW4_Voskoboinikov/HW4_Voskoboinikov/Wonderful_Team.jpg)

Aleksandr Voskoboinikov – Team Leader (wwoskie@gmail.com)

Artem Vasilev (artem_vasilev_01@list.ru)

Tatiana Lisitsa (ttnlsc@gmail.com)