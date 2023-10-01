# Ultimate Protein Tools

### Overview

This project contains a `ultimate_protein_tools.py` program, which implements the `run_ultimate_protein_tools()` function, the input of which is protein sequences and the action that needs to be applied to them. 

### Installation

```python
import ultimate_protein_tools as upt
```

Make sure the path to the directory with `ultimate_protein_tools.py` is added to the PATH so that Python can find it when importing.

### Usage

To run the script, just call it from the directory where the program is located:

```
python ultimate_protein_tools.py
```

To exit the program, type `exit` in the console.

While running the program, you can use next arguments for ***run_ultimate_protein_tools*** function (next – *main function*):

1. `read_seq_from_fasta`: Reads sequences from fasta file and returns dictionary.
    
    Arguments:
    - path_to_seq (str): path to file

    Return:
    - dict: dict of sequences names as keys and sequences themselves as values {'seq_name': 'sequence',}

2. `find_sites`: Finds indexes of given sites.

    Arguments:
    - seq (str): seq to be checked
    - *args (str): sites to be found
    - is_one_based (bool): whether result should be 0- (False) or 1-indexed (True). Default False

    Return:
    - dict: dictionary of sites as keys and lists of indexes for the site where it's been found

3. `get_protein_rnas`: Returns list of all possible RNA's from which can serve as matrix for protein synthesis.

    WARNING: can be computationally intensive on longer sequences, will NOT start unless check_if_user_conscious is True!

    Arguments:
    - seq (str): seq to be checked
    - check_if_user_conscious (bool): checks user's consciousness. Default False

    Return:
    - list: list of possible RNA's as str

4. `get_protein_rnas_number`: Get number of all possible RNA's for a given protein.

    Arguments:
    - seq (str): seq to be checked

    Return:
    - int: number of possible RNA's for seq

5. `get_frameshift_proteins`: Returns list of all possible proteins from all possible frames in peptide.

    WARNING: can be computationally intensive on longer sequences, will NOT start unless check_if_user_conscious is True!
    
    Arguments:
    - seq (str): seq to be checked
    - check_if_user_conscious (bool): checks user's consciousness. Default False
    - is_stop_codon_termination_enabled (bool): terminate translation when reached stop-codon. Default False.

    Return:
    - dict: dict of lists of all possible frames proteins:
    {frame_0: ['protein_seqs'], frame_1: ['protein_seqs'], frame_2: ['protein_seqs']}

6. `get_length_of_protein`: Calculates the length of a protein.

   Argument:
   - seq (str): sequence to calculate the length

   Return:
   - int: sequence length

7. `count_aa`: Counts the number of given or all amino acids in a protein sequence.

   Arguments:
   - seq (str): sequence to count amino acids
   - aminoacids (str): which amino acids to count in sequence. If you want to count all amino acids in the whole sequence, you can provide empty string to this argument or just don't provide this keyword

   Return:
   - dict: a dictionary with amino acids and its count

8. `get_fracture_of_aa`: Calculates the fracture or percentage of amino acids in a protein sequence.

    Arguments:
    - seq (str): sequence in which you need to calculate the fracture of amino acids
    - show_as_percentage (bool): change it to True, if you want to get results with percentages
    - aminoacids (str): the fracture of which amino acids to count in the sequence

    Return:
    - dict: a dictionary with amino acids and its fracture or percentage

9. `calculate_protein_mass`: Calculates the molecular mass of a protein based on its amino acid sequence and a dictionary of amino acid masses.

    Arguments:
    - sequence(str or list): A string or list of characters representing the amino acid sequence.
    - aa_atomic_mass(dict): A dictionary linking amino acids to their masses in atomic mass units.
    
    Return:
    - float: The molecular mass of a protein in atomic mass units, rounded to the third decimal place.

10. `get_atomic_mass`: Calculates the molecular mass of a biological molecule, primarily an amino acid, based on a simple chemical formula.

    Arguments:
    - chem (str): String representing a simple chemical formula, e.g. C2H5OH
    - atomic_mass (dict[str, float], optional): A dictionary linking the chemical elements Carbon, Hydrogen, Oxygen,
    Nitrogen, and Sulfur with their masses in atomic mass units.

    Return:
    - float: Molecular mass of a biological molecule in atomic mass units.

11. `convert_aa_name`: Converts a sequence of one-letter amino acid codes to three-letter designations.

    Arguments:
    - sequence (str): String with one-letter amino acid codes.
    - name_dict (dict[str, str], optional): A dictionary linking one-letter codes to three-letter designations.
    If not provided, the standard AA_NAME_DICT dictionary is used.
    - sep (str, optional): Separator between three-letter amino acid designations. There is no delimiter by default.
    - use_default_register(bool, optional): Determines whether to preserve letter case in three-letter designations.
    If True, the letters will be converted to upper or lower case depending on the case of the depending
    on the case of the one-letter code. The default is False.

    Return:
    - str: A string of three-letter amino acid designations separated by the specified delimiter.

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
![Wonderful Team](https://github.com/ArtemVaska/HW4_Functions2/blob/HW4_Vasilev/HW4_Voskoboinikov/Wonderful_team.jpg)

Aleksandr Voskoboinikov – Team Leader (wwoskie@gmail.com)

Artem Vasilev (artem_vasilev_01@list.ru)

Tatiana Lisitsa (ttnlsc@gmail.com)
