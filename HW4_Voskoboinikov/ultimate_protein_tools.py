AMINOACID_DICT = {
    'A': 'Alanine', 'a': 'alanine',
    'C': 'Cysteine', 'c': 'cysteine',
    'D': 'Aspartic acid', 'd': 'aspartic acid',
    'E': 'Glutamic acid', 'e': 'glutamic acid',
    'F': 'Phenylalanine', 'f': 'Phenylalanine',
    'G': 'Glycine', 'g': 'glycine',
    'H': 'Histidine', 'h': 'histidine',
    'I': 'Isoleucine', 'i': 'isoleucine',
    'K': 'Lysine', 'k': 'lysine',
    'L': 'Leucine', 'l': 'leucine',
    'M': 'Methionine', 'm': 'methionine',
    'N': 'Asparagine', 'n': 'asparagine',
    'P': 'Proline', 'p': 'proline',
    'Q': 'Glutamine', 'q': 'glutamine',
    'R': 'Arginine', 'r': 'arginine',
    'S': 'Serine',  's': 'serine',
    'T': 'Threonine', 't': 'threonine',
    'V': 'Valine', 'v': 'valine',
    'W': 'Tryptophan', 'w': 'tryptophan',
    'Y': 'Tyrosine', 'y': 'tyrosine'
    }


H2O_WEIGHT: float = 18.01468


AA_MASS_DICT: dict[str, float] = {
    'G': 75.0659, 'g': 75.0659,
    'L': 131.17262, 'l': 131.17262,
    'Y': 181.18894, 'y': 181.18894,
    'S': 105.09158, 's': 105.09158,
    'E': 147.12826, 'e': 147.12826,
    'Q': 146.1438, 'q': 146.1438,
    'D': 133.10158, 'd': 133.10158,
    'N': 132.11712, 'n': 132.11712,
    'F': 165.18994, 'f': 165.18994,
    'A': 89.09258, 'a': 89.09258,
    'K': 146.18716, 'k': 146.18716,
    'R': 174.20056, 'r': 174.20056,
    'H': 155.15466, 'h': 155.15466,
    'C': 121.15758, 'c': 121.15758,
    'V': 117.14594, 'v': 117.14594,
    'P': 115.13026, 'p': 115.13026,
    'W': 204.22648, 'w': 204.22648,
    'I': 131.17262, 'i': 131.17262,
    'M': 149.21094, 'm': 149.21094,
    'T': 119.11826, 't': 119.11826,
    }


ATOMIC_MASS: dict[str, float] = {
    'C': 12.011,
    'H': 1.00784,
    'O': 15.999,
    'N': 14.0067,
    'S': 32.065
}


def length_of_protein(seq: str) -> int:
    """
    Calculates the length of a protein.

    Argument:
    - seq (str): sequence to calculate the length

    Return:
    - int: sequence length
    """

    return len(seq)


def count_aa(seq: str, *, aminoacids: str = None) -> dict:
    """
    Counts the number of given or all amino acids in a protein sequence.

    Arguments:
    - seq (str): sequence to count amino acids
    - aminoacids (str): which amino acids to count in sequence

    Return:
    - dict: a dictionary with amino acids and its count
    """

    aa_dict_count = {}
    if (aminoacids is None) or (aminoacids == ''):
        '''
        I added an additional condition for user-friendly experience.
        E.g., we can want to find specific aminoacid, look on result and then look on all aminoacids.
        Without this condition we have to delete keyword argument, but with it we can only make it empty.
        '''
        aminoacids = ''.join(set(seq))
    for aa in aminoacids:
        aa_dict_count[aa] = seq.count(aa)
    return aa_dict_count


def get_fracture_of_aa(seq: str, *, show_as_percentage: bool = False, aminoacids: str = None) -> dict:
    """
    Returns the fracture or percentage of amino acids in a protein sequence.

    Arguments:
    - seq (str): sequence in which you need to calculate the fracture of amino acids
    - show_as_percentage (bool): change it to True, if you want to get results with percentages
    - aminoacids (str): the fracture of which amino acids to count in the sequence

    Return:
    - dict: a dictionary with amino acids and its fracture or percentage
    """

    if show_as_percentage:
        mult = 100
        round_var = 2
    else:
        mult = 1
        round_var = 4
    aa_dict_count = count_aa(seq, aminoacids=aminoacids)
    aa_dict_percent = {}
    len_of_protein = length_of_protein(seq)
    for aa, count in aa_dict_count.items():
        aa_dict_percent[aa] = round(count / len_of_protein * mult, round_var)
    return aa_dict_percent


def calculate_protein_mass(sequence: str, aa_atomic_mass: dict[str, float] = None) -> float:
    """

    Calculates the molecular mass of a protein based on its amino acid sequence and a dictionary of amino acid masses.

    Arguments / Args:
    - sequence(str or list): A string or list of characters representing the amino acid sequence.
    - aa_atomic_mass(dict): A dictionary linking amino acids to their masses in atomic mass units.
    
    Return:
    - float: The molecular mass of a protein in atomic mass units, rounded to the third decimal place.
    """

    total_mass = 0.0
    if aa_atomic_mass is None:
        aa_atomic_mass = AA_MASS_DICT

    for aa in sequence:
        if aa in aa_atomic_mass:
            total_mass += aa_atomic_mass[aa]
        else:
            raise ValueError(f'Unknown amino acid: {aa}')
    total_mass = total_mass - H2O_WEIGHT * (len(sequence) - 1)

    return round(total_mass, 3)
