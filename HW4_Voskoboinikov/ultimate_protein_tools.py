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


def length_of_protein(seq: str) -> int:
    """
    Calculates the length of a protein.

    Argument:
    - seq (str): sequence to calculate the length

    Return:
    - int: sequence length
    """

    return len(seq)