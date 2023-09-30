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