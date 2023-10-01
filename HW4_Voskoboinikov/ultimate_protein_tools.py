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


AA_NAME_DICT: dict[str, str] = {
    'G': 'Gly', 'g': 'Gly',
    'L': 'Leu', 'l': 'Leu',
    'Y': 'Tyr', 'y': 'Tyr',
    'S': 'Ser', 's': 'Ser',
    'E': 'Glu', 'e': 'Glu',
    'Q': 'Gln', 'q': 'Gln',
    'D': 'Asp', 'd': 'Asp',
    'N': 'Asn', 'n': 'Asn',
    'F': 'Phe', 'f': 'Phe',
    'A': 'Ala', 'a': 'Ala',
    'K': 'Lys', 'k': 'Lys',
    'R': 'Arg', 'r': 'Arg',
    'H': 'His', 'h': 'His',
    'C': 'Cys', 'c': 'Cys',
    'V': 'Val', 'v': 'Val',
    'P': 'Pro', 'p': 'Pro',
    'W': 'Trp', 'w': 'Trp',
    'I': 'Ile', 'i': 'Ile',
    'M': 'Met', 'm': 'Met',
    'T': 'Thr', 't': 'Thr'
    }


# TODO check if possible to rempve kwargs
RNA_AA_TABLE = {
'F': ['UUU', 'UUC'],
 'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
 'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
 'Y': ['UAU', 'UAC'],
 '*': ['UAA', 'UAG', 'UGA'],
 'C': ['UGU', 'UGC'],
 'W': ['UGG'],
 'P': ['CCU', 'CCC', 'CCA', 'CCG'],
 'H': ['CAU', 'CAC'],
 'Q': ['CAA', 'CAG'],
 'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
 'I': ['AUU', 'AUC', 'AUA'],
 'M': ['AUG'],
 'T': ['ACU', 'ACC', 'ACA', 'ACG'],
 'N': ['AAU', 'AAC'],
 'K': ['AAA', 'AAG'],
 'V': ['GUU', 'GUC', 'GUA', 'GUG'],
 'A': ['GCU', 'GCC', 'GCA', 'GCG'],
 'D': ['GAU', 'GAC'],
 'E': ['GAA', 'GAG'],
 'G': ['GGU', 'GGC', 'GGA', 'GGG'],
}


RNA_CODON_TABLE = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def read_seq_from_fasta(path_to_seq: str, use_full_name=False, **kwargs):
    with open(path_to_seq) as f:
        out_dct = {}
        for line in f:
            line = line.strip()
            if line.startswith('>'): # check for first line in seq
                if use_full_name: # check if user set full name in fasta
                    name = line[1:] # take whole fasta properties (e.g. if names not unique)
                else:
                    name = line[1:].split()[0]
            else:
                out_dct[name] = out_dct.get(name, '') + line # get value from dict (return '' if empty) and append str
    return out_dct


def get_sites_lengths(sites):
    sites_length_dct = {}
    for site in sites:
        sites_length_dct[site] = len(site)
    return sites_length_dct


def invert_dct(dct):
    inv_dct = {}
    for k, v in dct.items():
        inv_dct[v] = inv_dct.get(v, []) + [k] # get value from dict (return []) and append key
    return inv_dct


def find_sites(seq, *sites, is_one_based = False, **kwargs):
    window_sizes = invert_dct(get_sites_lengths(sites)) # get lengths of all sites and stick them together to avoid passing through seq multiple times if possible
    found_sites = {} 
    for window_size in window_sizes: # perform iteration for all given lengths of sites
        for i in range(len(seq) - window_size + 1): # iterate through seq with step one and consider window of site length each iteration 
            scatter = seq[i:i + window_size] # get fragment of sequence with length of window i.e. scatter
            for site in window_sizes[window_size]:
                if scatter == site: # check if scatter is site
                    found_sites[site] = (
                        found_sites.get(site, []) # get 
                        + [i + is_one_based]
                        ) # append index to list in dict
    return found_sites


def get_protein_rnas(seq, i_absolutely_fucking_know_what_im_doing = False):
    if i_absolutely_fucking_know_what_im_doing:
        kmers = [''] # set initial kmers
        for amino_acid in seq: # iterate AAs
            current_kmers = []
            codons = RNA_AA_TABLE[amino_acid] # get list of codons for AA
            for codon in codons:
                for kmer in kmers:
                    current_kmers.append(kmer + codon) # append every codon to existing kmers
            kmers = current_kmers # re-write k-mers for next iteration

        return kmers

    return "You don't fucking know what you're doing!" # politely ask user to reconsider their actions


def get_protein_rnas_number(seq):
    rnas_num = 1
    for amino_acid in seq:
        rnas_num *= len(RNA_AA_TABLE[amino_acid])
    return rnas_num


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


def get_atomic_mass(chem: str, atomic_mass: dict[str, float] = None) -> float:
    """

    Calculates the molecular mass of a biological molecule, primarily an amino acid, based on a simple chemical formula.

    Arguments / Args:
    - chem (str): String representing a simple chemical formula, e.g. C2H5OH
    - atomic_mass (dict[str, float], optional): A dictionary linking the chemical elements Carbon, Hydrogen, Oxygen,
    Nitrogen, and Sulfur with their masses in atomic mass units.

    Return:
    - float: Molecular mass of a biological molecule in atomic mass units.
    """

    total_mass = 0
    char = 0  # idx init
    if atomic_mass is None:
        atomic_mass = ATOMIC_MASS
    while char < len(chem):
        if chem[char].isalpha():
            element = chem[char]
            char += 1  # очень надо, а то я опять бесконечный цикл сделала
            if char < len(chem) and chem[char].isdigit():
                number = ''
                while char < len(chem) and chem[char].isdigit():
                    number += chem[char]
                    char += 1  # очень надо
                total_mass += atomic_mass[element] * int(number)
            else:
                total_mass += atomic_mass[element]
        else:
            raise ValueError(f'Unknown elem: {chem[char]}')

    return total_mass


def convert_aa_name(sequence: str, name_dict: dict[str, str] = None, sep: str = '',
                    use_default_register: bool = True) -> str:
    """

    Converts a sequence of one-letter amino acid codes to three-letter designations.

    Arguments / Args:
    - sequence (str): String with one-letter amino acid codes.
    - name_dict (dict[str, str], optional): A dictionary linking one-letter codes to three-letter designations.
    If not provided, the standard AA_NAME_DICT dictionary is used.
    - sep (str, optional): Separator between three-letter amino acid designations. There is no delimiter by default.
    - use_default_register(bool, optional): Determines whether to preserve letter case in three-letter designations.
    If True, the letters will be converted to upper or lower case depending on the case of the depending
    on the case of the one-letter code. The default is False.

    Return:
    - str: A string of three-letter amino acid designations separated by the specified delimiter.
    """
    
    new_name = ''
    if name_dict is None:
        name_dict = AA_NAME_DICT
    for i, aa in enumerate(sequence):
        if aa in name_dict:
            if use_default_register is False:
                new_name += name_dict[aa]
            elif use_default_register is True:
                if aa.isupper():
                    new_name += name_dict[aa].upper()
                else:
                    new_name += name_dict[aa].lower()
            else:
                if aa.isupper():
                    new_name += name_dict[aa].lower()
                else:
                    new_name += name_dict[aa].upper()
            if sep and (i + 1) < len(sequence):
                new_name += sep
        else:
            raise ValueError(f'Unknown amino acid: {aa}')
    return new_name
