
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from pandas import DataFrame


def read_clustal(filename, data_dict):
    """
    Reads a clustal protein alignment file. Raises an exception if more than
    two sequences are in the alignment.

    Parameters
    ----------
    filepath: str
        path to the clustal file

    data_dict: dict
        The data dictionary to be populated. Defined above.

    TODO
    ----
    Add a try/except to check that the path and file exists. Import os stuff.
    """

    for i, s in enumerate(SeqIO.parse(filename, format='clustal')):
        try:
            assert i < 2
        except:
            raise Exception("More than two sequences.")
        data['protein{}'.format(i+1)]['id'].append(s.id)
        data['protein{}'.format(i+1)]['seq_string'] = str(s.seq)
        data['protein{}'.format(i+1)]['description'].append(s.description)
        data['protein{}'.format(i+1)]['name'].append(s.name)
        data['protein{}'.format(i+1)]['seq_list'] = list(str(s.seq))



def read_dna(filename, data):
    """
    Reads a fasta file that contains a single DNA sequence. Raises an exception if
    more than one sequence is in the file.

    Parameters
    ----------
    filepath: str
        path to the fasta file

    data_dict: dict
        The data dictionary to be populated. Defined above.

    TODO
    ----
    Add a try/except to check that the path and file exists. Import os stuff.

    """


    for i, s in enumerate(SeqIO.parse(filename, format='fasta')):
        try:
            assert i < 1
        except:
            raise Exception("More than one sequence in file.")
        data['dna{}'.format(i+1)]['id'].append(s.id)
        data['dna{}'.format(i+1)]['seq_string'] = str(s.seq)
        data['dna{}'.format(i+1)]['description'].append(s.description)
        data['dna{}'.format(i+1)]['name'].append(s.name)
        s = data['dna1']['seq_string']
        codon_list = [s[i:i+3].upper() for i in range(0, len(s), 3)]
        data['dna{}'.format(i+1)]['codon_list'] = codon_list
        t = CodonTable.standard_dna_table
        data['dna{}'.format(i+1)]['translated_list'] = [t.forward_table[i] for i in codon_list]



data = {
    'protein1': {
                 'id': [],
                 'seq_string':[],
                 'seq_list':[],
                 'reduced_list':[],
                 'description': [],
                 'name':[]
                },
    'protein2': {
                 'id': [],
                 'seq_string':[],
                 'seq_list':[],
                 'reduced_list':[],
                 'description': [],
                 'name':[]
                },
    'dna1': {
             'id': [],
             'seq_string':[],
             'description': [],
             'name':[],
             'codon_list':[],
             'translated_list':[]
            },
    'diff_dict': {
                  'swap_index': [],
                  'swap_old_aa': [],
                  'swap_new_aa': [],
                  'swap_new_codon': [],
                  'indel_index': [],
                  'indel_new_aa': [],
                  'indel_new_codon': []
                 },
    'meta': {'time':'',
            }
}



def find_indels(data):
    """
    Finds indices and replacement amino acids for indels in the
    template protein sequence.

    Parameters
    ----------
    data: dict
        A dictionary defined previously that contains protein
        sequence lists.

    Returns
    -------
    Nothing. Populates the data dictionary with reduced protein
    sequences, the indel list, and insertion list.
    """

    # Copying new instances of the protein lists
    temp_prot1_list = data['protein1']['seq_list'][:]
    temp_prot2_list = data['protein2']['seq_list'][:]

    # Getting a list of indels in the first sequence
    indel_list = [i for i,res in enumerate(temp_prot1_list) if res=='-']

    # Getting the inserted residue at the indel postions
    insert_list = [data['protein2']['seq_list'][i] for i in indel_list]

    # Remove indels, in reverse order to maintain indices,
    # from temp_seq_list
    for index in sorted(indel_list, reverse=True):
        del temp_prot1_list[index]
        del temp_prot2_list[index]

    # Set the reduced protein sequence lists, indel list, and insert list
    # as elements in the data dictionary.
    # NOT SURE WHAT THE STRUCTURE OF THIS DICTIONARY SHOULD BE
    data['protein1']['reduced_list'] = temp_prot1_list
    data['protein2']['reduced_list'] = temp_prot2_list
    data['diff_dict']['indel_index'] = indel_list
    data['diff_dict']['indel_new_aa'] = insert_list



def find_differences(data):
    """
    Find the differences between two clustal aligned protein sequences.

    Parameters
    ----------
    data: dict
        A nested dictionary defined by the isor_primer program.

    Returns
    -------
    Nothing. This populates the diff_dict in that larger supplied data dictionary.

    TODO
    ----
    MAKE SURE THIS FUNCTION CALL COMES AFTER THE FIND_INDELS CALL.
    """

    # I SHOULD PUT THIS AS A GENERAL CHECK AT THE BEGINNING OF MY SCRIPT
    """try:
        assert len(prot1) > 0 and len(prot2) > 0
    except:
        raise Exception("Make sure you have two protein sequences loaded.")"""

    # USE THE LIST VERSION OF THE PROTEIN SEQUENCE
    # LIST INDEXING IS FASTER
    prot1 = data['protein1']['reduced_list']
    prot2 = data['protein2']['reduced_list']

    # Comparing the two sequences pairwise
    for i,res in enumerate(prot1):
        if prot1[i] != prot2[i]:
            data['diff_dict']['swap_index'].append(i)
            data['diff_dict']['swap_old_aa'].append(prot1[i])
            data['diff_dict']['swap_new_aa'].append(prot2[i])



def codon_swap(data):
    """
    Grabs new codons given a one-letter amino acid code.

    Parameters
    ----------
    data: dict
        A nested dictionary defined by the isor_primer program.

    Returns
    -------
    Nothing. This adds to the supplied data dictionary.

    TODO
    ----
    MAKE SURE THIS FUNCTION CALL COMES AFTER THE FIND_INDELS CALL.
    """

    # The back translation table for NNS codon with human preference
    back_table = {
        'A':'GCC',
        'C':'TGC',
        'D':'GAC',
        'E':'GAG',
        'F':'TTC',
        'G':'GGC',
        'H':'CAC',
        'I':'ATC',
        'K':'AAG',
        'L':'CTG',
        'M':'ATG',
        'N':'AAC',
        'P':'CCC',
        'Q':'CAG',
        'R':'CGG',
        'S':'AGC',
        'T':'ACC',
        'V':'GTG',
        'W':'TGG',
        'Y':'TAC',
        '-':''
       }

    # New codons for indel insertion
    data['diff_dict']['indel_new_codon'] = [back_table[res] for res in data['diff_dict']['indel_new_aa']]
    # New codons for amino acid swaps
    data['diff_dict']['swap_new_codon'] = [back_table[res] for res in data['diff_dict']['swap_new_aa']]
    return



def build_primers(data):
    """
    Pulls 5 codons upstream and downstream of the mutation site and builds
    a 30 or 33 base pair primer.

    Populates the data['primers'] nested dictionary

    TODO
    ----
    Deal with the issue of needing more sequence upstream and downstream of the
    protein coding sequence.
    """

    primers = {'indel':[], 'swap':[]}

    if len(data['diff_dict']['indel_index']) > 0:
        for i in data['diff_dict']['indel_index']:
            upstream = data['dna1']['codon_list'][(i-1)-5:(i-1)]
            downstream = data['dna1']['codon_list'][i:i+5]
            primers['indel'].append(''.join(upstream + downstream))

    if len(data['diff_dict']['swap_index']) > 0:
        for n,site in enumerate(data['diff_dict']['swap_index']):
            # Pull the 5 codons upstream of the swap site from the codon list
            upstream = data['dna1']['codon_list'][(site-1)-5:(site-1)]
            # Pull the 5 codons downstream of the swap site from the codon list
            downstream = data['dna1']['codon_list'][(site+1):(site+1)+5]
            # Add the upstream codons + the new site codon + the downstream codons
            primers['swap'].append(''.join(upstream + list(data['diff_dict']['swap_new_codon'][n]) + downstream))

    data['primers'] = primers['swap'] + primers['indel']
    return
