import gzip
import numpy as np

import pandas as pd

from .constants import codon_dict


watson_crick = {'A': 'T',
                'T': 'A',
                'C': 'G',
                'G': 'C',
                'U': 'A',
                'N': 'N'}

watson_crick.update({k.lower(): v.lower()
                     for k, v in watson_crick.items()})

codon_maps = {('e_coli', 0.12): 
 {'*': ['TAA', 'TGA'],
  'A': ['GCA', 'GCC', 'GCG', 'GCT'],
  'C': ['TGC', 'TGT'],
  'D': ['GAC', 'GAT'],
  'E': ['GAA', 'GAG'],
  'F': ['TTC', 'TTT'],
  'G': ['GGC', 'GGG', 'GGT'],
  'H': ['CAC', 'CAT'],
  'I': ['ATC', 'ATT'],
  'K': ['AAA', 'AAG'],
  'L': ['CTG', 'TTA', 'TTG'],
  'M': ['ATG'],
  'N': ['AAC', 'AAT'],
  'P': ['CCA', 'CCG', 'CCT'],
  'Q': ['CAA', 'CAG'],
  'R': ['CGC', 'CGT'],
  'S': ['AGC', 'AGT', 'TCC', 'TCG', 'TCT'],
  'T': ['ACA', 'ACC', 'ACG', 'ACT'],
  'V': ['GTA', 'GTC', 'GTG', 'GTT'],
  'W': ['TGG'],
  'Y': ['TAC', 'TAT']}}


def read_fasta(f, as_df=False):
    if f.endswith('gz'):
        fh = gzip.open(f)
        txt = fh.read().decode()
    else:
        fh = open(f, 'r')
        txt = fh.read()
    fh.close()
    records = parse_fasta(txt)
    if as_df:
        return pd.DataFrame(records, columns=('name', 'seq'))
    else:
        return records


def parse_fasta(txt):
    entries = []
    txt = '\n' + txt.strip()
    for raw in txt.split('\n>'):
        name = raw.split('\n')[0].strip()
        seq = ''.join(raw.split('\n')[1:]).replace(' ', '')
        if name:
            entries += [(name, seq)]
    return entries


def write_fasta(filename, list_or_records):
    if isinstance(list_or_records, pd.DataFrame) and list_or_records.shape[1] == 2:
        list_or_records = list_or_records.values
    list_or_records = list(list_or_records)
    with open(filename, 'w') as fh:
        fh.write(format_fasta(list_or_records))


def format_fasta(list_or_records):
    if len(list_or_records) == 0:
        records = []
    elif isinstance(list_or_records[0], str):
        records = list_to_records(list_or_records)
    else:
        records = list_or_records
    
    lines = []
    for name, seq in records:
        lines.extend([f'>{name}', str(seq)])
    return '\n'.join(lines)


def list_to_records(xs):
    n = len(xs)
    width = int(np.ceil(np.log10(n)))
    fmt = '{' + f':0{width}d' + '}'
    records = []
    for i, s in enumerate(xs):
        records += [(fmt.format(i), s)]
    return records


def translate_dna(s):
    assert len(s) % 3 == 0, 'length must be a multiple of 3'
    return ''.join([codon_dict[s[i*3:(i+1)*3]] for i in range(int(len(s)/3))])


def reverse_complement(seq):
    return ''.join(watson_crick[x] for x in seq)[::-1]


def reverse_translate_random(aa_seq, organism='e_coli', rs='input', cutoff=0.12):
    if rs == 'input':
        seed = hash(aa_seq) % 10**8
        rs = np.random.RandomState(seed=seed)
    if (organism, cutoff) not in codon_maps:
        raise NotImplementedError
    codon_map = codon_maps[(organism, cutoff)]
    return ''.join([rs.choice(codon_map[x]) for x in aa_seq])

