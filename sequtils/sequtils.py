import sys
from collections import Counter
codonTable = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

aa3_to1_dict = {
    'Ala': 'A',  # Alanine
    'Arg': 'R',  # Arginine
    'Asn': 'N',  # Asparagine
    'Asp': 'D',  # Aspartic Acid
    'Cys': 'C',  # Cysteine
    'Gln': 'Q',  # Glutamine
    'Glu': 'E',  # Glutamic Acid
    'Gly': 'G',  # Glycine
    'His': 'H',  # Histidine
    'Ile': 'I',  # Isoleucine
    'Leu': 'L',  # Leucine
    'Lys': 'K',  # Lysine
    'Met': 'M',  # Methionine
    'Phe': 'F',  # Phenylalanine
    'Pro': 'P',  # Proline
    'Ser': 'S',  # Serine
    'Thr': 'T',  # Threonine
    'Trp': 'W',  # Tryptophan
    'Tyr': 'Y',  # Tyrosine
    'Val': 'V',  # Valine
}

full_amino_acid_name = {
    'Alanine': 'Ala',
    'Arginine': 'Arg',
    'Asparagine': 'Asn',
    'Aspartic Acid': 'Asp',
    'Cysteine': 'Cys',
    'Glutamine': 'Gln',
    'Glutamic Acid': 'Glu',
    'Glycine': 'Gly',
    'Histidine': 'His',
    'Isoleucine': 'Ile',
    'Leucine': 'Leu',
    'Lysine': 'Lys',
    'Methionine': 'Met',
    'Phenylalanine': 'Phe',
    'Proline': 'Pro',
    'Serine': 'Ser',
    'Threonine': 'Thr',
    'Tryptophan': 'Trp',
    'Tyrosine': 'Tyr',
    'Valine': 'Val',
}

def gc_content(self):
          result = float(str(self.seq).count('G') + str(self.seq).count('C'))/len(self.seq) * 100
          return result
    
def at_content(self):
          result = float(str(self.seq).count('A') + str(self.seq).count('T'))/len(self.seq) * 100
          return result  

def __get_key(val,my_dict):
     for key,value in my_dict.items():
        if val == value:
            return key
                
def __get_value(val,my_dict):
     for key,value in my_dict.items():
          if val == key:
               return value 
    
def convert1to3(seq):
    term_list = []
    for i in seq:
        res =__get_key(i,aa3_to1_dict)
        term_list.append(res)
    return "".join(term_list)
    
def __kmers(seq,k=2):
    pair_list = []
    for i in range (0,len(seq),k): 
        pair_list.append(seq[i:i+k])
    return pair_list

def convert3to1(seq):
    term_list = []
    for i in __kmers(seq,k=3):
        res = __get_value(i,aa3_to1_dict)
        term_list.append(res)

    return "".join(term_list)
    
def get_kmers(seq,k=2):
    pair_list = []
    for i in range(0,len(seq),k):
        pair_list.append(str(seq)[i:i+k])
        return pair_list
            
    
    
def hamming(lhs,rhs):
     return len([(x,y) for x,y in zip(lhs,rhs) if x!=y])
    
def occurence(main_seq,sub_seq):
    start = 0
    indices = []
    while True:
        start = main_seq.find(sub_seq,start)
        if start > 0:
            indices.append(start)
        else:
            break
        start +=1
    return indices
          
def get_acid_name(seq):
    term_list = []
    for i in __kmers(seq,k=3):
        res = __get_key(i,full_amino_acid_name)
        term_list.append(res)
        return "".join(term_list)
    
def codon_frequency(seq, aminoacid):
    tmpList = []
    for i in range(0, len(seq) - 2, 3):
          if codonTable[seq[i:i +3]] == aminoacid:
                tmpList.append(seq[i:i + 3])
                
    freqDict = dict(Counter(tmpList))
    totalScore = sum(freqDict.values())
    for seq in freqDict:
        freqDict[seq] = round(freqDict[seq] / totalScore, 2)
    return freqDict