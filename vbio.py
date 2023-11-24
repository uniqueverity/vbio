import sys
from collections import Counter

codonTable = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

aa3_to1_dict = {
    "Ala": "A",  # Alanine
    "Arg": "R",  # Arginine
    "Asn": "N",  # Asparagine
    "Asp": "D",  # Aspartic Acid
    "Cys": "C",  # Cysteine
    "Gln": "Q",  # Glutamine
    "Glu": "E",  # Glutamic Acid
    "Gly": "G",  # Glycine
    "His": "H",  # Histidine
    "Ile": "I",  # Isoleucine
    "Leu": "L",  # Leucine
    "Lys": "K",  # Lysine
    "Met": "M",  # Methionine
    "Phe": "F",  # Phenylalanine
    "Pro": "P",  # Proline
    "Ser": "S",  # Serine
    "Thr": "T",  # Threonine
    "Trp": "W",  # Tryptophan
    "Tyr": "Y",  # Tyrosine
    "Val": "V",  # Valine
}

full_amino_acid_name = {
    "Alanine": "Ala",
    "Arginine": "Arg",
    "Asparagine": "Asn",
    "Aspartic Acid": "Asp",
    "Cysteine": "Cys",
    "Glutamine": "Gln",
    "Glutamic Acid": "Glu",
    "Glycine": "Gly",
    "Histidine": "His",
    "Isoleucine": "Ile",
    "Leucine": "Leu",
    "Lysine": "Lys",
    "Methionine": "Met",
    "Phenylalanine": "Phe",
    "Proline": "Pro",
    "Serine": "Ser",
    "Threonine": "Thr",
    "Tryptophan": "Trp",
    "Tyrosine": "Tyr",
    "Valine": "Val",
}


class Sequence(object):
    """create a valid DNA sequence, for DNA, RNA
    example: seq1 =  Sequence('ATGC')
    """

    def __init__(self, seq=None):
        super(Sequence, self).__init__()
        self.seq = seq

        # To enforce a string storage
        if not isinstance(self.__validate_seq(seq), str):
            raise TypeError(
                "The sequence data given to a Sequence object should"
                "be a sring (not another Sequence)",
                "nor a Non Valid Nucleotide [A,T,G,C,U]",
            )

    def __repr__(self):
        return "Sequence(seq='{}')".format(self.seq)

    def __str__(self):
        return self.seq

    def __validate_seq(self, seq):
        base_nucleotide = ["A", "T", "G", "C", "U"]
        real_seq = seq.upper()
        for base in real_seq:
            if base not in base_nucleotide:
                return False
        return real_seq

    def __len__(self):
        return len(self.seq)

    def __contains__(self, sub_char):
        return sub_char in str(self)

    def __getitem__(self, index):
        if isinstance(index, int):
            return self.seq[index]

    def count(self, subseq, start=0, end=sys.maxsize):
        """return the number of nucleotides in a sequencw"""
        return str(self).count(subseq, start, end)

    def find(self, subseq, start=0, end=sys.maxsize):
        return str(self).find(subseq, start, end)

    def rfind(self, subseq, start=0, end=sys.maxsize):
        return str(self).rfind(subseq, start, end)

    def index(self, subseq, start=0, end=sys.maxsize):
        return str(self).index(subseq, start, end)

    def rindex(self, subseq, start=0, end=sys.maxsize):
        return str(self).rindex(subseq, start, end)

    ####Main function

    def get_symbol_frequency(self):
        base_dict = {"A": 0, "T": 0, "G": 0, "C": 0}
        for base in self.seq:
            if self.__validate_seq != False:
                base_dict[base] += 1
            else:
                return "NucleotideError: {} not a nucleotide ['A, T, G, C']".format(
                    base
                )
            return base_dict

    @property
    def gc(self):
        result = (
            float(str(self.seq).count("G") + str(self.seq).count("C"))
            / len(self.seq)
            * 100
        )
        return result

    @property
    def at(self):
        result = (
            float(str(self.seq).count("A") + str(self.seq).count("T"))
            / len(self.seq)
            * 100
        )
        return result

    def complement(self):
        base_pairs = {"A": "T", "T": "A", "G": "C", "C": "G"}
        comp_pairs = [base_pairs[a] for a in self.seq if a in base_pairs.keys()]
        return "".join(comp_pairs)

    def reverse_complement(self):
        base_pairs = {"A": "T", "T": "A", "G": "C", "C": "G"}
        comp_pairs = [base_pairs[a] for a in self.seq if a in base_pairs.keys()]
        reverse_pairs = "".join(comp_pairs)[::-1]
        return reverse_pairs

    def transcribe(self):
        mrna_result = self.seq.replace("T", "U")
        return mrna_result

    def translate(self, start_pos=0):
        amino_acids_list = [
            codonTable[self.seq[pos : pos + 3]]
            for pos in range(start_pos, len(self.seq) - 2, 3)
        ]
        return "".join(amino_acids_list)
