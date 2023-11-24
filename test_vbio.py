from vbio import __version__
import vbio
from vbio.sequtils import gc_content,at_content


def test_version():
    assert __version__ == "0.0.1"

def test_issequence():
    seq1 = vbio.Sequence("ATGCTATGCTT")
    assert isinstance(seq1,type(seq1))

def test_complement():
    seq1 = vbio.Sequence("ATGCTATGCTT")
    res_complement = seq1.complement()
    assert res_complement == 'TACGATACGAA'

def test_reverse_complement():
    seq1 = vbio.Sequence("ATGCTATGCTT")
    res_complement = seq1.reverse_complement()
    assert res_complement =="AAGCATAGCAT"

def test_translate():
    seq1 = vbio.Sequence("ATGCTATGCTT")
    result = seq1.translate()
    assert result == 'MLC'

def test_transcription():
    seq1 = vbio.Sequence("ATGCTATGCTT")
    result = seq1.transcribe()
    assert result == 'AUGCUAUGCUU'

def test_nucleotide_freq():
    seq1 = vbio.Sequence('ATGCTATGCTT')
    result = seq1.get_symbol_frequency()
    assert result == {'A':1, 'T': 0, 'G': 0, 'C': 0}

def test_gc_content():
    seq1 = vbio.Sequence('ATGCTATGCTT')
    result = gc_content(seq1)
    assert result == 36.36363636363637
    

