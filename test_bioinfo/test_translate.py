
from __future__ import print_function

from joker.bioinfo.translate import CodonTable


def remove_whitespaces(s):
    return ''.join(s.split())


dna = """
    TTGATGGCTAAGAGTAAAATCTTAAAAAACACACTGGTTCTATATTTTCGTCAAGTTTTG
    ATTGTATTAATTACTCTCTATTCAATGAGAGTTGTATTAAATGAATTAGGTGTGGATGAT
    TTTGGTATTTATAGCGTTGTGGCTGGTTTTGTAACTTTACTTGCATTTTTACCCGGAAGC
    ATGGCGAGTGCAACGCAGCGGTTTTTCTCTTTTGCGATGGGGAAATCGGATATAGTAAAA
    TTAAAGCAAACCTTCAGTGTTAATTTAGTTATGTATACTGGCATAGCCTTGTTAGCATAT
    ATAACATTTCAAACTATCGGATTTTGGTATGTTGATGAATATCTAAAAATACCTCATAAC
    CGCTTTCATGCAGCCTTGGAATTATATCACTATGTGTCTTTATCATTTATTTTTTCAATT
    TTTTCTGCGCCTTTTATCGCGATTTTAATTGCGCACGAAGATATGCACATTTATGCGATC
    TTCTCGGTTTTTGATGCATTTTTAAAACTAGTAGCCGCAATTTCTTTAGACTATGTGAAC
    TATGATTTGTTAGCTTATTATGGAGTGGCTCTTTTGATTGTATCTGGATTGCTTGCTTTC
    GCGTATATTTTTATATGTATAAAGAAATATCCAGAGTGTCAAATGAAAAAGCTTTATTGG
    AGTTCGAGTATACTGAAAGAAATTATTGGTTTCACGATATGGACTTTGCTAGGTCAATTG
    AGCACTGTTTTTAGAAATCAGGCAGTAACTGTTCTTGTAAACCAAATGTTTAATCCTTCA
    ATTGTGGCAGCTCGTGCAATTGCCTTGAATGTGGCTAGTCAAGTTGGAATTTTTTCGAAT
    AATTTAAATACAGGGTTATATCCACCAATTATAAAAGCTTACGCAGCAAATCAAAAAGAG
    GAAATGCTGAGTTTAGTTTTTAATGGTTCTAAACTGACTTTCTTCTTGATGTGGGTATGT
    GCATTACCCATGTTGCTTGAAATGGAAACGATATTAACACTTTGGCTAAAAACACCACCA
    TCAGAAGCGATATTATTTACTCAGTTAGCGATTGTTGAATCCTTGATACTGGCTATAAGC
    ATGCCTTTAACTACTGCGGCAAGAGCACCAGGGAAAATGGCGTTGTATGAAATAATTCTA
    GGCGCCATTCAAATAGCAATATTTTTTGTTTCATGGTGGTTTCTGAAATTTGGTTATTCT
    GCTGAATGGGTTTTCTACATAGCAATAGCGGCTAATATTATTATGTTTAAAGTTCGCTTG
    TTATTAGTAAGCTATTTAACTGATCTTCCTGTAATGGCTTATTATCAAAGAGTTGTAGTT
    CCGGTTTTATCTGTTCTGGTTATTTCATCATTATCTAGTATCTTGCTATTGAACAACTTA
    CCAAGAAATTTAGGTGCATCTTTATTGGTGATTATTTTTTCTATCGGTGTTTCGATATTG
    ACAATGTACTTCTTAGGCTTAGATAAGTACTGGCGTGAAAAAGTGGTCGGTGTGCTAACC
    AGTAAATTTTTAAAATCTAGAGAGGTGCGATGA"""

protein = """
    LMAKSKILKNTLVLYFRQVLIVLITLYSMRVVLNELGVDDFGIYSVVAGFVTLLAFLPGSMASATQRFFS
    FAMGKSDIVKLKQTFSVNLVMYTGIALLAYITFQTIGFWYVDEYLKIPHNRFHAALELYHYVSLSFIFSI
    FSAPFIAILIAHEDMHIYAIFSVFDAFLKLVAAISLDYVNYDLLAYYGVALLIVSGLLAFAYIFICIKKY
    PECQMKKLYWSSSILKEIIGFTIWTLLGQLSTVFRNQAVTVLVNQMFNPSIVAARAIALNVASQVGIFSN
    NLNTGLYPPIIKAYAANQKEEMLSLVFNGSKLTFFLMWVCALPMLLEMETILTLWLKTPPSEAILFTQLA
    IVESLILAISMPLTTAARAPGKMALYEIILGAIQIAIFFVSWWFLKFGYSAEWVFYIAIAANIIMFKVRL
    LLVSYLTDLPVMAYYQRVVVPVLSVLVISSLSSILLLNNLPRNLGASLLVIIFSIGVSILTMYFLGLDKY
    WREKVVGVLTSKFLKSREVR"""


dna = remove_whitespaces(dna)
protein = remove_whitespaces(protein)


def test_translate():
    ctab = CodonTable.load()
    assert len(dna) % 3 == 0
    result = ctab.translate(dna)
    print(result)
    print(protein)

    assert result[:-1] == protein
    assert result[-1] == '*'


if __name__ == '__main__':
    test_translate()
