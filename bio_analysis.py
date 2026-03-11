import re

# Code génétique pour traduction ADN -> protéine
genetic_code = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

# ----------------------------
# Fonctions ADN
# ----------------------------
def reverse_complement(dna):
    complement = str.maketrans('ATCG', 'TAGC')
    return dna.translate(complement)[::-1]

def translate_dna(dna_seq):
    protein = ''
    for i in range(0, len(dna_seq)-2, 3):
        codon = dna_seq[i:i+3]
        protein += genetic_code.get(codon,'X')
    return protein

def find_orfs(dna):
    """Trouver tous les ORFs sur les 2 brins et 3 cadres par brin"""
    orfs = []
    strands = {'Forward': dna, 'Reverse': reverse_complement(dna)}

    for strand_name, seq in strands.items():
        for frame in range(3):
            for start in range(frame, len(seq)-2, 3):
                if seq[start:start+3] == 'ATG':  # codon start
                    for i in range(start, len(seq)-2, 3):
                        codon = seq[i:i+3]
                        if codon in ['TAA','TAG','TGA']:
                            orf_seq = seq[start:i+3]
                            protein = translate_dna(orf_seq)
                            orfs.append({
                                'strand': strand_name,
                                'frame': frame,
                                'start': start,
                                'stop': i+3,
                                'length': len(orf_seq),
                                'sequence': orf_seq,
                                'protein': protein
                            })
                            break
    return orfs

# ----------------------------
# Motifs
# ----------------------------
def find_promoters(dna):
    """Recherche de Box -35 et Box -10"""
    promoters = []
    for m35 in re.finditer(r'TTGACA', dna):
        for m10 in re.finditer(r'TATAAT', dna[m35.end():]):
            if 15 <= (m10.start() - m35.end()) <= 20:
                promoters.append({
                    'box35': (m35.start(), m35.group()),
                    'box10': (m10.start()+m35.end(), m10.group())
                })
    return promoters

def find_rho_independent_terminators(dna):
    """Recherche de séquences palindromiques suivies de polyT"""
    terminators = []
    for i in range(len(dna)-8):
        for l in range(4, 8):  # longueur palindrome
            seq1 = dna[i:i+l]
            seq2 = reverse_complement(seq1)
            if dna[i+l:i+2*l] == seq2:
                # vérifier si suivi de polyT
                j = i+2*l
                polyT = re.match(r'T{4,}', dna[j:j+10])
                if polyT:
                    terminators.append({
                        'start': i,
                        'sequence': dna[i:j+len(polyT.group())]
                    })
    return terminators

def find_shine_dalgarno(dna):
    """Recherche de sites Shine-Dalgarno avant ATG (7-9 nt)"""
    sd_sites = []
    for m in re.finditer('ATG', dna):
        start_atg = m.start()
        upstream = dna[max(0, start_atg-12):start_atg]  # zone possible
        for sd in re.finditer('AGGAGG', upstream):
            distance = start_atg - (max(0,start_atg-12)+sd.start())
            if 7 <= distance <= 9:
                sd_sites.append({
                    'position': max(0,start_atg-12)+sd.start(),
                    'sequence': sd.group(),
                    'atg_position': start_atg
                })
    return sd_sites