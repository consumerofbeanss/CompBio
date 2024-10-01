def transcription(dna_seq):
    dna_seq = dna_seq.upper()
    list = []
    for i in dna_seq:
        list.append(i)

    for i in range(len(list)):
        if list[i] == 'A':
            list[i] = 'T'
        elif list[i] == 'T':
            list[i] ='A'
        elif list[i] == 'G':
            list[i] = 'C'
        elif list[i] == 'C':
            list[i] = 'G'

    complement = "".join(list)

    mlist = list.copy()
    for i in range(len(mlist)):
        if mlist[i] == 'T':
            mlist[i] = 'U'

    mRNA = "".join(mlist)

    return complement, mRNA

def translation(mRNA):
    codonDict = {
        'UUU': ['Phe', 'F'], 'UUC': ['Phe', 'F'],  # Phenylalanine
        'UUA': ['Leu', 'L'], 'UUG': ['Leu', 'L'],  # Leucine
        'CUU': ['Leu', 'L'], 'CUC': ['Leu', 'L'],
        'CUA': ['Leu', 'L'], 'CUG': ['Leu', 'L'],  # Leucine (more codons)
        'AUU': ['Ile', 'I'], 'AUC': ['Ile', 'I'], 'AUA': ['Ile', 'I'],  # Isoleucine
        'AUG': ['Met', 'M'],  # Methionine (start codon)
        'GUU': ['Val', 'V'], 'GUC': ['Val', 'V'], 'GUA': ['Val', 'V'], 'GUG': ['Val', 'V'],  # Valine
        'UCU': ['Ser', 'S'], 'UCC': ['Ser', 'S'], 'UCA': ['Ser', 'S'], 'UCG': ['Ser', 'S'],  # Serine
        'AGU': ['Ser', 'S'], 'AGC': ['Ser', 'S'],  # Serine (more codons)
        'CCU': ['Pro', 'P'], 'CCC': ['Pro', 'P'], 'CCA': ['Pro', 'P'], 'CCG': ['Pro', 'P'],  # Proline
        'ACU': ['Thr', 'T'], 'ACC': ['Thr', 'T'], 'ACA': ['Thr', 'T'], 'ACG': ['Thr', 'T'],  # Threonine
        'GCU': ['Ala', 'A'], 'GCC': ['Ala', 'A'], 'GCA': ['Ala', 'A'], 'GCG': ['Ala', 'A'],  # Alanine
        'UAU': ['Tyr', 'Y'], 'UAC': ['Tyr', 'Y'],  # Tyrosine
        'UAA': ['Stop', 'Stop'], 'UAG': ['Stop', 'Stop'], 'UGA': ['Stop', 'Stop'],  # Stop codons
        'CAU': ['His', 'H'], 'CAC': ['His', 'H'],  # Histidine
        'CAA': ['Gln', 'Q'], 'CAG': ['Gln', 'Q'],  # Glutamine
        'AAU': ['Asn', 'N'], 'AAC': ['Asn', 'N'],  # Asparagine
        'AAA': ['Lys', 'K'], 'AAG': ['Lys', 'K'],  # Lysine
        'GAU': ['Asp', 'D'], 'GAC': ['Asp', 'D'],  # Aspartic Acid
        'GAA': ['Glu', 'E'], 'GAG': ['Glu', 'E'],  # Glutamic Acid
        'UGU': ['Cys', 'C'], 'UGC': ['Cys', 'C'],  # Cysteine
        'UGG': ['Trp', 'W'],  # Tryptophan
        'CGU': ['Arg', 'R'], 'CGC': ['Arg', 'R'], 'CGA': ['Arg', 'R'], 'CGG': ['Arg', 'R'],  # Arginine
        'AGA': ['Arg', 'R'], 'AGG': ['Arg', 'R'],  # Arginine (more codons)
        'GGU': ['Gly', 'G'], 'GGC': ['Gly', 'G'], 'GGA': ['Gly', 'G'], 'GGG': ['Gly', 'G']  # Glycine
    }

    alist = []
    for i in range(0, len(mRNA), 3):
        codon = mRNA[i:i+3]
        aminoAcidCodon = codonDict.get(codon)

        if aminoAcidCodon:
            aminoAcid, symbol = aminoAcidCodon
            alist.append(f'{aminoAcid} ({symbol})')

        else:
            alist.append(f'Invalid codon: ({codon})')

    return ' - '.join(alist)

dnaSeq = input("Input DNA = ")

while len(dnaSeq) % 3 != 0:
    print("\nDNA should be multiple of 3")
    dnaSeq = input("Reinput DNA = ")

complement, mRNA = transcription(dnaSeq)
aminoAcidSeq = translation(mRNA)

print("Complement = ", complement)
print("mRNA = ", mRNA)
print("Amino Acid Sequence = ", aminoAcidSeq)






