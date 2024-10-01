def aminoAcidCombination(aminoAcid):
    aminoAcidDict = {
        'F': ['UUU', 'UUC'],                                       # Phenylalanine
        'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],           # Leucine
        'I': ['AUU', 'AUC', 'AUA'],                                # Isoleucine
        'M': ['AUG'],                                              # Methionine (start codon)
        'V': ['GUU', 'GUC', 'GUA', 'GUG'],                         # Valine
        'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],           # Serine
        'P': ['CCU', 'CCC', 'CCA', 'CCG'],                         # Proline
        'T': ['ACU', 'ACC', 'ACA', 'ACG'],                         # Threonine
        'A': ['GCU', 'GCC', 'GCA', 'GCG'],                         # Alanine
        'Y': ['UAU', 'UAC'],                                       # Tyrosine
        'Stop': ['UAA', 'UAG', 'UGA'],                             # Stop codons
        'H': ['CAU', 'CAC'],                                       # Histidine
        'Q': ['CAA', 'CAG'],                                       # Glutamine
        'N': ['AAU', 'AAC'],                                       # Asparagine
        'K': ['AAA', 'AAG'],                                       # Lysine
        'D': ['GAU', 'GAC'],                                       # Aspartic Acid
        'E': ['GAA', 'GAG'],                                       # Glutamic Acid
        'C': ['UGU', 'UGC'],                                       # Cysteine
        'W': ['UGG'],                                              # Tryptophan
        'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],           # Arginine
        'G': ['GGU', 'GGC', 'GGA', 'GGG']                          # Glycine
    }

    codonComb = []
    for i in aminoAcid:
        codonComb.append(aminoAcidDict[i])

    def cartesianProduct(list):
        if len(list) == 0:
            return [[]]  # Base case: return empty list for the Cartesian product of 0 lists
        result = []
        for item in list[0]:
            for rest in cartesianProduct(list[1:]):
                result.append([item] + rest)
        return result

    combinations = cartesianProduct(codonComb)

    for combination in combinations:
        mRNA = "".join(combination)
        codons = " - ".join(combination)
        count = {}

        for codon in combination:
            if codon in count:
                count[codon] += 1
            else:
                count[codon] = 1

        print("mRNA = ", mRNA)
        print("codons = ", codons)
        for codon, count in count.items():
            print(codon, " = ", count)
        print("\n")


aminoAcid = input("Input Amino Acid Letters (max 3) = ")

while len(aminoAcid) > 3:
    print("\nMax 3 amino acids ")
    aminoAcid = input("Reinput Amino Acid Letters = ")

aminoAcidCombination(aminoAcid)


