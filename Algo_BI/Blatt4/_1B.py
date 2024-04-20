genetic_code = {
    'UUU': 'F', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V',
    'UUC': 'F', 'CUC': 'L', 'AUC': 'I', 'GUC': 'V',
    'UUA': 'L', 'CUA': 'L', 'AUA': 'I', 'GUA': 'V',
    'UUG': 'L', 'CUG': 'L', 'AUG': 'M', 'GUG': 'V',
    'UCU': 'S', 'CCU': 'P', 'ACU': 'T', 'GCU': 'A',
    'UCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
    'UCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A',
    'UCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A',
    'UAU': 'Y', 'CAU': 'H', 'AAU': 'N', 'GAU': 'D',
    'UAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
    'UAA': 'Stop', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E',
    'UAG': 'Stop', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
    'UGU': 'C', 'CGU': 'R', 'AGU': 'S', 'GGU': 'G',
    'UGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G',
    'UGA': 'Stop', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
    'UGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'
}

def Protein_translation(rna_pattern, genetic_code):
    Protein = []
    for i in range(0, len(rna_pattern)-2, 3):
        codon = rna_pattern[i:i+3]
        if codon in genetic_code:
            if genetic_code[codon] == 'Stop':
                Protein.append('*')
            else:
                Protein.append(genetic_code[codon])
    return ''.join(Protein)


def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, '') for base in reversed(dna))
#get method of the dictionary, which returns the value for a given key if it exists in the dictionary, and a default value otherwise

def translate_genome(file_path, genetic_code):
    with open(file_path, 'r') as file:
        next(file)  # Skip the first line
        dna_text = ''
        for line in file:
            if line.startswith('>'):
                break  # Stop reading when the second '>' is encountered
            dna_text += line.strip()  # Add the line to dna_text
        rna_text = dna_text.replace('T', 'U')
        reverse_rna_text = reverse_complement(dna_text).replace('T', 'U')
        
    with open('_1Btranslated_genome.txt', 'w') as output:
        for k in range(3):
            protein = Protein_translation(rna_text[k:], genetic_code)
            output.write(f"Reading frame {k+1}:\n{protein}\n")
        for k in range(3):
            protein = Protein_translation(reverse_rna_text[k:], genetic_code)
            output.write(f"Reading frame {k+4} (reverse):\n{protein}\n")

translate_genome('GCF_027925565.1_ASM2792556v1_genomic.fna', genetic_code)

#all operations runtime complexity is O(n)