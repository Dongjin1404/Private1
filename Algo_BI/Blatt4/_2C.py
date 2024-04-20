amino_acid_mass = {
    'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259,
    'F': 147.06841, 'G': 57.02146, 'H': 137.05891, 'I': 113.08406,
    'K': 128.09496, 'L': 113.08406, 'M': 131.04049, 'N': 114.04293,
    'P': 97.05276, 'Q': 128.05858, 'R': 156.10111, 'S': 87.03203,
    'T': 101.04768, 'V': 99.06841, 'W': 186.07931, 'Y': 163.06333,
}

amino_acid_mass = {k: round(v) for k, v in amino_acid_mass.items()}

def cyclospectrum(peptide, amino_acid_mass):
    n = len(peptide)
    extended_peptide = peptide + peptide
    masses = [0]
    for k in range(1, n):          #PEPTIDE LENGTH: Exclude the full peptide
        for i in range(0, n):        #STARTING POSITION: peptide fragments can start at any position in the peptide.
            mass = 0
            for j in range(i, i+k):                         #iterates over a slice of the peptide string
                mass += amino_acid_mass[extended_peptide[j]]#each amino acid is looked up in the amino_acid_mass dictionary
            masses.append(mass)         #the mass of the peptide fragment is stored in the masses dictionary
    masses.append(sum(amino_acid_mass[aa] for aa in peptide))  # Include the mass of the full peptide once
    return sorted(masses)

def linearspectrum(peptide, amino_acid_mass):
    masses = [0]
    for k in range(1, len(peptide)):          #peptide fragments can have lengths from 1 to the length of the peptide-1.
        for i in range(0, len(peptide)-k+1):    #peptide fragments can start at any position in the peptide, but do not wrap around.
            mass = 0
            for j in range(i, i+k):             #iterates over a slice of the peptide string
                mass += amino_acid_mass[peptide[j]]#each amino acid is looked up in the amino_acid_mass dictionary
            masses.append(mass)      #the mass of the peptide fragment is stored in the masses dictionary
    masses.append(sum(amino_acid_mass[aa] for aa in peptide))  # Include the mass of the full peptide once
    return sorted(masses)


    
peptide = 'RVLKA'

cyclo_spectrum = cyclospectrum(peptide, amino_acid_mass)
linear_spectrum = linearspectrum(peptide, amino_acid_mass)

print("Cyclospectrum: ", cyclo_spectrum)
print("Linearspectrum: ", linear_spectrum)
