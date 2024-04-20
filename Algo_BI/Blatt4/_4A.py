amino_acid_mass = {
    'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259,
    'F': 147.06841, 'G': 57.02146, 'H': 137.05891, 'I': 113.08406,
    'K': 128.09496, 'L': 113.08406, 'M': 131.04049, 'N': 114.04293,
    'P': 97.05276, 'Q': 128.05858, 'R': 156.10111, 'S': 87.03203,
    'T': 101.04768, 'V': 99.06841, 'W': 186.07931, 'Y': 163.06333,
}
amino_acid_mass = {k: round(v) for k, v in amino_acid_mass.items()}

def spectral_convolution(spectrum):
    convolution = {}
    for i in range(len(spectrum)):
        for j in range(i):
            diff = spectrum[i] - spectrum[j]
            if 57 <= diff <= 200:  # the mass of the amino acids is between 57 and 200 Da
                convolution[diff] = 1 if diff not in convolution else convolution[diff] + 1
    # Create a list of tuples (mass, multiplicity)
    result = [(mass, multiplicity) for mass, multiplicity in sorted(convolution.items(), key=lambda x: x[1], reverse=True)]
    return result[:5]  # Return only the first 5 elements (A,N,D,E,R)


def cyclopeptide_sequencing(spectrum, amino_acid_mass, use_20_amino_acids=False):

    valid_amino_acids = spectral_convolution(spectrum)

    def generate_new_candidates(candidates):      #list of peptides, list of masses
        if use_20_amino_acids:
            return [peptide + [mass] for peptide in candidates for mass in amino_acid_mass.values()]
        else:
            return [peptide + [mass] for peptide in candidates for mass in valid_amino_acids.values()]
        #extends each peptide in candidates with each possible amino acid mass
        #[[57, 57], [57, 71], [57, 87], [57, 97], [71, 57], [71, 71], [71, 87],...

    def calculate_total_mass(peptide):          #list of masses,
        return sum(peptide)

    def theoretical_cyclospectrum(peptide, amino_acid_mass):
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

    candidates = [[]]
    solutions = []

    while candidates:
        print(f"Before expand step, number of peptides: {len(candidates)}")

        candidates = generate_new_candidates(candidates) #extends with each possible amino acid mass
        print(f"After expand step, number of peptides: {len(candidates)}")

        for peptide in candidates[:]:                   #slicing a list without specifying, creates a copy of the list
            # Define a reverse dictionary mapping masses to amino acids
            mass_to_amino_acid = {mass: aa for aa, mass in amino_acid_mass.items()}

            # Convert the masses in peptide back to amino acids
            peptide_amino_acids = [mass_to_amino_acid[mass] for mass in peptide]

            if calculate_total_mass(peptide) == max(spectrum):
                if theoretical_cyclospectrum(peptide_amino_acids, amino_acid_mass) == sorted(spectrum):#theoretical spectrum of current peptide matches given spectrum
                    solutions.append(peptide)
                candidates.remove(peptide)
            elif not set(theoretical_cyclospectrum(peptide_amino_acids, amino_acid_mass)).issubset(spectrum):
                candidates.remove(peptide) #checks if the theoretical spectrum of the current peptide is not a subset of the given spectrum

        print(f"After bound step, number of peptides: {len(candidates)}")
        print(f"Number of solutions found in this iteration: {len(solutions)}")

    return solutions

#branch-and-bound algorithm
#a branching step to increase the number of candidate solutions
#bounding step to remove hopeless candidates
#this algorithm has not been proven to be polynomial


#The CYCLOPEPTIDESEQUENCING algorithm, which is used to reconstruct a peptide sequence from its spectrum.

#Spectral Convolution:
#This function calculates the differences between all pairs of masses in the spectrum that are within the range of 57 to 200,
#and counts the occurrences of each difference. 

#The generate_new_candidates function:
#generates new candidates by extending each current candidate with each possible amino acid mass.

#The theoretical_cyclospectrum function:
#generates all possible subpeptides of the peptide, calculates their masses, and returns a sorted list of these masses.

#The calculate_total_mass function:
#calculates the total mass of a peptide by summing its amino acid masses.

#INITIALIZE: 
#an empty list of candidates and an empty list of solutions.
#Each candidate is a list of amino acid masses representing a possible peptide
#each solution is a list of amino acid masses representing a peptide that matches the given spectrum.

#MAIN LOOP: continues until there are no more candidates. 

#Expand Step: The function generates new candidates by extending each current candidate with each possible amino acid mass.

#Bound Step: checks each candidate to see if it matches the given spectrum
#If the total mass of the candidate is equal to the maximum mass in the spectrum, 
#and the theoretical cyclospectrum of the candidate matches the given spectrum
#candidate is added to the solutions and removed from the candidates.

#If the theoretical cyclospectrum of the candidate is not a subset of the given spectrum, 
#the candidate is removed from the candidates.

#SO candidate peptide is kept if the theoretical cyclospectrum of the candidate is a subset of the given spectrum,
#and the total mass of the candidate is not equal to the maximum mass in the spectrum.