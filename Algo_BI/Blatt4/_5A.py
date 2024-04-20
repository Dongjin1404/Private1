import pandas as pd
from collections import Counter

theoretical_spectrum = [0, 71, 114, 115, 129, 156, 185, 227, 243, 244, 271, 314, 341, 342, 358, 400, 429, 456, 470, 471, 514, 585]
experimental_spectrum = [0, 71, 114, 115, 128, 156, 185, 227, 243, 244, 271, 300, 314, 341, 342, 358, 400, 429, 456, 470, 471, 514, 585]

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


def cyclopeptide_sequencing(spectrum, amino_acid_mass, spectrum_type, use_20_amino_acids=False):

    df = pd.DataFrame(columns=['Iteration', 'Before Expand', 'After Expand', 'After Bound', 'Solutions'])

    valid_amino_acids = spectral_convolution(spectrum)
    print(f"{spectrum_type} - Valid Amino Acids: {valid_amino_acids}")
    print(f"{spectrum_type} - Spectrum: {spectrum}")
    print(f"{spectrum_type} - Amino Acid Masses: {amino_acid_mass}")

    def generate_new_candidates(candidates):      #list of peptides, list of masses
        if use_20_amino_acids:
            return [peptide + [mass] for peptide in candidates for mass in amino_acid_mass.values()]
        else:
            return [peptide + [mass] for peptide in candidates for mass, _ in valid_amino_acids]
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

    def theoretical_linearspectrum(peptide, amino_acid_mass):
        masses = [0]
        for k in range(1, len(peptide)):          #peptide fragments can have lengths from 1 to the length of the peptide-1.
            for i in range(0, len(peptide)-k+1):    #peptide fragments can start at any position in the peptide, but do not wrap around.
                mass = 0
                for j in range(i, i+k):             #iterates over a slice of the peptide string
                    mass += amino_acid_mass[peptide[j]]#each amino acid is looked up in the amino_acid_mass dictionary
                masses.append(mass)      #the mass of the peptide fragment is stored in the masses dictionary
        masses.append(sum(amino_acid_mass[aa] for aa in peptide))  # Include the mass of the full peptide once
        return sorted(masses)
    
    candidates = [[]]
    solutions = []

    


    iteration = 0
    while candidates:
        iteration += 1
        before_expand = len(candidates)
        print(f"{spectrum_type} - Before expand step, number of peptides: {len(candidates)}")

        candidates = generate_new_candidates(candidates) #extends with each possible amino acid mass
        after_expand = len(candidates)
        print(f"{spectrum_type} - After expand step, number of peptides: {len(candidates)}")

        for peptide in candidates[:]:                   #slicing a list without specifying, creates a copy of the list
            # Define a reverse dictionary mapping masses to amino acids
            mass_to_amino_acid = {mass: aa for aa, mass in amino_acid_mass.items()}

            # Convert the masses in peptide back to amino acids
            peptide_amino_acids = [mass_to_amino_acid[mass] for mass in peptide]

            if calculate_total_mass(peptide) == max(spectrum):
                if theoretical_cyclospectrum(peptide_amino_acids, amino_acid_mass) == sorted(spectrum):#theoretical spectrum of current peptide matches given spectrum
                    solutions.append(peptide)
                candidates.remove(peptide)
                print(f"{spectrum_type} - Removed Peptide (Total Mass == Max Spectrum): {peptide}")
            else:
                linear_spectrum = theoretical_linearspectrum(peptide_amino_acids, amino_acid_mass)
                if not all(count <= spectrum.count(mass) for mass, count in Counter(linear_spectrum).items()):
                    candidates.remove(peptide)#checks if the count of a mass in the theoretical linear spectrum is less/equal in the input spectrum.
                    #means that the mass is contained at least as often in the input spectrum as in the theoretical linear spectrum.
                
        after_bound = len(candidates)
        solutions_found = len(solutions)
        print(f"{spectrum_type} - After bound step, number of peptides: {len(candidates)}")
        print(f"{spectrum_type} - Number of solutions found in this iteration: {len(solutions)}")

        # Add a new row to the DataFrame
        df = pd.concat([df, pd.DataFrame([{'Iteration': iteration, 'Before Expand': before_expand, 'After Expand': after_expand, 'After Bound': after_bound, 'Solutions': solutions_found}])], ignore_index=True)


    # Convert the solutions to strings of amino acids
    mass_to_amino_acid = {mass: aa for aa, mass in amino_acid_mass.items()}
    solutions_amino_acids = [''.join(mass_to_amino_acid[mass] for mass in peptide) for peptide in solutions]

    return solutions_amino_acids, df
    

    


import os

# Run the CYCLOPEPTIDESEQUENCING algorithm on the theoretical spectrum
solutions_theoretical, df_theoretical = cyclopeptide_sequencing(theoretical_spectrum, amino_acid_mass, "Theoretical", use_20_amino_acids=False)
print(f"Theoretical Solutions: {solutions_theoretical}")
# Check if the file exists
if not os.path.isfile('theoretical_spectrum.csv'):
    df_theoretical.to_csv('theoretical_spectrum.csv', index=False)
else:  # else it exists so append without writing the header
    df_theoretical.to_csv('theoretical_spectrum.csv', mode='a', header=False, index=False)

# Run the CYCLOPEPTIDESEQUENCING algorithm on the experimental spectrum
solutions_experimental, df_experimental = cyclopeptide_sequencing(experimental_spectrum, amino_acid_mass, "Experimental", use_20_amino_acids=False)
print(f"Experimental Solutions: {solutions_experimental}")
# Check if the file exists
if not os.path.isfile('experimental_spectrum.csv'):
    df_experimental.to_csv('experimental_spectrum.csv', index=False)
else:  # else it exists so append without writing the header
    df_experimental.to_csv('experimental_spectrum.csv', mode='a', header=False, index=False)
