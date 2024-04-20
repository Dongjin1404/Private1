amino_acid_mass = {
    'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259,
    'F': 147.06841, 'G': 57.02146, 'H': 137.05891, 'I': 113.08406,
    'K': 128.09496, 'L': 113.08406, 'M': 131.04049, 'N': 114.04293,
    'P': 97.05276, 'Q': 128.05858, 'R': 156.10111, 'S': 87.03203,
    'T': 101.04768, 'V': 99.06841, 'W': 186.07931, 'Y': 163.06333,
}
amino_acid_mass = {k: round(v) for k, v in amino_acid_mass.items()}

theoretical_spectrum = [0, 71, 114, 115, 129, 156, 185, 227, 243, 244, 271, 314, 341, 342, 358, 400, 429, 456, 470, 471, 514, 585]
experimental_spectrum = [0, 71, 114, 115, 128, 156, 185, 227, 243, 244, 271, 300, 314, 341, 342, 358, 400, 429, 456, 470, 471, 514, 585]


def spectral_convolution(spectrum):
    convolution = {}
    for i in range(len(spectrum)):
        for j in range(i):
            diff = spectrum[i] - spectrum[j]
            if 57 <= diff <= 200:  # the mass of the amino acids is between 57 and 200 Da
                convolution[diff] = 1 if diff not in convolution else convolution[diff] + 1
    # Create a list of tuples (mass, multiplicity)
    result = [(mass, multiplicity) for mass, multiplicity in sorted(convolution.items(), key=lambda x: x[1], reverse=True)]
    return result

# Calculate and print spectral convolution of theoretical_spectrum
convolution_t = spectral_convolution(theoretical_spectrum)
for mass, multiplicity in convolution_t:
    print(f"Theoretical - Mass: {mass}, Multiplicity: {multiplicity}")

# Calculate and print spectral convolution of experimental_spectrum
convolution_e = spectral_convolution(experimental_spectrum)
for mass, multiplicity in convolution_e:
    print(f"Experimental - Mass: {mass}, Multiplicity: {multiplicity}") 

#A higher multiplicity indicates that the corresponding amino acid is more likely to be part of the original peptide sequence.
    
#For the theoretical spectrum, the amino acids with the highest multiplicity are A, N, D, E, and R, all with a multiplicity of 8.
#These would be the primary candidates for further reconstruction.
    
#For the experimental spectrum, the amino acids with the highest multiplicity are D and R, both with a multiplicity of 9,
#followed by A and E with a multiplicity of 8. 
#These would be the primary candidates for further reconstruction.
    
#In reality, it involves considering all possible sequences of the candidate amino acids and comparing them to the experimental spectrum. 
#Other factors, such as the possibility of modifications to the amino acids, can also complicate the problem.