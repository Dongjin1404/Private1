def profile_entropy(profile):
    """Task 2c: Calculate the entropy of a profile

    Arguments:
        profile (matrix): A matrix where entry ij is the probability of
            nucleotide i being at position j.

    Returns:
        float: The entropy
    """
    return 0.0

import numpy as np

def profile_entropy(profile):
    # Calculate the entropy for each position and sum them
    entropy = -np.sum(profile * np.log2(profile + np.finfo(float).eps)) # Add epsilon to avoid log(0) which is undefined
    # If a probability is zero, it becomes np.finfo(float).eps before taking the log, so there's no error or undefined value.
    return entropy

# The entropy of a profile matrix is a measure of the uncertainty or randomness of the nucleotides at each position.
# For each position (column), the entropy is the sum of -p*log2(p) for all nucleotides

profile = np.array([[0.2, 0.6, 0.2, 0.0],
                    [0.1, 0.1, 0.7, 0.1],
                    [0.0, 0.2, 0.0, 0.8],
                    [0.7, 0.1, 0.1, 0.1]])

print(profile_entropy(profile))