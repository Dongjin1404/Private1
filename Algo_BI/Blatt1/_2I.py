from _2A import calculate_hamming_distance


def approximatePatternCount(text, pattern, k):
    count = 0
    for i in range(len(text) - len(pattern) + 1):
        if calculate_hamming_distance(text[i : i + len(pattern)], pattern) <= k:
            count += 1
    return count

# function takes in text, pattern, and hammind distance k
# count is initialized to 0
# for loop begins at 0 and ends at len(text) - len(pattern) + 1
# if the hamming distance between the pattern and the text is less than or equal to k, count is incremented by 1
# Yes, allowing for mismatches can significantly change the results
