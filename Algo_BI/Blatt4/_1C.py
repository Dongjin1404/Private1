def count_instance(file_path, pattern):
    with open(file_path, 'r') as file:
        lines = file.readlines()       #reads the file into a list of lines

    instances = {}
    for i in range(0, len(lines), 2):  #two lines at a time: access the reading frame and the sequence
        reading_frame = lines[i].strip().rstrip(':')
        sequence = lines[i+1].strip()
        positions = [j for j in range(len(sequence)) if sequence.startswith(pattern, j)]
        if positions:
            instances[reading_frame] = positions

    return instances

instances = count_instance('_1Btranslated_genome.txt', 'RVLKA')
for reading_frame, positions in instances.items():
    print(f"{reading_frame}: {len(positions)} instances at positions {positions}")

#Reading the file: O(n)
#Iterating over the lines: O(n/2)
#Finding the pattern positions: O(n*k)
    
#For each line, the code iterates over the sequence and checks if the pattern starts at each position.
#For each position in the sequence, the code checks the next k characters to see if they match the pattern.