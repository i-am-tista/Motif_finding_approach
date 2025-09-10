import random
import time

# Function to generate a random DNA sequence
def generate_dna_sequence(length=60):
    return ''.join(random.choices("ACGT", k=length))

# Function to calculate count matrix and determine probable motif
def get_consensus_motif(motif_list, length):
    if not motif_list:
        return None, None
    
    count_matrix = {'A': [0]*length, 'C': [0]*length, 'G': [0]*length, 'T': [0]*length}
    
    for motif in motif_list:
        for i, nucleotide in enumerate(motif):
            if nucleotide in count_matrix:   # ensure valid nucleotide
                count_matrix[nucleotide][i] += 1
    
    probable_motif = "".join(max(count_matrix, key=lambda x: count_matrix[x][i]) for i in range(length))
    return probable_motif, count_matrix

# Function to find motifs with mutation tolerance
def find_motifs(dna_seq, min_length=4, max_length=10):
    potential_motifs = {}
    log_output = ""

    # Loop through lengths 4 to 10
    for length in range(min_length, max_length + 1):
        motif_count = {}

        # Step 3: Generate all substrings of 'length' and count occurrences
        for i in range(len(dna_seq) - length + 1):
            motif = dna_seq[i:i + length]
            motif_count[motif] = motif_count.get(motif, 0) + 1

        # Step 4: Store pot_motif dictionary (all counted motifs)
        log_output += f"\nStep 4: pot_motif dictionary for length {length}:\n{motif_count}\n"

        # Step 5: Apply mutation rules for filtering
        filtered_motifs = {k: v for k, v in motif_count.items() if v > 3}
        if length > 4:
            mutation_threshold = 0.2 if length <= 6 else 0.4
            filtered_motifs = {
                k: v for k, v in motif_count.items()
                if v > 3 or random.random() < mutation_threshold
            }

        # Store motifs that passed filtering
        potential_motifs[length] = list(filtered_motifs.keys())

        # Step 6: Store potential_motif dictionary
        log_output += f"\nStep 6: potential_motif dictionary for length {length}:\n{potential_motifs[length]}\n"

        # Continue processing until we reach length 10 or no new matches are found
        if not potential_motifs[length]:
            continue

        # Expand motifs beyond current length to length+1
        if length < max_length:
            expanded_motifs = []
            for motif in potential_motifs[length]:
                for nucleotide in "ACGT":
                    expanded_motifs.append(nucleotide + motif)  # Prefix expansion
                    expanded_motifs.append(motif + nucleotide)  # Suffix expansion

            # Step 7: Check expanded motifs and continue processing
            new_potential_motifs = []
            for expanded_motif in expanded_motifs:
                if expanded_motif in dna_seq:
                    new_potential_motifs.append(expanded_motif)

            if new_potential_motifs:
                potential_motifs[length + 1] = new_potential_motifs

    # Final output of motifs and counts
    final_motif = None
    final_count_matrix = None
    if potential_motifs:
        for length, motifs in potential_motifs.items():
            probable_motif, count_matrix = get_consensus_motif(motifs, length)
            if probable_motif:  # ensure valid result
                final_motif = probable_motif
                final_count_matrix = count_matrix

    return final_motif, potential_motifs, final_count_matrix, log_output

# Main Execution
if __name__ == "__main__":
    start_time = time.time()

    random_dna = generate_dna_sequence(random.randint(50, 60))
    probable_motif, all_motifs, count_matrix, log_output = find_motifs(random_dna)

    # Stop clock and calculate execution time
    end_time = time.time()
    execution_time = end_time - start_time

    # Prepare final output for file
    final_output = (
        f"\n\n========== New Run ==========\n"
        f"Generated DNA Sequence:\n{random_dna}\n"
        f"{log_output}"
        f"\nStep 8: Potential motif from given sequence: {probable_motif}\n"
        f"Step 9: All motifs found: {all_motifs}\n"
        f"\nStep 10: Execution Time: {execution_time:.5f} seconds\n"
        f"===================================\n"
    )

    # Append all outputs to text file
    with open("motif_results.txt", "a") as file:
        file.write(final_output)

    print("All results have been saved to 'motif_results.txt'.")
    print(final_output)  # Also print on console for quick check
