# Motif_finding_approach

Project Overview

This project implements a novel algorithm for motif discovery in DNA sequences, focusing on detecting recurring nucleotide patterns (motifs) along with their mutated variants. Unlike traditional approaches, this method balances accuracy and efficiency by handling substring lengths and mutation thresholds dynamically.

Key Features
Motif Extraction: Substrings of length 4 to 10 are considered.
Mutation Tolerance:
Length = 4 → allows up to 20% mutation (1 nucleotide)
Length = 5–8 → allows up to 20% or 40% mutation (1–3 nucleotides)
Length = 9–10 → allows up to 20% or 40% mutation

Count Matrix–based Tracking: Each substring and its mutated variants are stored and updated using a dynamic count matrix.
Efficient Memory Use: The algorithm maintains only necessary rows for computation (reducing space complexity).
Backtracking for Alignment: Direction matrix is used to reconstruct all optimal global alignments.

Complexity

Time Complexity (T(n)) ≈ O(n × m)

Space Complexity (S(n)) ≈ O(2 × (m+2)) for scoreback + O((n+2) × (m+2)) for direction matrix
n = length of sequence 1
m = length of sequence 2
