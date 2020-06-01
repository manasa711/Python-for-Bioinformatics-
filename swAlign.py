#!/usr/bin/env python3

import sys

#storing the arguments into variables

seq_file1 = sys.argv[1]

seq_file2 = sys.argv[2]

#Retrieving sequences from the fasta files
seq1 = ""

seq2 = ""

with open(seq_file1) as fhand1:
    for line in fhand1:
        line.rstrip()
        if line.startswith(">"):
            continue
        seq1 = seq1 + line.strip("\n")

with open(seq_file2) as fhand2:
    for line in fhand2:
        line.rstrip()
        if line.startswith(">"):
            continue
        seq2 = seq2 + line.strip("\n")

seq1 = list(seq1)
seq2 = list(seq2)

#print(seq1,seq2)

len_seq1 = len(seq1) #getting the length of sequence 1
len_seq2 = len(seq2) #getting the lenght of sequence 2

m = len_seq1+1
n = len_seq2+1

#1. initializing a matrix using lists
#and filling it up with 0s
score_mat = [[0 for i in range(n)] for i in range(m)]

#filling the first row with 0s
for i in range(m):
    score_mat[i][0] = 0

#filling the first column with 0s
for j in range(n):
    score_mat[0][j] = 0

match = 1
mismatch = -1
gap = -1

#Function for returning match or mismatch
def match_vs_mismatch(n1,n2):
    match = 1
    mismatch = -1
    if n1 == n2:
        return match
    else:
        return mismatch

#Function for printing the matrix
def print_mat(matrix_name):
    for row in score_mat:
        for i in row:
            print(i, end = "\t")
        print("\n")

#2. Matrix Filling
#filling up the cells based on matches, mismatches and gaps
for i in range(1,m):
    for j in range(1,n):
        diagonal_score = score_mat[i-1][j-1] + match_vs_mismatch(seq1[i-1],seq2[j-1])
        leftcell_score = score_mat[i][j-1] + gap
        topcell_score = score_mat[i-1][j] + gap
        great_1 = max(diagonal_score,leftcell_score,topcell_score)
        if great_1 < 0: #checking if the maximum value from the left, top and diagonal cell are greater than zero, if not, making the value to be filled as 0
            score_mat[i][j] = 0
        else:
            score_mat[i][j] = great_1

#print_mat(score_mat)

#3. Back-tracking

#finding the cell with the highest value

x = len_seq1
y = len_seq2

start_cell = score_mat[x][y]

for x in range(m):
    for y in range(n):
         if score_mat[x][y] > start_cell:
                start_cell = score_mat[x][y]

#print(start_cell)
#print(x,y)

aligned_seq1 = ""
aligned_seq2 = ""

while score_mat[x][y] != 0: #trace-back till the diagonal cell is 0

        nucleotide_1 = seq1[x-1]
        nucleotide_2  = seq2[y-1]

        if nucleotide_1 == nucleotide_2:
            aligned_seq1 = aligned_seq1 + nucleotide_1
            aligned_seq2 = aligned_seq2 + nucleotide_2
            x = x - 1
            y = y - 1

        else:
            top_cell = score_mat[x-1][y]
            left_cell = score_mat[x][y-1]
            di_cell = score_mat[x-1][y-1]

            great_2 = max(top_cell,left_cell,di_cell)

            if di_cell == great_2:
                aligned_seq1 = aligned_seq1 + nucleotide_1
                aligned_seq2 = aligned_seq2 + nucleotide_2
                x = x - 1
                y = y - 1
            elif top_cell == great_2:
                aligned_seq1 = aligned_seq1 + nucleotide_1
                aligned_seq2 = aligned_seq2 + "_"
                x = x - 1
            elif left_cell == great_2:
                aligned_seq1 = aligned_seq1 + "_"
                aligned_seq2 = aligned_seq2 + nucleotide_2
                y = y - 1

#reversing the aligned sequences:
aligned_seq1 = aligned_seq1[::-1]
aligned_seq2 = aligned_seq2[::-1]

#getting the alignment score
align_score = 0
align_sym = ''
for i in range(len(aligned_seq1)):
    if aligned_seq1[i] == aligned_seq2[i]:
        align_score = align_score + 1
        align_sym = align_sym + '|'
    elif aligned_seq1[i] == "-" or aligned_seq2[i] == "-":
        align_score = align_score - 1
        align_sym = align_sym + ' '
    elif aligned_seq1[i] != aligned_seq2[i]:
        align_score = align_score - 1
        align_sym = align_sym + "*"

print(aligned_seq1)
print(align_sym)
print(aligned_seq2)
print("Alignment Score: ",align_score)
