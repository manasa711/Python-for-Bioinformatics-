#!/usr/bin/env python3

import sys

#storing the arguments into variables

seq1_file = sys.argv[1]

seq2_file = sys.argv[2]

#Retrieving sequences from the fasta files
seq1 = ""

seq2 = ""

with open(seq1_file) as fhand1:
    for line in fhand1:
        line.rstrip()
        if line.startswith(">"):
            continue
        seq1 = seq1 + line.strip("\n")

with open(seq2_file) as fhand2:
    for line in fhand2:
        line.rstrip()
        if line.startswith(">"):
            continue
        seq2 = seq2 + line.strip("\n")

#print(seq1)
#print(seq2)

seq1 = list(seq1)
seq2 = list(seq2)

len_seq1 = len(seq1) #getting the length of sequence 1
len_seq2 = len(seq2) #getting the lenght of sequence 2

m = len_seq1+1
n = len_seq2+1

#1. initializing a matrix using lists
#and filling it up with 0s
score_mat = [[0 for i in range(n)] for i in range(m)]

#print_mat(score_mat)

#assigning scores
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

#Function for printing out the matrix
def print_mat(matrix_name):
    for row in matrix_name:
        for i in row:
            print(i, end = "\t")
        print("\n")

#filling the first row with gap penalties
for i in range(m):
    score_mat[i][0] = gap * i

#filling the first column with gap penalties
for j in range(n):
    score_mat[0][j] = gap * j

#score_mat [0][0] = 0 #making the first cell zero

#2. Matrix Filling
#filling up the cells based on matches, mismatches and gaps
for i in range(1,m):
    for j in range(1,n):
        diagonal_score = score_mat[i-1][j-1] + match_vs_mismatch(seq1[i-1],seq2[j-1])
        leftcell_score = score_mat[i][j-1] + gap
        topcell_score = score_mat[i-1][j] + gap
        score_mat[i][j] = max(diagonal_score,leftcell_score,topcell_score)

#printing the matrix
#print_mat(score_mat)

#3. Backtracking Step

x = len(seq1)
y = len(seq2)

aligned_seq1 = ""
aligned_seq2 = ""

while x>0 and y>0: #starting with the bottomcell

    nucleotide_1 = seq1[x-1]
    nucleotide_2  = seq2[y-1]

    current_cell = score_mat[x][y] #storing the value of the current cell
    #print(current_cell)

    if nucleotide_1 == nucleotide_2: #if its a match, the path moves diagonally
        aligned_seq1 = aligned_seq1 + nucleotide_1
        aligned_seq2 = aligned_seq2 + nucleotide_2
        x = x - 1
        y = y - 1

    else: #if its a mismatch
        top_cell = score_mat[x-1][y] #storing the value of the top cell
        left_cell = score_mat[x][y-1] #storing the value of the left cell
        di_cell = score_mat[x-1][y-1] #storing the value of the diagonal cell

        great = max(top_cell,left_cell,di_cell) #finding the maximum value

        if di_cell == great: #if the diagonal cell value is greater, the path moves there
            aligned_seq1 = aligned_seq1 + nucleotide_1
            aligned_seq2 = aligned_seq2 + nucleotide_2
            x = x - 1
            y = y - 1
        elif top_cell == great: #if the top cell value is greater, the path moves up
            aligned_seq1 = aligned_seq1 + nucleotide_1
            aligned_seq2 = aligned_seq2 + "-"
            x = x - 1
        elif left_cell == great: #if the left cell value is greater, the path moves to the left
            aligned_seq1 = aligned_seq1 + "-"
            aligned_seq2 = aligned_seq2 + nucleotide_2
            y = y - 1

#reversing the aligned sequences:
aligned_seq1 = aligned_seq1[::-1]
aligned_seq2 = aligned_seq2[::-1]

#Alignment Score
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
