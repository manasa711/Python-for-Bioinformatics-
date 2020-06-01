#!/usr/bin/env python3

#Script for finding orthologs using reciprocal BLAST hits.

import argparse

import subprocess

def get_reciprocal_hits(file_one, file_two, input_sequence_type):

    op1 = []

    op2 = []

    output_list = []

    if input_sequence_type == 'n':

        #creating the nucleotide database for sequence A
        dbAn = "makeblastdb -in "+str(file_one)+" -dbtype nucl -out temp/dbAn"
        dbAn_subp = subprocess.check_output(dbAn.split())

        #performing blastn on sequence B with sequence A as the database
        bnA = "blastn -query "+str(file_two)+" -db temp/dbAn -out temp/blastn_AtoB.txt -max_target_seqs 1 -max_hsps 1 -outfmt 6"
        bnA_subp = subprocess.check_output(bnA.split())

        #creating the nucleotide database for sequence B
        dbBn = "makeblastdb -in "+str(file_two)+" -dbtype nucl -out temp/dbBn"
        dbBn_subp = subprocess.check_output(dbBn.split())

        #performing blastn on sequence A with sequence B as the database
        bnB = "blastn -query "+str(file_one)+" -db temp/dbBn -out temp/blastn_BtoA.txt -max_target_seqs 1 -max_hsps 1 -outfmt 6"
        bnB_subp = subprocess.check_output(bnB.split())


        with open("temp/blastn_AtoB.txt") as fhand:
            for line in fhand:
                a = line.split()[0] #gene B
                b = line.split()[1] #gene A
                op1.append(a+"\t"+b)

        #print("Length of output1:", len(op1))


        with open("temp/blastn_BtoA.txt") as fhand:
            for line in fhand:
                c = line.split()[0] #gene A
                d = line.split()[1] #gene B
                op2.append(d+"\t"+c)

        #print("Length of output2:", len(op2))

    elif input_sequence_type == 'p':

        #creating the protein database for sequence A
        dbAp = "makeblastdb -in "+str(file_one)+" -dbtype prot -out temp/dbAp"
        dbAp_subp = subprocess.check_output(dbAp.split())

        #performing blastp on sequence B with sequence A as the database
        bpA = "blastp -query "+str(file_two)+" -db temp/dbAp -out temp/blastp_AtoB.txt -max_target_seqs 1 -max_hsps 1 -outfmt 6"
        bpA_subp = subprocess.check_output(bpA.split())

        #creating the nucleotide database for sequence B
        dbBp = "makeblastdb -in "+str(file_two)+" -dbtype prot -out temp/dbBp"
        dbBp_subp = subprocess.check_output(dbBp.split())

        #performing blastn on sequence A with sequence B as the database
        bpB = "blastn -query "+str(file_one)+" -db temp/dbBp -out temp/blastp_BtoA.txt -max_target_seqs 1 -max_hsps 1 -outfmt 6"
        bpB_subp = subprocess.check_output(bpB.split())

        with open("temp/blastp_AtoB.txt") as fhand:
            for line in fhand:
                a = line.split()[0] #gene B
                b = line.split()[1] #gene A
                op1.append(str(a)+"\t"+str(b))

        with open("temp/blastp_BtoA.txt") as fhand:
            for line in fhand:
                c = line.split()[0] #gene A
                d = line.split()[1] #gene B
                op2.append(str(d)+"\t"+str(c))

    for i in op1:
        if i in op2:
            output_list.append(i)

    #print("Length of output list:", len(output_list))

    return(output_list)


def main():

# Argparse code
parser = argparse.ArgumentParser()
parser.add_argument("-i1", help= "This argument is for input file 1")
parser.add_argument("-i2" , help = "This argument is for input file 2")
parser.add_argument("-o", help="This argument takes in the name of the output file")
parser.add_argument("-t", help="This argument can be used to specify the type of sequence: nucleotide or protein")
args = parser.parse_args()

if args.i1:
    print("The first input file is:"+str(args.i1))

if args.i2:
    print("The second input file is:"+str(args.i2))

if args.o:
    print("The name of the output file specified:"+str(args.o))

if args.t:
    print("Type of sequence specified:"+str(args.t))

file_one = args.i1
file_two = args.i2
output_file = args.o
input_sequence_type = args.t


make_temp = "mkdir temp"
subprocess.call(make_temp.split())

output_list = get_reciprocal_hits(file_one, file_two, input_sequence_type)

print(len(output_list))

#print(output_list)

with open(output_file, 'w') as output_fh:
    for ortholog_pair in output_list:
        output_fh.write(ortholog_pair)
        output_fh.write("\n")

del_temp = "rm -r temp"
subprocess.call(del_temp.split())

if __name__ == "__main__":
    main()
