#Hazelyn Cates
#Started 3/13/26
#This is part of my Genomic Data Science Specialization class, for Module 2 final project
from re import match

import Bio
from Bio import SeqIO
import re #regex

#for storing number of sequences in file:
seq_nums = 0

#Read in FASTA file
#Get name of file:
print("Enter FASTA file name including .fasta file extension:")
fasta_file = input()

sample = SeqIO.parse(fasta_file, "fasta")

seqs = [] #empty list to store each sequence
identifiers = [] #stores identifiers of each sequence
lengths = [] #stores lengths of each sequence
#max_seq_lens = [] #list to hold max sequence length, duplicates possible
max_index = [] #list to hold index/indices of longest sequence(s)
#min_seq_lens = [] #list to hold min sequence length, duplicates possible
min_index = [] #list to hold index/indices of shortest sequence(s)

#Parse through FASTA file:
for record in sample: #where "record" is each entry in the file
    #print(record.seq)
    seqs.append(record.seq) #add each sequence to the list
    seq_nums += 1 #count how many sequences are in the file
    identifiers.append(record.id) #get the identifier of each sequence
    #print(record.id)

    lengths.append(len(record.seq)) #get the length of each sequence
    #print(len(record.seq))

print("The length of each sequence are as follows:")
print(lengths)

#Find the longest sequence
longest_length = max(lengths)

longest_index = lengths.index(longest_length)
#max_index.append(longest_index)
#print(max_index)

#Find the shortest sequence
shortest_length = min(lengths)

shortest_index = lengths.index(shortest_length)
#min_index.append(shortest_index)
#print(min_index)

#check for max or min seq length DUPLICATES
for i in range(seq_nums): #iterate for as many seqs are in the file
    if lengths[i] == longest_length: #if the value matches, add its index
        max_index.append(i) #add the index in which the duplicate was found
    elif lengths[i] == shortest_length: #and same for min
        min_index.append(i)
    else:
        continue

print("\nThe longest sequence is %d base pairs long" % longest_length)
print("The shortest sequence is %d base pairs long" % shortest_length)
#print(identifiers)
#print(max_index)
#print(min_index)

#Store the identifiers of the longest and shortest sequences
longest_id = []
shortest_id = []

#NOW,
#Want to match the max(s) and min(s) to the correct sequence identifier
#If there were multiple, would have to use a for loop and iterate the length of the max and min indices
#Or might just do that anyway under the assumption there are multiple
#Since this program has to work with any FASTA file

for i in range(len(identifiers)): #go through all 25 identifiers
    for j in range(len(max_index)): #go through however many max lengths there are
        if i == max_index[j]: #where index 0 = whatever value, index 1 = other value, etc.
            longest_id.append(identifiers[i]) #if they match, add the corresponding identfier to the list
        elif i == min_index[j]: # same thing with the min index/indices
            shortest_id.append(identifiers[i]) #same with the min

print("\nThe identifier(s) of the longest sequence(s) is/are:")
print(longest_id)

print("\nThe identifier(s) of the shortest sequence(s) is/are:")
print(shortest_id)

print("\n%d sequences in file" % seq_nums)

#Convert list of sequences in "seqs" to strings so I can parse through them
seqs_strings = list(map(str, seqs))
#print(len(seqs_strings))
#type(seqs_strings)

#Now the get the open reading frames, I think I'm gonna have to do another nested for loop
#I need to go sequence by sequence in the "seq" list and look for and "atg" to a "taa", "tga", or "tag"
#Can maybe use regex????

rows, cols = (25, longest_length) #want to make sure it's big enough to hold the longest sequence
#Create 2D array to store the ORFs for each of the 25 sequences
#orfs = [[0 for i in range(cols)] for j in range(rows)] #stores ORFs from each sequence

orf1, orf2, orf3 = [], [], []
orf_nums = 0 #all ORFs
#Storing SEQUENCE index (which sequence the ORF occured in)
orf1_seq_index = []
orf2_seq_index = []
orf3_seq_index = []

"""#Storing POSITION index
position_orfs1 = []
position_orfs2 = []
position_orfs3 = []"""


#First reading frame, starting at the beginning (index = 0)
for i in range(len(seqs_strings)): #goes through all 25 sequences
    seq = seqs_strings[i] #"seq" is going to act as a temporary hold through each iteration (redundant, but helps readability)
    #Where "seq" holds each of the 25 sequences one at a time. So for ex., seq[0] is sequence 1
    if re.findall(r"ATG.*TAA", seq): #first possible ORF
        #print(re.findall(r"ATG.*TAA", seq))
        orf1.append(re.findall(r"ATG.*TAA", seq)) #add it to a list containing this kind of ORF
        orf_nums += 1 #add up total number of ORFs in file (all sequences)
        orf1_seq_index.append(i) #which sequence ORF1 was found in
        #print(seq[i[)
    if re.findall(r"ATG.*TAG", seq): #second possible ORF
        #print("ORF2 found ")
        orf2.append(re.findall(r"ATG.*TAG", seq))
        orf_nums += 1
        orf2_seq_index.append(i)
        #print(seq[i])
    if re.findall(r"ATG.*TGA", seq): #third possible ORF
        #print("ORF3 found ")
        orf3.append(re.findall(r"ATG.*TGA", seq))
        orf_nums += 1
        orf3_seq_index.append(i)
        #print(seq[i])
    else:
        print("No ORFs found ")


#Now want to get the start and end positions of each ORF in each sequence for reading frame 1
#NOTE this was a test (which worked!!) but the start and end positions aren't stored in any variables
#I'll get to that next
for i in range(len(seqs_strings)):
    seq = seqs_strings[i]
    for m in re.finditer(r"ATG.*TAA", seq): #ORF1
        start1 = m.start() #gets start position of ORF1 in given sequence (starting with sequence 1)
        end1 = m.end() #gets end position of ORF1 in a given sequence (starting with sequence 1)
        #print(start1, "\n", end1)
    for m in re.finditer(r"ATG.*TAG", seq):
        start2 = m.start()
        end2 = m.end()
        #print(start2, "\n", end2)
    for m in re.finditer("ATG.*TGA", seq):
        start3 = m.start()
        end3 = m.end()
        #print(start3, "\n", end3)


#print(orf1_seq_index)

#For a given sequence identifier, what is the longest ORF contained in the sequence represented by that identifier??



seq_orfs = [[0 for i in range(3)] for j in range(len(seqs_strings))]
#this variable creates a list of lists with 3 rows (one for each ORF) and 25 columns (one for each sequence)
#The idea is to store the ORFs from each sequence IN the index for that sequence if the ORF was FOUND in that sequence
#So like, seq_orfs[0][0] would hold the first ORF for sequence 1
#And seq_orfs[[1][0] would hold the second ORF for sequence 1
#This is me attempting to match each ORF to its sequence in a 3 row x 25 column list of lists
#It's giving me an error at the moment, but I'm working on it
#*NOTE this is for READING FRAME 1!

#WAIT A MINUTE I know why this is a problem and giving me like, a list of lists of lists
#The variable "orf1" is a list ALREADY. So when I add an index of it to "seq_orfs", it's adding a list to a list,
#And "seq_orfs" was created as a list of lists.
#So then, it becomes a list of lists of lists
#Shit
#OK IMPORTANT!!! This is JUST for ORF 1!!!!!!!
for i in range(len(orf1)):
    #print(i)
    for j in range(0, 3):
        if i in orf1_seq_index: #means there's an ORF in that sequence
            #print(orf1[i])
            seq_orfs[j][0] = orf1[i] #put the corresponding ORF in THAT SEQUENCE into the row j at column i
            #seq_orfs.append(re.findall(r"ATG.*TAA", seq))
            #print(seq_orfs)
        elif i not in orf1_seq_index:
            continue

print("\nFor reading frame 1, with sequences starting at index 0:")
print("%d ORFs found within all %d FASTA sequences" % (orf_nums, len(seqs_strings)))
#print("\n", orf1, "\n", orf2, "\n", orf3)
print("\nSequences ORF1 found in: ", orf1_seq_index, "\nSequences ORF2 found in: ", orf2_seq_index, "\nSequences ORF3 found in: ", orf3_seq_index)

#SO now I know which sequences have which ORFs
#I need to find the POSITIONS of each


#This checks the lengths of the ORF LISTS, not the ORFs themsevles (should be max of 25 (for 25 sequences))
#print("\n", len(orf1), "\n", len(orf2), "\n", len(orf3))
#NOW to find the length of each ORF in each ORF list:
orf1_lengths = []
orf2_lengths = []
orf3_lengths = []

#This is horribly inefficient but I cannot come up with a better way at the moment
#ORF 1
for i in range(len(orf1)): #23 times
    for j in range(len(orf1[i])): #iterates through the length of ORFs in the list of lists (orf1 in this loop)
        orf1_lengths.append(len(orf1[i][j])) #finds the length of each ORF and adds it to the list
    max_orf1 = max(orf1_lengths) #gets max ORF length as it goes through each one, updates accordingly
    #Want to go by VALUE in "orf1_seq_index" list to make sure the INDEX is correct for the longest ORF (not by index alone)
    if i == orf1_seq_index[i]: #this ensures that the correct index is assigned
        max_ORF_index1 = orf1_lengths.index(max_orf1) #gets index of longest ORF
    #print(max_ORF_index1)
print("\nLongest ORF1 in reading frame 1 is %d bases long in sequence %d" % (max_orf1, max_ORF_index1 + 1))

#ORF 2
for i in range(len(orf2)):
    for j in range(len(orf2[i])):
        orf2_lengths.append(len(orf2[i][j]))
    max_orf2 = max(orf2_lengths)
    if i == orf2_seq_index[i]:
        max_ORF_index2 = orf2_lengths.index(max_orf2)
    #print(max_ORF_index2)
print("Longest ORF2 in reading frame 1 is %d bases long in sequence %d" % (max_orf2, max_ORF_index2 + 1))

#ORF 3
for i in range(len(orf3)):
    for j in range(len(orf3[i])):
        orf3_lengths.append(len(orf3[i][j]))
    max_orf3 = max(orf3_lengths)
    if i == orf3_seq_index[i]:
        max_ORF_index3 = orf3_lengths.index(max_orf3)
    #print(max_ORF_index3)
print("Longest ORF3 in reading frame 1 is %d bases long in sequence %d" % (max_orf3, max_ORF_index3 + 1))

"""print("\nORF1 lengths:", orf1_lengths)
print("\nORF2 lengths:", orf2_lengths)
print("\nORF3 lengths:", orf3_lengths)"""

#Find longest ORF of all of them:
longest_orf = max(max_orf1, max_orf2, max_orf3)
#Check to see which ORF group has the longest overall ORF and assign correct index
if longest_orf == max_orf1:
    longest_orf_index = max_ORF_index1
if longest_orf == max_orf2:
    longest_orf_index = max_ORF_index2
if longest_orf == max_orf3:
    longest_orf_index = max_ORF_index3

print("\nLongest ORF in file with reading frame 1 is %d base pairs long in sequence %d" % (longest_orf, longest_orf_index +1))
print("Sequence identifier of sequence %d:" % (longest_orf_index+1))
print(identifiers[longest_orf_index])

#FINSIH THIS!!!!!!!!!!!!
#FILL IN WITH OTHER READING FRAMES
print("...............................................................................................................")



#NOW need to find longest ORF in EACH SEQUENCE
#This is the hard part b/c I need to associate each ORF to the proper sequence
#But there are a lot of sequences so I CANNOT do the brute-force method

#ALSO need to get start and end positions for the longest ORF in each sequence.
#I think that's going to be a bitch. Might have to start again and use the "finditer" function to actually get
#indicies of matches. We'll see. (AHHHHHHHHHHHHHH!!!!) Ok breathe. Cool. Yeah. Let's do it.

#Could go by SEQUENCE index and compare length values from each one
#So like, sequence 1 (index 0), you could compare all ORFs in that sequence
#I know that sequences 7 and 8 only have ORF2's and ORF3's, so just 2 to compare


#***NOTE:
#COULD POTENTIALLY USE A "WHILE LOOP" TO FIND THE STARTING POSITIONS OF THE ORFS u
#OR TO FIND THE ORFS AND GO UNTIL WE REACH THE END OF THE SEQUENCE (position > -1)
#Use the "find" function WITHIN the WHILE loop to get ALL the positions!!!
#See 4.2 lecture
