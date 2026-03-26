#Hazelyn Cates
#Started 3/13/26
#This is part of my Genomic Data Science Specialization class, for Module 2 final project
from re import match

import Bio
from Bio import SeqIO
import re #regex
from random import randrange

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
            longest_id.append(identifiers[i]) #if they match, add the corresponding identifier to the list
        elif i == min_index[j]: # same thing with the min index/indices
            shortest_id.append(identifiers[i]) #same with the min

print("\nThe identifier(s) of the longest sequence(s) is/are:")
print(longest_id)

print("\nThe identifier(s) of the shortest sequence(s) is/are:")
print(shortest_id)

print("\n%d sequences in file" % seq_nums)
print("\n...............................................................................................................")

print("\nReading frame 1:")
#Convert list of sequences in "seqs" to strings so I can parse through them
seqs_strings = list(map(str, seqs))
#print(len(seqs_strings))
#type(seqs_strings)


#For a given sequence identifier, what is the longest ORF contained in the sequence represented by that identifier?
#What is the starting position of the longest ORF in the sequence that contains it?
#Pick sequence identifier at random:
random = randrange(25) #pick a number between 0 and 25 at random
random_id = identifiers[random]
random_seq = 0 #default to zero, stores correct sequence index

for i in range(len(identifiers)): #loop through all identifiers
    #print(identifiers[i])
    if random_id == identifiers[i]: #find the one that matches
        random_seq = i #assign the proper sequence number to "random_seq"

#print(random_seq)

orf1_random = [] #store ORF1 of randomly picked sequence
orf2_random = [] #store ORF2 of randomly picked sequence
orf3_random = [] #stores ORF3 of randomly picked sequence

orf1, orf2, orf3 = [], [], []
orf_nums = 0 #all ORFs
#Storing SEQUENCE index (which sequence the ORF occurred in)
orf1_seq_index = []
orf2_seq_index = []
orf3_seq_index = []

#just capture length of ORFs of randomly picked sequence
random_orf_len = []

#Finding ORFs in first reading frame, starting at the beginning (index = 0)
for i in range(len(seqs_strings)): #goes through all 25 sequences
    seq = seqs_strings[i] #"seq" is going to act as a temporary hold through each iteration (redundant, but helps readability)
    #Where "seq" holds each of the 25 sequences one at a time. So for ex., seq[0] is sequence 1
    if re.findall(r"ATG.*TAA", seq): #first possible ORF
        #print(re.findall(r"ATG.*TAA", seq))
        orf1.append(re.findall(r"ATG.*TAA", seq)) #add it to a list containing this kind of ORF
        orf_nums += 1 #add up total number of ORFs in file (all sequences)
        orf1_seq_index.append(i) #which sequence ORF1 was found in
        #To do two things at once, find the ORFs in the sequence that was randomly picked and add it to it's corresponding ORF list
        if i == random_seq: #meaning, that's the sequence that was picked at random
            orf1_random.append(re.findall(r"ATG.*TAA", seq)) #add that sequence to a list containing ORF1 for that randomly picked sequence
        #print(seq[i[)
    if re.findall(r"ATG.*TAG", seq): #second possible ORF
        #print("ORF2 found ")
        orf2.append(re.findall(r"ATG.*TAG", seq))
        orf_nums += 1
        orf2_seq_index.append(i)
        if i == random_seq:
            orf2_random.append(re.findall(r"ATG.*TAG", seq))
        #print(seq[i])
    if re.findall(r"ATG.*TGA", seq): #third possible ORF
        #print("ORF3 found ")
        orf3.append(re.findall(r"ATG.*TGA", seq))
        orf_nums += 1
        orf3_seq_index.append(i)
        if i == random_seq:
            orf3_random.append(re.findall(r"ATG.*TGA", seq))
        #print(seq[i])
    else:
        print("No ORFs found ")

#print(orf1_random)
#print(orf2_random)
#print(orf3_random)
max1 = []
max2 = []
max3 = []

#print(random_seq)

for i in range(len(orf1_random)):
    if not orf1_random: #check if list is empty
        print("Sequence %d does not have any ORF1" % random_seq)
    else: #if it's not
        for j in range(len(orf1_random)):
            max1.append(len(orf1_random[i][j]))

for i in range(len(orf2_random)):
    if not orf2_random:
        print("Sequence %d does not have any ORF2" % random_seq)
    else:
        for j in range(len(orf2_random)):
            max2.append(len(orf2_random[i][j]))

for i in range(len(orf3_random)):
    if not orf3_random:
        print("Sequence %d does not haev any ORF3" % random_seq)
    else:
        for j in range(len(orf3_random)):
            max3.append(len(orf3_random[i][j]))

#print(max1)
#print(max2)
#print(max3)

max_rand_len = max(max1 + max2 + max3)

print("Longest ORF in sequence %d is %d base pairs long\n" % (random_seq, max_rand_len))

start_pos1 = [] #stores ORF1 start positions
end_pos1 = [] #stores ORF1 end positions

start_pos2 = []
end_pos2 = []

start_pos3 = []
end_pos3 = []

#These store the positions of the ORFs in the randomly selected sequence
rand_start = []
rand_end = []

#I'll get to that next
for i in range(len(seqs_strings)):
    if i == random_seq:
        seq = seqs_strings[i]
        for m in re.finditer(r"ATG.*TAA", seq):  # ORF1
            rand_start1 = m.start() #Get start position
            rand_end1 = m.end() #Get end position
            print("ORF1 sequence %d:" % random_seq)
            print("Start:", rand_start1)
            print("End:", rand_end1)
            print(rand_end1 - rand_start1, "\n")

        for m in re.finditer(r"ATG.*TAG", seq): #ORF2
            rand_start2 = m.start()
            rand_end2 = m.end()
            print("ORF2 sequence %d:" % random_seq)
            print("Start:", rand_start2)
            print("End:", rand_end2)
            print(rand_end2-rand_start2, "\n")

        for m in re.finditer("ATG.*TGA", seq): #ORF 3
            rand_start3 = m.start()
            rand_end3 = m.end()
            print("ORF3 sequence %d:" % random_seq)
            print("Start:", rand_start3)
            print("End:", rand_end3)
            print(rand_end3-rand_start3, "\n")
    else:
        continue


print("\nFor reading frame 1, with sequences starting at index 0:")
print("%d ORFs found within all %d FASTA sequences" % (orf_nums, len(seqs_strings)))
#print("\n", orf1, "\n", orf2, "\n", orf3)
#print("\nSequences ORF1 found in: ", orf1_seq_index, "\nSequences ORF2 found in: ", orf2_seq_index, "\nSequences ORF3 found in: ", orf3_seq_index)

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
print("\nLongest ORF1 in reading frame 1 is %d bases long in sequence %d" % (max_orf1, max_ORF_index1))

#ORF 2
for i in range(len(orf2)):
    for j in range(len(orf2[i])):
        orf2_lengths.append(len(orf2[i][j]))
    max_orf2 = max(orf2_lengths)
    if i == orf2_seq_index[i]:
        max_ORF_index2 = orf2_lengths.index(max_orf2)
    #print(max_ORF_index2)
print("Longest ORF2 in reading frame 1 is %d bases long in sequence %d" % (max_orf2, max_ORF_index2))

#ORF 3
for i in range(len(orf3)):
    for j in range(len(orf3[i])):
        orf3_lengths.append(len(orf3[i][j]))
    max_orf3 = max(orf3_lengths)
    if i == orf3_seq_index[i]:
        max_ORF_index3 = orf3_lengths.index(max_orf3)
    #print(max_ORF_index3)
print("Longest ORF3 in reading frame 1 is %d bases long in sequence %d" % (max_orf3, max_ORF_index3))


#Find longest ORF of all of them:
longest_orf = max(max_orf1, max_orf2, max_orf3)
#Check to see which ORF group has the longest overall ORF and assign correct index
if longest_orf == max_orf1:
    longest_orf_index = max_ORF_index1
if longest_orf == max_orf2:
    longest_orf_index = max_ORF_index2
if longest_orf == max_orf3:
    longest_orf_index = max_ORF_index3

print("\nLongest ORF in file with reading frame 1 is %d base pairs long in sequence %d" % (longest_orf, longest_orf_index))
print("Sequence identifier of sequence %d:" % (longest_orf_index))
print(identifiers[longest_orf_index])

#FINSIH THIS!!!!!!!!!!!!
#FILL IN WITH OTHER READING FRAMES
print("\n\n...............................................................................................................")
print("\nReading frame 2 (+1)")
i = 1
#OK so for shifting the reading frame, it only matters in the loops where I'm searching for the ORFs and their positions
#All I have to do is set i = 1 for this block and change it to i=2 for the third, and I can just copy and paste the code






print("\n\n...............................................................................................................")
print("\nReading frame 3 (+2)")
i = 2



#4.) Find all repeats of length n
#Use regex (most likely)
#SO any sequence of length n, NOT a SPECIFIC repeat
#I know these FASTA files are correct (i.e. no other bases besides the usual 4)
#So it would be something like "/w{n}" for all alphanumeric values n times
#So like, if you want to find repeats where n = 3
#Say you have the sequence: ATTTAAGAGAGAGTTTTATTTATTTATTTTG
#TTT is one, GAG is another, ATT is another, etc.
#HAVE to make sure it allows for OVERLAPPING repeats!!!!

