#Hazelyn Cates
#Started 3/13/26
#This is part of my Genomic Data Science Specialization class, for Module 2 final project
#from re import match

#import Bio
from Bio import SeqIO
from collections import Counter
#import re #regex
#from random import randrange

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
longest_lengths = [] #stores lengths of sequences that are the longest (if there are multiple)
shortest_lengths = [] #stores lengths of sequences that are the shortest (if there are multiple)
longest_indices = [] #stores indices of longest sequences (if there are multiple)
shortest_indices = [] #stores indices of shortest sequences (if there are multiple)

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

longest_length = max(lengths) #calculates longest sequence in all sequences in file
shortest_length = min(lengths) #calculates shortest sequence in all sequences in file
#Have to see if there is more than 1 sequence w/ the longest length (i.e. they're the same)
#Find the longest sequence
for i in range(len(lengths)):
    if lengths[i] == longest_length: #this will check the length of all sequences, including "longest_length" itself
        longest_lengths.append(lengths[i]) #all longest lengths will be added
        longest_indices.append(lengths.index(lengths[i])) #along with there indices
    elif lengths[i] == shortest_length: #checks shortest sequences, including "shortest_length" itself
        shortest_lengths.append(lengths[i])
        shortest_indices.append(lengths.index(lengths[i]))
#print(longest_lengths)


print("\nLongest sequence length(s):")
print(longest_lengths)
print("Longest sequence index/indices:")
print(longest_indices)
print("\nShortest sequence length(s):")
print(shortest_lengths)
print("Shortest sequence index/indices:")
print(shortest_indices)

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
    for j in range(len(longest_indices)): #go through however many longest lengths there are
        if i == longest_indices[j]: #where index 0 = whatever value, index 1 = other value, etc.
            longest_id.append(identifiers[i]) #if they match, add the corresponding identifier to the list
    for h in range(len(shortest_indices)): #go through however many shortest lengths there are
        if i == shortest_indices[h]: # same thing with the min index/indices
            shortest_id.append(identifiers[i]) #same with the min

print("\nThe identifier(s) of the longest sequence(s) is/are:")
print(longest_id)

print("\nThe identifier(s) of the shortest sequence(s) is/are:")
print(shortest_id)

print("\n%d sequences in file" % seq_nums)
print("\n\n")

#Convert list of sequences in "seqs" to strings so I can parse through them
seqs_strings = list(map(str, seqs))
#print(len(seqs_strings))
#type(seqs_strings)


#Create function that finds ORFs given the input frame
def get_orfs(seq, frame):
    print("\nThis is READING FRAME %d" % frame)
    stop_codons = ["TAA", "TGA", "TAG"] #three possible stop codons to terminate ORF
    orfs1 = [] #stores ORFs
    lengths = [] #stores lengths of ORFs
    j = frame-1 #start at either 0, +1, or +2

    for i in range(len(seqs_strings)): #go through all 25 sequences
        seq = seqs_strings[i] #holder for the current sequence, helps readability

        while j <= len(seq) - 3: #go through sequence making sure at least one codon is left in the current sequence
            #print(i)
            codon = seq[j:j+3] #iterate through bases in sets of three, which is a codon
            if codon == "ATG": #check to see if current codon is a start codon
                h = j + 3
                while h <= len(seq) - 3: #meaning go through current sequence
                    stop = seq[h: h+3] #keeps track of each codon
                    if stop in stop_codons: #if "stop" hits one of the three stop codons:
                        start1 = j + 1
                        end1 = h + 3
                        length1 = end1 - start1
                        orfs1.append((start1, end1, length1, i))
                        #print(i)
                        #print(orfs)
                        break #break out inner loop, find more in current sequence
                    h += 3
                j = h if h > j else j + 3
            else:
                j += 3

    length_orfs1 = []
    orf_indices1 = []

    #print(orfs)

    for i in range(len(orfs1)):
        length_orfs1.append(orfs1[i][2])
        orf_indices1.append(orfs1[i][3])

    print("\nReading frame %d ORF lengths:" % frame)
    print(length_orfs1)
    print("\nSequences of reading frame %d ORFs:" % frame)
    print(orf_indices1)


    #Now need to get max ORF in the file and it's identifier
    max_orf1 = max(length_orfs1)
    for i in range(len(orf_indices1)):
        if max_orf1 == length_orfs1[i]:
            print("\nSequence %d has the longest ORF at %d base pairs long" % (orf_indices1[i]+1, max_orf1))
            print("Sequence %d identifier:" % (orf_indices1[i]+1))
            print(identifiers[orf_indices1[i]])

    #For a given sequence identifier, what is the longest ORF contained in the sequence represented by that identifier?
    #What is the starting position of the longest ORF in the sequence that contains it?

    #this is for like if the user actually puts in the entire string for the identifier name, not just a sequence number
    """user_identifier = string(input())

    for i in range(len(identifiers)):
        if user_identifier == identifiers[i]:
            picked_identifier = i"""

    #Let user pick the identifier
    print("\nPick a number between 0 and %d" % (seq_nums-1))
    picked_identifier = int(input())
    #print(picked_identifier)
    while picked_identifier < 0 or picked_identifier > (seq_nums-1):
        print("Pick a number between 0 and %d" % (seq_nums-1))
        picked_identifier = int(input())

    if picked_identifier in orf_indices1:
        print("Sequence identifier of sequence %d" % (picked_identifier+1))
        print(identifiers[picked_identifier])
    else:
        print("Pick a number from the following list:")
        print(orf_indices1)
        picked_identifier = int(input())
        while picked_identifier not in orf_indices1:
            print("Pick a number from the following list:")
            print(orf_indices1)
            picked_identifier = int(input())
        print("Sequence identifier of sequence %d" % (picked_identifier + 1))
        print(identifiers[picked_identifier])


    #With the identifier chosen and valid, have to find longest ORF in that sequence
    #So if sequence 1 (or 0 in this case) was chosen, have to find the ORFs that have their fourth item in the list match
    #The chosen identifier
    #So like:
    picked_orf1_len = []
    picked_orf1_start = 0 #stores start position of chosen sequence ORF
    picked_orf1_end = 0 #stores end position of chosen sequence ORF
    for i in range(len(orfs1)):
        if orfs1[i][3] == picked_identifier: #where orfs1[i][3] is the sequence index where that ORF was found (1-25)
            picked_orf1_len.append(orfs1[i][2]) #add lengths of correct sequence ORFs to a list
            if max(picked_orf1_len) == orfs1[i][2]:  # find the entry that has the matching length
                picked_orf1_start = orfs1[i][0] #store the start position
                picked_orf1_end = orfs1[i][1] #store the end position

    #THIS IS CORRECT
    #BUT
    #I'm thinking I have to add 1 to the start and end positions since indexing starts at the 0th position
    #So I'm pretty sure they're a base behind
    return(print("\nSequence %d longest ORF is %d base pairs long starting at position %d and ending at position %d"
          % ((picked_identifier+1), max(picked_orf1_len), (picked_orf1_start+1), (picked_orf1_end+1) )))


#call ORF finder function, ask for user input
print("What reading frame do you want? 1, 2, or 3")
frame = int(input())

get_orfs(seqs_strings ,frame)

print("Enter repeat length:")
repeat_len = int(input())


#Write function that finds ALL repeats (overlaps allowed) of length n in each sequence in the FASTA file
#Given a length n, your program should be able to identify all repeats of length n in all sequences in the FASTA file.
#your program should also determine how many times each repeat occurs in the file,
#and which is the most frequent repeat of a given length.
def get_repeats(seqs, n):
    repeats = Counter() #stores repeats of each variety, recording how many times that repeat occurs in all sequences

    for i in seqs: #go through all sequences in file
        if n <= 0 or n > len(seqs): #check that n (repeat length) is valid
            continue
        for j in range(0, (len(seqs) - n + 1)): #go through that current sequence not exceeding the length of the repeat
            repeat = i[j:j+n] #Where "i" is the current sequence and [j:j+n] is the repeat ("ATG, "GGTC", whatever)
            repeats[repeat] += 1 #add the repeat to the "repeats" dictionary along with its count

    max_repeat = max(repeats, key = repeats.get)
    print("Most frequent repeat is:")
    print(max_repeat)
    return(repeats, max_repeat)

get_repeats(seqs_strings, repeat_len)
