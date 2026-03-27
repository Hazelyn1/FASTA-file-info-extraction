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

#Convert list of sequences in "seqs" to strings so I can parse through them
seqs_strings = list(map(str, seqs))
#print(len(seqs_strings))
#type(seqs_strings)


print("\n...............................................................................................................")
#This is REALLY going to be finding the ORFs for reading frame 1 :(
print("Reading frame 1:")
stop_codons = ["TAA", "TGA", "TAG"] #three possible stop codons to terminate ORF
orfs1 = [] #stores ORFs
lengths = [] #stores lengths of ORFs
j = 0 #start at first index (reading frame 1)

for i in range(len(seqs_strings)): #go through all 25 sequences
    seq = seqs_strings[i] #holder for the current sequence, helps readability
    #print(i)
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

print("\nReading frame 1 ORF lengths:")
print(length_orfs1)
print("\nSequences of reading frame 1 ORFs:")
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

#Let user pick the identifier
print("Pick a number between 0 and %d" % (seq_nums-1))
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
    while (picked_identifier not in orf_indices1):
        print("Pick a number from the following list:")
        print(orf_indices1)
        picked_identifier = int(input())
    print("Sequence identifier of sequence %d" % (picked_identifier + 1))
    print(identifiers[picked_identifier])


#With the identifier chosen and valid, have to find longest ORF in that sequence
#So if sequence 1 (or 0 in this case) was chosen, have to find the ORFs that have their fourth item in the list match
#The chosen identifier
#So like:
picked_orf_len = []
for i in range(len(orfs1)):
    if orfs1[i][3] == picked_identifier: #where orfs1[i][3] is the sequence where that ORF was found (1-25)
        picked_orf_len.append(orfs1[i][2]) #add lengths of correct sequence ORFs to a list

print("Sequence %d longest ORF is %d base pairs long" % ((picked_identifier+1), max(picked_orf_len)))

#Now have to find the start and end positions of THAT LONGEST ORF
#Could approach it the same way using the "max" function and matching it up with the proper entries in "orf1"
#In this case, have to find the entry of the correct sequence with the longest length entry (third position, index 2)
#and then extract the start and end positions of that SAME entry (indices 0 and 1, respectively)
#But that's a problem for tomorrow b/c I'm tired rn



print("\n\n\n................................................................................................................................")
print("Reading frame 2:")
#start j = 1
j = 1
orfs2 = []

for i in range(len(seqs_strings)): #go through all 25 sequences
    seq = seqs_strings[i] #holder for the current sequence, helps readability
    while j <= len(seq) - 3: #go through sequence making sure at least one codon is left in the current sequence
        codon = seq[j:j+3] #iterate through bases in sets of three, which is a codon
        if codon == "ATG": #check to see if current codon is a start codon
            h = j + 3
            while h <= len(seq) - 3: #meaning go through current sequence
                stop = seq[h: h+3] #keeps track of each codon
                if stop in stop_codons: #if "stop" hits one of the three stop codons:
                    start2 = j + 1
                    end2 = h + 3
                    length2 = end2 - start2
                    orfs2.append((start2, end2, length2, i))
                    #print(i)
                    #print(orfs)
                    break #break out inner loop, find more in current sequence
                h += 3
            j = h if h > j else j + 3
        else:
            j += 3

length_orfs2 = []
orf_indices2 = []
#print(len(orfs2))

for i in range(len(orfs2)):
    length_orfs2.append(orfs2[i][2])
    orf_indices2.append(orfs2[i][3])

#FINISH THIS!!!!!!!
print("Reading frame 2 ORF lengths:")
print(length_orfs2)
print("Sequence of reading frame 2 ORFs:")
print(orf_indices2)



print("\n\n\n................................................................................................................................")

#FINISH THIS!!!!!
print("Reading frame 3:")
j = 2
orfs3 = []
for i in range(len(seqs_strings)): #go through all 25 sequences
    seq = seqs_strings[i] #holder for the current sequence, helps readability
    while j <= len(seq) - 3: #go through sequence making sure at least one codon is left in the current sequence
        codon = seq[j:j+3] #iterate through bases in sets of three, which is a codon
        if codon == "ATG": #check to see if current codon is a start codon
            h = j + 3
            while h <= len(seq) - 3: #meaning go through current sequence
                stop = seq[h: h+3] #keeps track of each codon
                if stop in stop_codons: #if "stop" hits one of the three stop codons:
                    start3 = j + 1
                    end3 = h + 3
                    length3 = end3 - start3
                    orfs3.append((start3, end3, length3, i))
                    #print(i)
                    #print(orfs)
                    break #break out inner loop, find more in current sequence
                h += 3
            j = h if h > j else j + 3
        else:
            j += 3

length_orfs3 = []
orf_indices3 = []
#print(orfs)

for i in range(len(orfs3)):
    length_orfs3.append(orfs3[i][2])
    orf_indices3.append(orfs3[i][3])

print("Reading frame 3 ORF lengths:")
print(length_orfs3)
print("Sequence of reading frame 3 ORFs:")
print(orf_indices3)




















print("\n\n\n\n\n\n\n\n\n...............................................................................................................")

print("\nReading frame 1:")

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

#Finding ORFs in first reading frame, starting at the beginning (index i = 0)
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

#Store max length values of the ORFs in the sequence that was randomly chosen
max1 = []
max2 = []
max3 = []

#print(random_seq)
#This index i will NOT change with the reading frame
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
        print("Sequence %d does not have any ORF3" % random_seq)
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

#Now find the start and end positions of the ORFs
#FIX THIS!! IT ISN'T CORRECT!! IT'S NOT GOING BY TRIPLETS, IT'S JUST PATTERN MATCHING!!!!!
#DAMN IT!!! I spent so much time on this WAAAAAAAHHHHHHHHHHHHHHHHH!!!!!
for i in range(len(seqs_strings)): #goes through all 25 sequences
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

#OK so for shifting the reading frame, it only matters in the loops where I'm searching for the ORFs and their positions
#All I have to do is skip ahead one in the sequence being parsed (using [1:], where it starts at position 1, not 0

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

#Finding ORFs in second reading frame, starting at the 1 base forward ([1:])
for i in range(len(seqs_strings)): #goes through all 25 sequences
    seq = seqs_strings[i] #"seq" is going to act as a temporary hold through each iteration (redundant, but helps readability)
    #Where "seq" holds each of the 25 sequences one at a time. So for ex., seq[0] is sequence 1

    #FIX THIS!!!! I think I have this wrong. Starting at one position forward does NOT change how the ORFs are detected
    #Meaning, it's not putting the reading frame off by 1
    #ORF start and end codons are THREE bases in length
    #A shift in the reading frame would DISPLACE these codons, either creating "new" ones and disrupted old ones (i.e. in reading frame 1)
    #So I don't know how to mimic this displacement using regex
    #Might have to use a more brute force method???
    #Like what's described here???? --> https://python.plainenglish.io/fasta-essentials-with-python-records-lengths-orfs-and-repeats-fcf6f3ff8c31

    #OK WAIT I SEE WHAT THE PROBLEM IS AND WHY THIS ISN'T WORKING!!!!
    #EVEN THOUGH I'M STARTING ONE BASE PAIR FORWARD, IT DOESN'T SHIFT THE READING FRAME
    #CODONS ARE READ IN TRIPLETS
    #MY USE OF REGEX DOESN'T CARE ABOUT TRIPLETS IT JUST FINDS THOSE ORFS that match the pattern!!!
    #Which works when the reading frame is 1 (i.e. starts at index 0)
    #But NOT for when the reading frame shifts BECAUSE it's still just pattern matching
    #SHOOT!!!!!!!!!!!!!!!!
    if re.findall(r"ATG.*TAA", seq[1:]): #first possible ORF STARTING one base FORWARD (i=1)
        #print(seq[1:])
        #print(re.findall(r"ATG.*TAA", seq))
        orf1.append(re.findall(r"ATG.*TAA", seq[1:])) #add it to a list containing this kind of ORF
        orf_nums += 1 #add up total number of ORFs in file (all sequences)
        orf1_seq_index.append(i) #which sequence ORF1 was found in
        #To do two things at once, find the ORFs in the sequence that was randomly picked and add it to it's corresponding ORF list
        if i == random_seq: #meaning, that's the sequence that was picked at random
            orf1_random.append(re.findall(r"ATG.*TAA", seq[1:])) #add that sequence to a list containing ORF1 for that randomly picked sequence
        #print(seq[i[)
    if re.findall(r"ATG.*TAG", seq[1:]): #second possible ORF
        #print("ORF2 found ")
        orf2.append(re.findall(r"ATG.*TAG", seq[1:]))
        orf_nums += 1
        orf2_seq_index.append(i)
        if i == random_seq:
            orf2_random.append(re.findall(r"ATG.*TAG", seq[1:]))
        #print(seq[i])
    if re.findall(r"ATG.*TGA", seq[1:]): #third possible ORF
        #print("ORF3 found ")
        orf3.append(re.findall(r"ATG.*TGA", seq[1:]))
        orf_nums += 1
        orf3_seq_index.append(i)
        if i == random_seq:
            orf3_random.append(re.findall(r"ATG.*TGA", seq[1:]))
        #print(seq[i])
    else:
        print("No ORFs found ")

#Store max length values of the ORFs in the sequence that was randomly chosen
max1 = []
max2 = []
max3 = []

#print(random_seq)
#This index i will NOT change with the reading frame
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
        print("Sequence %d does not have any ORF3" % random_seq)
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

#Now find the start and end positions of the ORFs
#Change seq range to seq[1:] so the reading frame advances 1
for i in range(len(seqs_strings)): #goes through all 25 sequences
    if i == random_seq:
        seq = seqs_strings[i]
        for m in re.finditer(r"ATG.*TAA", seq[1:]):  # ORF1
            rand_start1 = m.start() #Get start position
            rand_end1 = m.end() #Get end position
            print("ORF1 sequence %d:" % random_seq)
            print("Start:", rand_start1)
            print("End:", rand_end1)
            print(rand_end1 - rand_start1, "\n")

        for m in re.finditer(r"ATG.*TAG", seq[1:]): #ORF2
            rand_start2 = m.start()
            rand_end2 = m.end()
            print("ORF2 sequence %d:" % random_seq)
            print("Start:", rand_start2)
            print("End:", rand_end2)
            print(rand_end2-rand_start2, "\n")

        for m in re.finditer("ATG.*TGA", seq[1:]): #ORF 3
            rand_start3 = m.start()
            rand_end3 = m.end()
            print("ORF3 sequence %d:" % random_seq)
            print("Start:", rand_start3)
            print("End:", rand_end3)
            print(rand_end3-rand_start3, "\n")
    else:
        continue


print("\nFor reading frame 2, with sequences starting at index 1:")
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
print("\nLongest ORF1 in reading frame 2 is %d bases long in sequence %d" % (max_orf1, max_ORF_index1))

#ORF 2
for i in range(len(orf2)):
    for j in range(len(orf2[i])):
        orf2_lengths.append(len(orf2[i][j]))
    max_orf2 = max(orf2_lengths)
    if i == orf2_seq_index[i]:
        max_ORF_index2 = orf2_lengths.index(max_orf2)
    #print(max_ORF_index2)
print("Longest ORF2 in reading frame 2 is %d bases long in sequence %d" % (max_orf2, max_ORF_index2))

#ORF 3
for i in range(len(orf3)):
    for j in range(len(orf3[i])):
        orf3_lengths.append(len(orf3[i][j]))
    max_orf3 = max(orf3_lengths)
    if i == orf3_seq_index[i]:
        max_ORF_index3 = orf3_lengths.index(max_orf3)
    #print(max_ORF_index3)
print("Longest ORF3 in reading frame 2 is %d bases long in sequence %d" % (max_orf3, max_ORF_index3))


#Find longest ORF of all of them:
longest_orf = max(max_orf1, max_orf2, max_orf3)
#Check to see which ORF group has the longest overall ORF and assign correct index
if longest_orf == max_orf1:
    longest_orf_index = max_ORF_index1
if longest_orf == max_orf2:
    longest_orf_index = max_ORF_index2
if longest_orf == max_orf3:
    longest_orf_index = max_ORF_index3

print("\nLongest ORF in file with reading frame 2 is %d base pairs long in sequence %d" % (longest_orf, longest_orf_index))
print("Sequence identifier of sequence %d:" % (longest_orf_index))
print(identifiers[longest_orf_index])





print("\n\n...............................................................................................................")
print("\nReading frame 3 (+2)")




#4.) Find all repeats of length n
#Use regex (most likely)
#SO any sequence of length n, NOT a SPECIFIC repeat
#I know these FASTA files are correct (i.e. no other bases besides the usual 4)
#So it would be something like "/w{n}" for all alphanumeric values n times
#So like, if you want to find repeats where n = 3
#Say you have the sequence: ATTTAAGAGAGAGTTTTATTTATTTATTTTG
#TTT is one, GAG is another, ATT is another, etc.
#HAVE to make sure it allows for OVERLAPPING repeats!!!!
#"count_overlap" function in BioSeq might work.......????

