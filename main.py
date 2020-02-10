# Bioinformatics:: Motif Finder
# A transcription factor regulates a gene by binding to a specific short DNA interval called a regulatory motif, or transcription factor binding site, in the gene’s upstream region, a 600-1000 nucleotide-long region preceding the start of the gene.

import random

# Creates a matrix counting the number of occurrences of each nucleotide in each column of the motif matrix.
def Count(Motifs):
    count = {} 
    col = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(col):
             count[symbol].append(0)
    row = len(Motifs)
    for i in range(row):
        for j in range(col):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

# Takes a list of strings Motifs as input and returns the count matrix of Motifs with pseudocounts as a dictionary of lists.
def CountWithPseudocounts(Motifs):
    count = {} 
    col = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(col):
             count[symbol].append(1)
    row = len(Motifs)
    for i in range(row):
        for j in range(col):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

# Profile(Motifs) that takes Motifs as input and returns their profile matrix as a dictionary of lists.
def Profile(Motifs):
    profile = Count(Motifs)

    for key, entry in profile.items():
        for i in range(len(entry)):
            entry[i] = entry[i]/len(Motifs)

    return profile

# Takes a list of strings Motifs as input and returns the profile matrix of Motifs with pseudocounts as a dictionary of lists.
def ProfileWithPseudocounts(Motifs):
    profile = CountWithPseudocounts(Motifs)
    for key, entry in profile.items():
        for i in range(len(entry)):
            entry[i] = entry[i]/(len(Motifs)+4)
    return profile

# Input:  A set of kmers Motifs.
# Output: A consensus string of Motifs.
def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol

    return consensus

def HammingDistance(p, q):
    count = 0
    for i in range(0, len(p)):
        if p[i] != q[i]:
            count += 1
    return count

# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.  Smaller score indicate that the sequence is the more accurate motif.
def Score(Motifs):
    # Insert code here
    score = 0
    for i in range(len(Motifs)):
        score += HammingDistance(Motifs[i], Consensus(Motifs))
    return score

def compute(profile):
  v = 1
  for k,e in profile.items():
    print(sum(e))
    v = sum(e) * v
  print(v)

# Determines the overall probability of a sequence for a motif based on individual base probability profiles
def Pr(Text, Profile):
    prob = 1
    for i in range(len(Text)):
        prob *= Profile[Text[i]][i]
    return prob

# Given profile matrix, this code computes the probability of every k-mer in a string Text and find a Profile-most probable k-mer in Text.
def ProfileMostProbableKmer(text, k, profile):
        p = {}
        for i in range(0,len(text)-k):
            pattern = text[i:i+k]
            prob = Pr(pattern, profile)
            p[pattern] = prob
        return max(p,key=p.get)

# Dind the set of motifs across a number of DNA sequences that match each other  most closely.
# k = k-mer length; t = number of strings in Dna.
def GreedyMotifSearch(Dna,k,t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for m in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][m:m+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

# Note to change consensus function with pseudocount call.
def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = [] # output variable
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for m in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][m:m+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

# Takes a profile matrix Profile corresponding to a list of strings Dna as input and returns a list of the Profile-most probable k-mers in each string from Dna.
def Motifs(Profile, Dna, k):
  result = []
  for i in range(len(Dna)):
        result.append(ProfileMostProbableKmer(Dna[i], k, Profile))
  return result

# Random method of motif fnding, uses random.randint to choose a random k-mer from each of t different strings Dna, and returns a list of t strings.
def RandomMotifs(Dna, k, t):
    result = []
    for i in range(t):
        idx = random.randint(0,len(Dna[0])-k)
        result.append(Dna[i][idx:idx+k])
    return result

# Conducts a randomized motif search, where Dna is a series of genomic sequence.
def RandomizedMotifSearch(Dna, k, t):
    M = RandomMotifs(Dna, k, t) 
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna, k)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs 

# Input: A dictionary Probabilities, where keys are k-mers and values are the probabilities of these k-mers (which do not necessarily sum up to 1) Output: A normalized dictionary where the probability of each k-mer was divided by the sum of all k-mers' probabilities
def Normalize(Probabilities):
    total = sum(Probabilities.values())
    for key, entry in Probabilities.items():
        Probabilities[key] = entry / total
    return Probabilities

# Input:  A dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these kmers
# Output: A randomly chosen k-mer with respect to the values in Probabilities
def WeightedDie(Probabilities):
    n = random.uniform(0, 1)
    for p in Probabilities:
        n -= Probabilities[p]
        if n <= 0:
            return p

# Takes a string Text, a profile matrix profile , and an integer k as input.  Returns a randomly generated k-mer from Text whose probabilities are generated from profile.
def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {} 
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)

# Another method of motif finding using random sampling to create the probabilities profile
def GibbsSampler(Dna, k, t, N):
    Motifs = RandomMotifs(Dna, k, t)
    BestMotifs = Motifs
    for i in range(1,N):
      num = random.randint(0,t-1)
      m = []
      for j in Motifs:
            if j != Motifs[num]:
                m.append(j)
      profile = ProfileWithPseudocounts(m)
      Motifs[num] = ProfileGeneratedString(Dna[num], profile, k)
      if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

