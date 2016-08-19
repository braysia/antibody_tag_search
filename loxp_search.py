from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import re
# coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", generic_dna)
# coding_dna.translate()
from itertools import product
from operator import itemgetter
import numpy as np
import matplotlib.pyplot as plt

left = ['TACCGTTCGTATA', 'AATTATTCGTATA']
spacer = ['GCATACAT', 'GTATACAT']
right = ['TATACGAACGGTA', 'TATACGAATACCT']

products = product(left, spacer, right)
possible_seq = [''.join(i) for i in products]

REF_FP = 'MVSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTFGYGVACFSRYPDHMKQHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSHQSALSKDPNEKRDHMVLLEFVTAA'


def translate_seq(seq):
    sequences = Seq(seq, generic_dna)
    amino_seq = str(sequences.translate())
    # amino_seq = re.search('(?P<mag>.*)\*', amino_seq).group('mag') # remove the stop codon
    return amino_seq

def score_seq(ref_seq, num, cand_seq):
    scores = []
    acids = []
    for i in range(3):
        lox_acids = translate_seq(cand_seq[0:])
        ref = ref_seq[num:num+len(lox_acids)]
        scores.append(len([i for i, ii in zip(ref, lox_acids) if i==ii]))
        acids.append(lox_acids)
    max_idx = scores.index(max(scores))
    return scores[max_idx], acids[max_idx]


def calculate_score(ref_seq, num):
    scores = []
    acid_store = []
    for cand_seq in possible_seq:
        each_score, each_acid = score_seq(ref_seq, num, cand_seq)
        scores.append(each_score)
        acid_store.append(each_acid)
    max_idx = scores.index(max(scores))
    return scores[max_idx], acid_store[max_idx], possible_seq[max_idx]

# index, element = max(enumerate(items), key=itemgetter(1))

def search_fp(ref_fp):
    scores = []
    acid_store = []
    seq_store = []
    for num in range(len(ref_fp)-11):
        each_score, each_acid, each_seq = calculate_score(ref_fp, num)
        scores.append(each_score)
        acid_store.append(each_acid)
        seq_store.append(each_seq)
    arr_scores = np.array(scores)
    idx = np.where(arr_scores == arr_scores.max())[0]
    for i in idx:
        print '{0} matches, {1}, {2} at {3}'.format(scores[i], acid_store[i], seq_store[i], i)

search_fp(REF_FP)


# scores = []
# acid_store = []
# seq_store = []
# for num in range(len(REF_FP)-11):
#     each_score, each_acid, each_seq = calculate_score(REF_FP, num)
#     scores.append(each_score)
#     acid_store.append(each_acid)
#     seq_store.append(each_seq)
# arr_scores = np.array(scores)
# idx = np.where(arr_scores == arr_scores.max())[0]
# for i in idx:
#     print '{0} matches, {1}, {2} at {3}'.format(scores[i], acid_store[i], seq_store[i], i)
#
# plt.plot(scores)
# plt.show()
