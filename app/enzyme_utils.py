# app/enzyme_utils.py

import operator
import argparse
from io import StringIO, BytesIO
from socket import inet_aton
from Bio import SeqIO
import operator
import pandas as pd
import textwrap
import numpy as np
import re
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from collections import defaultdict
from tqdm import tqdm

def GC(dna_sequence):
    # GC content calculation
    dna_sequence = dna_sequence.upper()
    gc_count = dna_sequence.count('G') + dna_sequence.count('C')
    gc_content = (gc_count / len(dna_sequence)) * 100
    return gc_content

def not_in(motifs, sequence):
    for i in motifs:
        if i in Seq(sequence):
            return(False)
            break
        else:
            continue
    return(True)

def motif_remove(site_list, sequence):
    x = MutableSeq(sequence)
    for site in site_list:
        while site in x:
            start = x.find(site)
            locus = start-start%3
            code = str(x)[locus:locus+6]
            codons = textwrap.wrap(code,3)
            mutable_codons = [MutableSeq(items) for items in codons]
            dntps = ['G','C','A','T']
            for nt in dntps:
                mutable_codons[0][2] = nt
                new_code = "".join([str(codon) for codon in mutable_codons])
                if Seq(new_code).translate() == Seq(x[locus:locus+6]).translate() and nt != codons[0][2] and not_in(site_list, new_code):
                    break
            else:
                for nt in dntps:
                    mutable_codons[1][2] = nt
                    new_code = "".join([str(codon) for codon in mutable_codons])
                    if Seq(new_code).translate() == Seq(x[locus:locus+6]).translate() and nt != codons[1][2] and not_in(site_list, new_code):
                        break
            x[locus:locus+6] = new_code
    return(str(x))

def is_CpG_island(DNA):
    gc = GC(DNA)
    if (DNA.count('C') * DNA.count('G'))/ len(DNA) == 0:
        return False
    else:
        obs = DNA.count('CG')
        exp = (DNA.count('C') * DNA.count('G'))/ len(DNA) # Gardiner-Garden et al. 1987 PMID 3656447
        #exp = (DNA.count('C') + DNA.count('G')/2)**2 / len(DNA) # Saxonov et al. 2006 PMID 16432200
        return(bool(gc > 50 and obs/exp > 0.6))

# Below functions act to minimise prevalence of CpG islands
def CpG_check(sequence):
    for index, base in enumerate(sequence):
        g = sequence[index:index+210]
        if is_CpG_island(g):
            count += 1
    return count/(len(sequence)/210)

def CpG_remove(sequence,codon_dict, cpg_depletion_level):
    if cpg_depletion_level == 'Stringent':
        factor = 0.2
    elif cpg_depletion_level == 'Moderate':
        factor = 0.07
    x = MutableSeq(sequence)
    CpG_motif = re.compile("CG")  # Define the motif as a regular expression
    # Assuming x is a MutableSeq object
    matches = CpG_motif.findall(str(x).upper())
    # Check if there are any matches
    if matches:
        # Determine the number of matches to modify
        num_matches_to_modify = int(len(matches) * factor)
        matches_to_modify = np.random.choice(matches, num_matches_to_modify, replace=False)
        # Modify the DNA sequence for selected matches
        for match in matches_to_modify:
            start = str(x).upper().find(match)
            locus = start-start%3
            code = str(x)[locus:locus+6]
            aas = Seq(code).translate()
            for i in range(24):
                t = [np.random.choice(codon_dict[i][0]) for i in aas]
                new_code = ''.join(t)
                if 'CG' not in new_code and new_code == code:
                    break
            x[locus:locus+6] = new_code
    return(str(x))

def CpG_island_remove(seq, codon_dict, cpg_depletion_level):
    sequence = MutableSeq(seq)
    peptide = Seq(seq).translate()
    if len(seq) > 210:
        for index, residue in enumerate(peptide):
            g = sequence[index*3:(index+70)*3]
            if is_CpG_island(g):
                sequence[index*3:(index+70)*3] = CpG_remove(g, codon_dict, cpg_depletion_level)
    out_seq = str(sequence)
    return out_seq


"""
Copyright (c) 2024 VECTOR FUTURES LTD
All rights reserved.
This file is part of the Vector Futures Codon Optimization App and may not be copied, distributed, or modified without express written permission.
"""