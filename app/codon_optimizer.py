# app/codon_optimizer.py

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
from .sequence_utils import not_in, tuple_type, validate_input, GC, check_correct, find_expansion, opt_expansion, repeat_checker, find_hairpin_repeat, prime, back_translate, repeat_remove, hairpin_remove, screen_sequence
from .enzyme_utils import motif_remove, is_CpG_island, CpG_check, CpG_remove, CpG_island_remove
from .progress_utils import update_progress_bar

def codon_optimization(uploaded_seq_file, cpg_depletion_level, codon_bias_or_GC, pattern, pattern2, codon_df, selected_GC_opt, expedite=True):
    # Create codon_dict from codon_df within the function
    codon_dict = {}
    for index, row in codon_df.iterrows():
        aa = row['Single Letter']  # Access the amino acid letter
        codon = row['Codon']       # Access the codon
        ratio = row['Ratio']       # Access the ratio
        # If the amino acid is not yet in the dictionary, add it
        if aa not in codon_dict:
            codon_dict[aa] = ([], [])
        
        # Append the codon and its ratio to the corresponding lists
        codon_dict[aa][0].append(codon)
        codon_dict[aa][1].append(ratio)
    
    # Validate input.
    input_peptide = validate_input(uploaded_seq_file)

    # Optimise how many sequences are processed
    if expedite == True:
        if len(input_peptide) < 300:
            iteration_number = 24
        elif 300 < len(input_peptide) < 900:
            iteration_number = 12
        elif 900 < len(input_peptide) < 1500:
            iteration_number = 6
        else:
            iteration_number = 3
    else:
        iteration_number = 24
    
    if cpg_depletion_level != 'None':
        iteration_number = iteration_number//3
    
    total_steps = 4 + iteration_number-1

    print("Codon optimization initializing")
    if len(input_peptide) > 1200:
        print("Large sequences may take several minutes depending on complexity")
    
    result = None
    co_list = []  # Codon-optimized sequences
    scores = np.array([])  # To hold scores for each sequence
    size = 0  # Initialize size before the loop

    with tqdm(total=100, desc="Codon Optimization", unit="%", bar_format="{l_bar}{bar} [progress: {percentage:.1f}%]") as pbar:
        update_progress_bar(pbar, total_steps)  # Step 1
        while size < iteration_number:
            if len(input_peptide) < 300:
                x = back_translate(input_peptide, codon_dict, pattern, pattern2, selected_GC_opt)
                repeated_motifs = repeat_checker(x)
                modified_sequence = repeat_remove(repeated_motifs, x, codon_dict, pattern, pattern2)
                modified_sequence = hairpin_remove(modified_sequence, codon_dict, pattern, pattern2)
                if modified_sequence and cpg_depletion_level != "None":
                    modified_sequence = CpG_island_remove(modified_sequence, codon_dict, cpg_depletion_level)
                if Seq(modified_sequence).translate() != input_peptide:
                    modified_sequence = check_correct(input_peptide, Seq(str(modified_sequence)).translate(), modified_sequence, codon_dict, pattern, pattern2)
                if not_in(pattern, modified_sequence) == False:
                    modified_sequence = motif_remove(pattern, modified_sequence)
                boo, score = screen_sequence(modified_sequence)
                if boo == True:
                    co_list.append(modified_sequence)
                    scores = np.append(scores, score)
                    size+=1
                    update_progress_bar(pbar, total_steps)  # Step 2
            else:
                chunk_size = 300-(300%3)
                input_peptide_split = textwrap.wrap(str(input_peptide), chunk_size)
                input_peptide_unsplit = []
                for split in input_peptide_split:
                    x = back_translate(split, codon_dict, pattern, pattern2, selected_GC_opt)
                    repeated_motifs = repeat_checker(x)
                    modified_sequence = repeat_remove(repeated_motifs, x, codon_dict, pattern, pattern2)
                    modified_sequence = hairpin_remove(modified_sequence, codon_dict, pattern, pattern2)
                    if modified_sequence and cpg_depletion_level != "None":
                        modified_sequence = CpG_island_remove(modified_sequence, codon_dict, cpg_depletion_level)
                    if Seq(modified_sequence).translate() != split:
                        modified_sequence = check_correct(split, Seq(str(modified_sequence)).translate(), modified_sequence, codon_dict, pattern, pattern2)
                    if not_in(pattern, modified_sequence) == False:
                        modified_sequence = motif_remove(pattern, modified_sequence)
                    input_peptide_unsplit.append(modified_sequence)
                produce = ''.join(input_peptide_unsplit)
                boo, score = screen_sequence(produce, stringent=False, quiet=True)
                if boo == True:
                    co_list.append(produce)
                    scores = np.append(scores, score)
                    size+=1
                    update_progress_bar(pbar, total_steps)  # Step 2
                    

        # next populate a database with results, sort by scores and select highest value
        d = {'Sequence':co_list, 'Score':scores}
        df = pd.DataFrame(d, columns=['Sequence','Score'])
        df2 = df[df['Score'] < np.inf]
        df_sorted = df2.sort_values(by=['Score'], ascending=True)
        df_sorted = df_sorted.reset_index(drop=True)
        print(df_sorted)
        update_progress_bar(pbar, total_steps)  # Step 3

        result1 = ''
        for i in range(len(df_sorted)):
            test = df_sorted['Sequence'][i]
            if Seq(test).translate() == input_peptide and not_in(pattern,str(test)):
                result1 = test
                break
            elif Seq(test).translate() == input_peptide and not_in(pattern2,str(test)):
                result1 = test
                break
            
        if len(result1) > 0:
            update_progress_bar(pbar, total_steps)  # Step 4
            return result1
        else:
            update_progress_bar(pbar, total_steps)  # Step 4
            return None

"""
Copyright (c) 2024 VECTOR FUTURES LTD
All rights reserved.
This file is part of the Vector Futures Codon Optimization App and may not be copied, distributed, or modified without express written permission.
"""