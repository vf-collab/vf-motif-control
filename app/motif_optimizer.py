# app/motif_optimizer.py

from Bio.Seq import Seq
from itertools import product
from typing import List, Tuple
import re
import random
from .sequence_utils import *
from .progress_utils import update_progress_bar


def insert_motif(dna_sequence: str, iupac: str, max_diff: int = 2) -> Tuple[str, int, int]:
    dntps = ['G', 'C', 'A', 'T']
    
    # Create an IUPAC dictionary for generating permutations
    iupac_to_dna = {
        "A": ["A"],  # Adenine
        "C": ["C"],  # Cytosine
        "G": ["G"],  # Guanine
        "T": ["T"],  # Thymine
        "R": ["A", "G"],  # Purine (A or G)
        "Y": ["C", "T"],  # Pyrimidine (C or T)
        "S": ["G", "C"],  # Strong interaction (G or C)
        "W": ["A", "T"],  # Weak interaction (A or T)
        "K": ["G", "T"],  # Keto (G or T)
        "M": ["A", "C"],  # Amino (A or C)
        "B": ["C", "G", "T"],  # Not A (C or G or T)
        "D": ["A", "G", "T"],  # Not C (A or G or T)
        "H": ["A", "C", "T"],  # Not G (A or C or T)
        "V": ["A", "C", "G"],  # Not T (A or C or G)
        "N": ["A", "C", "G", "T"]  # Any Nucleotide (A or C or G or T)
    }

    # Generate all possible DNA permutations that fit the IUPAC motif
    motif_permutations = generate_iupac_permutations(iupac, iupac_to_dna)
    
    # Build regex pattern for the IUPAC motif
    motif_regex = ''.join([f"[{''.join(iupac_to_dna[base])}]" if len(iupac_to_dna[base]) > 1 else iupac_to_dna[base][0] for base in iupac])
    
    # Find exact matches for the motif (before any mutations)
    exact_matches = len(re.findall(motif_regex, dna_sequence))

    # Find near matches (that are up to `max_diff` mutations away)
    potential_matches = [match.start() for match in re.finditer('.{' + str(len(iupac)) + '}', dna_sequence)
                         if is_near_match(dna_sequence[match.start():match.start() + len(iupac)], iupac, iupac_to_dna, max_diff)]
    
    print(f'Found {len(potential_matches)} possible near-match sites to mutate')

    count = 0
    mutated_sites = []
    modified_sequence = dna_sequence

    # Iterate through potential matches
    for site in potential_matches:
        # Shuffle permutations each time to ensure no bias towards the first permutations in the list
        random.shuffle(motif_permutations)
        # Attempt to apply each permutation to the site
        for permutation in motif_permutations:
            new_seq = modified_sequence[:site] + permutation + modified_sequence[site+len(permutation):]
            if Seq(new_seq).translate() == Seq(dna_sequence).translate():
                modified_sequence = new_seq
                mutated_sites.append(f"Site {site} mutated to {permutation}")
                count += 1
                break

    # Count new exact matches after mutation (in the modified sequence)
    new_exact_matches = len(re.findall(motif_regex, modified_sequence))
    
    print(f"Mutated {count} sites")
    print(f"Mutated sites include: {', '.join(mutated_sites)}")
    print(f"Initial motifs: {exact_matches}, New motifs: {new_exact_matches}")
    
    return modified_sequence, exact_matches, new_exact_matches


def destroy_motif(dna_sequence, iupac):
    dna_sequence = dna_sequence.upper()
    dntps = ['G', 'C', 'A', 'T']
    # Create an iupac dictionary
    iupac_to_regex = {
        "A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]",
        "S": "[GC]", "W": "[AT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]",
        "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"
    }
    motif_regex = "".join([iupac_to_regex.get(base, "") for base in iupac])
    # Find all exact matches for the motif
    exact_matches = len(re.findall(motif_regex, dna_sequence))
    print('Found ' + str(exact_matches) + ' motif sites')
    # Find all potential matches
    potential_matches = [match.start() for match in re.finditer(motif_regex, dna_sequence)]
    if len(potential_matches) == 0:
        print("Motif not found")
        return dna_sequence, None, None
    count = 0
    mutated_sites = []
    # Mutate matches in reverse to avoid shifting issues
    for site in reversed(potential_matches):
        match = dna_sequence[site:site + len(iupac)]
        for i, base in enumerate(match):
            A_site = site + i
            base_regex = iupac_to_regex.get(base, "")
            alt_dntps = [i for i in dntps if i not in base_regex]
            if len(alt_dntps) == 0:
                continue
            else:
                for dntp in alt_dntps:
                    modified_sequence = dna_sequence[:A_site] + dntp + dna_sequence[A_site + 1:]   
                    # Return modified sequence if protein sequence matches original
                    if Seq(modified_sequence).translate() == Seq(dna_sequence).translate():
                        dna_sequence = modified_sequence
                        count += 1
                        mutated_sites.append(str(A_site))
                        break
    # Report on outcomes
    print("Mutated " + str(count) + " bases")
    print("Mutated sites include " + ", ".join(mutated_sites))
    new_exact_matches = len(re.findall(motif_regex, dna_sequence))
    print("Initial motifs: ", exact_matches, "Remaining motifs: ", new_exact_matches)
    return dna_sequence, exact_matches, new_exact_matches


def count_motif(dna_sequence, iupac):
    dna_sequence = dna_sequence.upper()
    dntps = ['G', 'C', 'A', 'T']
    # Create an iupac dictionary
    iupac_to_regex = {
        "A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]",
        "S": "[GC]", "W": "[AT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]",
        "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"
    }
    motif_regex = "".join([iupac_to_regex.get(base, "") for base in iupac])
    # Find all exact matches for the motif
    return len(re.findall(motif_regex, dna_sequence))


"""
Copyright (c) 2024 VECTOR FUTURES LTD
All rights reserved.
This file is part of the Vector Futures Codon Optimization App and may not be copied, distributed, or modified without express written permission.
"""