# app/sequence_utils.py

from Bio.Seq import Seq
from itertools import product
from typing import List, Tuple
import re
import random

def generate_iupac_permutations(iupac: str, iupac_to_dna: dict) -> List[str]:
    # For each IUPAC base, generate the list of possible nucleotides
    bases = [iupac_to_dna[base] for base in iupac]
    # Generate all possible combinations of nucleotides
    permutations = [''.join(p) for p in product(*bases)]
    return permutations

def is_near_match(sequence: str, motif: str, iupac_to_dna: dict, max_diff: int) -> bool:
    """Check if a sequence is within `max_diff` mutations from the motif."""
    diff = 0
    for i, base in enumerate(motif):
        if sequence[i] not in iupac_to_dna[base]:
            diff += 1
        if diff > max_diff:
            return False
    return True



"""
Copyright (c) 2024 VECTOR FUTURES LTD
All rights reserved.
This file is part of the Vector Futures Codon Optimization App and may not be copied, distributed, or modified without express written permission.
"""