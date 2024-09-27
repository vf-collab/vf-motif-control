# app/sequence_utils.py

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
from Bio.Seq import Seq
from .enzyme_utils import motif_remove, not_in


# create a set of functions for characterising sequences
def tuple_type(s):
    try:
        x, y = map(int, s.split(','))
        return (x, y)
    except:
        raise argparse.ArgumentTypeError("Tuples must be x,y")

def validate_input(input):
    na = set('ATGC')
    seq_object = Seq(input)
    if not set(str(seq_object)) - na:
        return(str(seq_object.translate()).upper())
    else:
        return(str(seq_object).upper())

def GC(dna_sequence):
    # GC content calculation
    dna_sequence = dna_sequence.upper()
    gc_count = dna_sequence.count('G') + dna_sequence.count('C')
    gc_content = (gc_count / len(dna_sequence)) * 100
    return gc_content

def check_correct(inpep, outpep, outDNA, codon_dict, pattern, pattern2):
    x = str(outDNA)
    if len(inpep) != len(outpep):
        n = len(inpep) - len(outpep)
        outpep = str(outpep) + '_'*n
    for i, item in enumerate(inpep):
        if inpep[i] != outpep[i]:
            aa = dict(zip(codon_dict[inpep[i]][0],codon_dict[inpep[i]][1]))
            s_aa = sorted(aa.items(), key=operator.itemgetter(1), reverse=True)
            for f in range(6):
                am = s_aa[f%len(s_aa)][0]
                if not_in(pattern, x[(i-1)*3:(i*3)] + am + x[(i+1)*3:(i+2)*3]):
                    x = x[:(i*3)] + am + x[(i+1)*3:]
                    break
            else:
                for t in range(6):
                    ap = np.random.choice(codon_dict[item][0])
                    if not_in(pattern, x[(i-1)*3:(i*3)] + ap + x[(i+1)*3:(i+2)*3]):
                        x = x[:(i*3)] + ap + x[(i+1)*3:]
                        break
    return(x)

    # Write a function to find any polymer expansions and reoptimize sequence to minimise repetition
def find_expansion(seq):
    pattn = re.compile(r'(\w)\1{4,}')  # matches a character that repeats at least 5 times
    for match in re.finditer(pattn, seq):
        start = match.start()
        end = match.end()
        yield (start, end, match.group())

def repeat_checker(sequence):
    min_repeat_length = 8
    repeats = rf"(.{{{min_repeat_length},}})(?=(.*\1))"
    repeated_sequences = []
    current_position = 0
    while True:
        match = re.search(repeats, sequence[current_position:])
        if not match:
            break
        match_start = current_position + match.start()
        repeated_sequences.append(match_start)
        current_position = match_start + 1
    return(repeated_sequences)

def find_hairpin_repeat(seq):
    for i, base in enumerate(seq):
        sh = seq[i:i+8]
        l = seq[i+9:]
        if str(sh) in str(l) or str(Seq(sh).reverse_complement()) in str(l):
            start = seq.find(sh)
            end = start+len(seq)
            match = sh
            yield (start, end, match)
    
def prime(x):
    if x < 2 or x%1 != 0:
        return False
    else:
        bool = [x % i == 0 for i in range(2, x)]
        return sum(bool) == 0

def opt_expansion(inpep, outDNA, codon_dict, pattern, pattern2):
    x = MutableSeq(outDNA)
    if find_expansion(inpep):
        for repeat in find_expansion(inpep):
            start_DNA = repeat[0]*3
            end_DNA = repeat[1]*3
            coding_seq = str(x)[start_DNA:end_DNA]
            peptide = Seq(coding_seq).translate()
            streeng = ''
            for i, item in enumerate(peptide):
                if i == 0:
                    t = 0
                elif i%4 == 0:
                    t = 3
                elif i%3 == 0:
                    t = 2
                elif prime(i):
                    t = 0
                else:
                    t = 1
                a = dict(zip(codon_dict[item][0],codon_dict[item][1]))
                sorted_a = sorted(a.items(), key=operator.itemgetter(1), reverse=True)
                if not_in(pattern, streeng + sorted_a[t%len(sorted_a)][0]):
                    streeng = streeng + sorted_a[t%len(sorted_a)][0]
                else:      
                    for list_item in sorted_a:
                        if not_in(pattern, streeng + list_item[0]):
                            streeng = streeng + list_item[0]
            x[start_DNA:end_DNA] = streeng
    return(str(x))    

# for back translation
def back_translate(peptide_seq, codon_dict, pattern, pattern2, selected_GC_opt):
    chunks = textwrap.wrap(peptide_seq,16)
    new = ''
    for chunk in chunks:
        t = [np.random.choice(codon_dict[i][0],1,p=codon_dict[i][1])[0] for i in chunk]
        t = ''.join(t)
        if not_in(pattern, str(t)) and Seq(t).translate() == chunk and (selected_GC_opt[0]-4) <= GC(str(t)) <= (selected_GC_opt[1]+6):
            new = [new, t]
            new = ''.join(new)
        else:
            for q in range(100):
                t = [np.random.choice(codon_dict[i][0],1,p=codon_dict[i][1])[0] for i in chunk]
                t = ''.join(t)
                if not_in(pattern, str(t)) and Seq(t).translate() == chunk and (selected_GC_opt[0]-12) <= GC(str(t)) <= (selected_GC_opt[1]+18):
                    new = [new, t]
                    new = ''.join(new)
                    break
    new = ''.join([str(item) for item in new])
    new = opt_expansion(peptide_seq, new, codon_dict, pattern, pattern2)
    if peptide_seq != Seq(new).translate():
        new = check_correct(peptide_seq, Seq(new).translate(), new, codon_dict, pattern, pattern2)
        new = new[:len(peptide_seq)*3]
    new = motif_remove(pattern, new)
    return(new)

def repeat_remove(repeats, sequence, codon_dict, pattern, pattern2):
    x = MutableSeq(sequence)
    for loc in repeats:
        locus = loc - loc % 3
        code = x[locus:locus+6]
        codons = textwrap.wrap(str(code), 3)
        dntps = ['G', 'C', 'A', 'T']
        mutable_codons = [MutableSeq(items) for items in codons]
        for nt in dntps:
            mutable_codons[0][2] = nt
            if Seq(str(mutable_codons[0])).translate() == Seq(codons[0]).translate() and nt !=codons[0][2]:
                x[locus:locus+6] = "".join(str(codon) for codon in mutable_codons)
                break
        else:
            for nt2 in dntps:
                mutable_codons[1][2] = nt2
                if Seq(str(mutable_codons[1])).translate() == Seq(codons[1]).translate() and nt2 != codons[1][2]:
                        x[locus:locus+6] = "".join(str(codon) for codon in mutable_codons)
                        break
    if Seq(sequence).translate() != Seq(str(x)).translate():
        x = check_correct(Seq(sequence).translate(), Seq(str(x)).translate(), x, codon_dict, pattern, pattern2)
        x = x[:len(Seq(sequence).translate())*3]
    x = motif_remove(pattern, x)
    return(str(x))

def hairpin_remove(sequence, codon_dict, pattern, pattern2):
    x = MutableSeq(sequence)
    for i, base in enumerate(x):
        sh = str(x[i:i+8])
        l = str(x[i+9:])
        if str(sh) in str(Seq(l).reverse_complement()):
            start = str(x).find(sh)
            locus = start - start % 3
            code = x[locus:locus+6]
            codons = textwrap.wrap(str(code), 3)
            dntps = ['G', 'C', 'A', 'T']
            mutable_codons = [MutableSeq(items) for items in codons]
            for nt in dntps:
                mutable_codons[0][2] = nt
                if Seq(str(mutable_codons[0])).translate() == Seq(codons[0]).translate() and nt !=codons[0][2]:
                    x[locus:locus+6] = "".join(str(codon) for codon in mutable_codons)
                    break
            else:
                for nt2 in dntps:
                    mutable_codons[1][2] = nt2
                    if Seq(str(mutable_codons[1])).translate() == Seq(codons[1]).translate() and nt2 != codons[1][2]:
                            x[locus:locus+6] = "".join(str(codon) for codon in mutable_codons)
                            break
    if Seq(sequence).translate() != Seq(str(x)).translate():
        x = check_correct(Seq(sequence).translate(), Seq(str(x)).translate(), x, codon_dict, pattern, pattern2)
        x = x[:len(Seq(sequence).translate())*3]
    x = motif_remove(pattern, x)
    return(str(x))

# below function scores each sequence
def screen_sequence(seq, window_size=100, stringent=True, quiet=True):

    if stringent:
        score_threshold = 10
    else:
        score_threshold = 40


    def check_repeat_coverage(sequence, min_repeat_length=8, threshold_percentage=40):
        repeat_score = 0
        repeat_issues = []
        repeat_positions = set()  # To store the positions of bases that are part of repeats
        sequence_length = len(sequence)
        
        # Dictionary to count occurrences of subsequences
        subseq_count = defaultdict(int)
        
        # Use a sliding window to find repeated sequences
        for repeat_length in range(min_repeat_length, sequence_length // 2):
            for i in range(sequence_length - repeat_length + 1):
                subseq = sequence[i:i + repeat_length]
                subseq_count[subseq] += 1  # Increment count of subsequence
                
                # Only consider subsequences that appear more than once
                if subseq_count[subseq] > 1:
                    for j in range(i, i + repeat_length):
                        repeat_positions.add(j)

        # Calculate the percentage of the sequence that is part of repeated regions
        repeat_base_count = len(repeat_positions)
        repeat_percentage = (repeat_base_count / sequence_length) * 100
        
        if repeat_percentage > threshold_percentage:
            repeat_score += 2
            repeat_issues.append(f"Approximately {repeat_percentage:.1f}% of the sequence is formed of repeats greater than {min_repeat_length} bases.")
        
        return repeat_score, repeat_issues

    # Function to check for short repeats in sequence
    def check_short_repeats(window):
        repeat_score = 0
        repeat_issues = []
        # Example: Find simple repeats of length 5 or more
        for i in range(len(window) - 4):
            subseq = window[i:i+5]
            if window.count(subseq) > 7:  # If the subsequence appears more than 7 times
                repeat_score += 0.2  # Add to the score for each instance
                repeat_issues.append(f"Short repeat {subseq} found")# in window: {window}")
        return repeat_score, repeat_issues

    # Function to check for short repeats in sequence
    def check_long_repeats(window):
        repeat_score = 0
        repeat_issues = []
        # Example: Find long repeats of length 10 or more
        for i in range(len(window) - 9):
            subseq = window[i:i+10]
            if window.count(subseq) > 1:  # If the subsequence appears more than once
                repeat_score += 2  # Add to the score for each instance
                repeat_issues.append(f"Long repeat {subseq} found")# in window: {window}")
        return repeat_score, repeat_issues

    # Function to check for hairpins in a window
    def check_hairpins(window):
        hairpin_score = 0
        hairpin_issues = []
        # Using a sliding window to find hairpins
        for i in range(len(window) - 9):
            sh = window[i:i+10]  # 10-base stem
            l = window[i+11:]    # Remainder of the sequence
            if str(sh) in str(l) or str(Seq(sh).reverse_complement()) in str(l):
                hairpin_score += 1  # Penalize more heavily for hairpins
                hairpin_issues.append(f"Hairpin found: {sh}")# in window: {window}")
        return hairpin_score, hairpin_issues

    # Function to calculate GC richness
    def GC(dna_sequence):
        # Convert the sequence to uppercase to handle any lowercase input
        dna_sequence = dna_sequence.upper()
        # Count the number of G's and C's in the sequence
        gc_count = dna_sequence.count('G') + dna_sequence.count('C')
        # Calculate GC content percentage
        gc_content = (gc_count / len(dna_sequence)) * 100
        return gc_content

    # Function to check for GC-rich regions
    def check_gc_content(window, gc_threshold=70):
        gc_score = 0
        gc_issues = []
        gc_content = GC(window)
        if gc_content > gc_threshold:
            gc_score += 0.5  # Penalty for GC-rich regions
            gc_issues.append(f"High GC content ({gc_content}%) in window")# in window: {window}")
        return gc_score, gc_issues

    total_score = 0
    issues = []

    # Accounting for large sequences which gets computationally expensive when handling as one chunk
    if len(seq) < 1000:
        # If sequence is small, make no adjustments
        slithers = [seq]
    else:
        # For big sequences, split sequence into sections to avoid a mammoth job
        factor = int(len(seq)/1000)
        one_size = int(len(seq)/factor)
        slithers = textwrap.wrap(seq, one_size)

    for slither in slithers:
    
        # Check for short repeats
        short_repeat_score, short_repeat_issues = check_short_repeats(slither)
        total_score += short_repeat_score
        issues.extend(short_repeat_issues)

        # Check for long repeats
        long_repeat_score, long_repeat_issues = check_long_repeats(slither)
        total_score += long_repeat_score
        issues.extend(long_repeat_issues)

        # Check for repeat coverage
        repeat_coverage_score, repeat_coverage_issues = check_repeat_coverage(slither)
        total_score += repeat_coverage_score
        issues.extend(repeat_coverage_issues)

        # Check for hairpins
        hairpin_score, hairpin_issues = check_hairpins(slither)
        total_score += hairpin_score
        issues.extend(hairpin_issues)

        # Sliding window over the sequence
        for i in range(0, len(slither) - window_size + 1, window_size):
            window = slither[i:i+window_size]
            # Check for GC-rich regions
            gc_score, gc_issues = check_gc_content(window)
            total_score += gc_score
            issues.extend(gc_issues)

    # If total score exceeds threshold, reject the sequence
    if total_score > score_threshold:
        # Print results unless verbosity is quiet
        if not quiet:
            print(f"Sequence rejected. Total score: {total_score}")
            print("List of issues:")
            for issue in issues:
                print(issue)
        return False, total_score

    if not quiet:
        print(f"Sequence passed. Total score: {total_score}")
    return True, total_score

def pattern_generator(enzymes):
    Enz_dict = {'AanI': 'TTATAA', 'AarI': 'CACCTGC', 'AasI': 'GACNNNNNNGTC', 'AatII': 'GACGTC', 'Acc65I': 'GGTACC', 'AccI': 'GTMKAC', 'AciI': 'CCGC', 'AclI': 'AACGTT', 'AcuI': 'CTGAAG', 'AdeI': 'CACNNNGTG', 'AfeI': 'AGCGCT', 'AflII': 'CTTAAG', 'AflIII': 'ACRYGT', 'AgeI': 'ACCGGT', 'AhdI': 'GACNNNNNGTC', 'AjiI': 'CACGTC', 'AjuI': 'GAANNNNNNNTTGG', 'AleI': 'CACNNNNGTG', 'AlfI': 'GCANNNNNNTGC', 'AloI': 'GAACNNNNNNTCC', 'AluI': 'AGCT', 'Alw21I': 'GWGCWC', 'Alw26I': 'GTCTC', 'Alw44I': 'GTGCAC', 'AlwI': 'GGATC', 'AlwNI': 'CAGNNNCTG', 'ApaI': 'GGGCCC', 'ApaLI': 'GTGCAC', 'ApeKI': 'GCWGC', 'ApoI': 'RAATTY', 'AscI': 'GGCGCGCC', 'AseI': 'ATTAAT', 'AsiSI': 'GCGATCGC', 'AvaI': 'CYCGRG', 'AvaII': 'GGWCC', 'AvrII': 'CCTAGG', 'BaeGI': 'GKGCMC', 'BaeI': 'ACNNNNGTAYC', 'BamHI': 'GGATCC', 'BanI': 'GGYRCC', 'BanII': 'GRGCYC', 'BauI': 'CACGAG', 'BbsI': 'GAAGAC', 'BbvCI': 'CCTCAGC', 'BbvI': 'GCAGC', 'BccI': 'CCATC', 'BceAI': 'ACGGC', 'BcgI': 'CGANNNNNNTGC', 'BciVI': 'GTATCC', 'BclI': 'TGATCA', 'BcnI': 'CCSGG', 'BcoDI': 'GTCTC', 'BcuI': 'ACTAGT', 'BdaI': 'TGANNNNNNTCA', 'BfaI': 'CTAG', 'BfiI': 'ACTGGG', 'BfmI': 'CTRYAG', 'BfuAI': 'ACCTGC', 'BfuCI': 'GATC', 'BfuI': 'GTATCC', 'BglI': 'GCCNNNNNGGC', 'BglII': 'AGATCT', 'BlpI': 'GCTNAGC', 'Bme1390I': 'CCNGG', 'Bme1580I': 'GKGCMC', 'BmgBI': 'CACGTC', 'BmrI': 'ACTGGG', 'BmsI': 'GCATC', 'BmtI': 'GCTAGC', 'BoxI': 'GACNNNNGTC', 'BpiI': 'GAAGAC', 'BplI': 'GAGNNNNNCTC', 'BpmI': 'CTGGAG', 'Bpu10I': 'CCTNAGC', 'Bpu1102I': 'GCTNAGC', 'BpuEI': 'CTTGAG', 'BsaAI': 'YACGTR', 'BsaBI': 'GATNNNNATC', 'BsaHI': 'GRCGYC', 'BsaI': 'GGTCTC', 'BsaJI': 'CCNNGG', 'BsaWI': 'WCCGGW', 'BsaXI': 'ACNNNNNCTCC', 'BseDI': 'CCNNGG', 'BseGI': 'GGATG', 'BseJI': 'GATNNNNATC', 'BseLI': 'CCNNNNNNNGG', 'BseMI': 'GCAATG', 'BseMII': 'CTCAG', 'BseNI': 'ACTGG', 'BseRI': 'GAGGAG', 'BseSI': 'GKGCMC', 'BseXI': 'GCAGC', 'BseYI': 'CCCAGC', 'BsgI': 'GTGCAG', 'Bsh1236I': 'CGCG', 'Bsh1285I': 'CGRYCG', 'BshNI': 'GGYRCC', 'BshTI': 'ACCGGT', 'BsiEI': 'CGRYCG', 'BsiHKAI': 'GWGCWC', 'BsiWI': 'CGTACG', 'BslI': 'CCNNNNNNNGG', 'BsmAI': 'GTCTC', 'BsmBI': 'CGTCTC', 'BsmBI-v2': 'CGTCTC', 'BsmFI': 'GGGAC', 'BsmI': 'GAATGC', 'BsoBI': 'CYCGRG', 'Bsp119I': 'TTCGAA', 'Bsp120I': 'GGGCCC', 'Bsp1286I': 'GDGCHC', 'Bsp1407I': 'TGTACA', 'Bsp143I': 'GATC', 'Bsp68I': 'TCGCGA', 'BspCNI': 'CTCAG', 'BspDI': 'ATCGAT', 'BspEI': 'TCCGGA', 'BspHI': 'TCATGA', 'BspLI': 'GGNNCC', 'BspMI': 'ACCTGC', 'BspOI': 'GCTAGC', 'BspPI': 'GGATC', 'BspQI': 'GCTCTTC', 'BspTI': 'CTTAAG', 'BsrBI': 'CCGCTC', 'BsrDI': 'GCAATG', 'BsrFI': 'RCCGGY', 'BsrGI': 'TGTACA', 'BsrI': 'ACTGG', 'BssHII': 'GCGCGC', 'BssKI': 'CCNGG', 'BssSI-v2': 'CACGAG', 'BssSαI': 'CACGAG', 'Bst1107I': 'GTATAC', 'BstAPI': 'GCANNNNNTGC', 'BstBI': 'TTCGAA', 'BstEII': 'GGTNACC', 'BstNI': 'CCWGG', 'BstUI': 'CGCG', 'BstXI': 'CCANNNNNNTGG', 'BstYI': 'RGATCY', 'BstZ17I': 'GTATAC', 'Bsu15I': 'ATCGAT', 'Bsu36I': 'CCTNAGG', 'BsuRI': 'GGCC', 'BtgI': 'CCRYGG', 'BtgZI': 'GCGATG', 'BtsCI': 'GGATG', 'BtsI-v2': 'GCAGTG', 'BtsIMutI': 'CAGTG', 'BtsαI': 'GCAGTG', 'BveI': 'ACCTGC', 'Cac8I': 'GCNNGC', 'CaiI': 'CAGNNNCTG', 'Cfr10I': 'RCCGGY', 'Cfr13I': 'GGNCC', 'Cfr42I': 'CCGCGG', 'Cfr9I': 'CCCGGG', 'CfrI': 'YGGCCR', 'ClaI': 'ATCGAT', 'CpoI': 'CGGWCCG', 'CseI': 'GACGC', 'CsiI': 'ACCWGGT', 'Csp6I': 'GTAC', 'CspCI': 'CAANNNNNGTGG', 'CviAII': 'CATG', 'CviKI-1': 'RGCY', 'CviQI': 'GTAC', 'DdeI': 'CTNAG', 'DpnI': 'GATC', 'DpnII': 'GATC', 'DraI': 'TTTAAA', 'DraIII': 'CACNNNGTG', 'DrdI': 'GACNNNNNNGTC', 'EaeI': 'YGGCCR', 'EagI': 'CGGCCG', 'Eam1104I': 'CTCTTC', 'Eam1105I': 'GACNNNNNGTC', 'EarI': 'CTCTTC', 'EciI': 'GGCGGA', 'Ecl136II': 'GAGCTC', 'Eco105I': 'TACGTA', 'Eco130I': 'CCWWGG', 'Eco147I': 'AGGCCT', 'Eco24I': 'GRGCYC', 'Eco31I': 'GGTCTC', 'Eco32I': 'GATATC', 'Eco47I': 'GGWCC', 'Eco47III': 'AGCGCT', 'Eco52I': 'CGGCCG', 'Eco53kI': 'GAGCTC', 'Eco57I': 'CTGAAG', 'Eco57MI': 'CTGRAG', 'Eco72I': 'CACGTG', 'Eco81I': 'CCTNAGG', 'Eco88I': 'CYCGRG', 'Eco91I': 'GGTNACC', 'EcoNI': 'CCTNNNNNAGG', 'EcoO109I': 'RGGNCCY', 'EcoP15I': 'CAGCAG', 'EcoRI': 'GAATTC', 'EcoRII': 'CCWGG', 'EcoRV': 'GATATC', 'EheI': 'GGCGCC', 'Esp3I': 'CGTCTC', 'FaqI': 'GGGAC', 'FatI': 'CATG', 'FauI': 'CCCGC', 'Fnu4HI': 'GCNGC', 'FokI': 'GGATG', 'FseI': 'GGCCGGCC', 'FspAI': 'RTGCGCAY', 'FspBI': 'CTAG', 'FspI': 'TGCGCA', 'GsuI': 'CTGGAG', 'HaeII': 'RGCGCY', 'HaeIII': 'GGCC', 'HgaI': 'GACGC', 'HhaI': 'GCGC', 'Hin1I': 'GRCGYC', 'Hin1II': 'CATG', 'Hin4I': 'GAYNNNNNVTC', 'Hin6I': 'GCGC', 'HincII': 'GTYRAC', 'HindIII': 'AAGCTT', 'HinfI': 'GANTC', 'HinP1I': 'GCGC', 'HpaI': 'GTTAAC', 'HpaII': 'CCGG', 'HphI': 'GGTGA', 'Hpy166II': 'GTNNAC', 'Hpy188I': 'TCNGA', 'Hpy188III': 'TCNNGA', 'Hpy8I': 'GTNNAC', 'Hpy99I': 'CGWCG', 'HpyAV': 'CCTTC', 'HpyCH4III': 'ACNGT', 'HpyCH4IV': 'ACGT', 'HpyCH4V': 'TGCA', 'HpyF10VI': 'GCNNNNNNNGC', 'HpyF3I': 'CTNAG', 'I-CeuI': 'TAACTATAACGGTCCTAAGGTAGCGA', 'I-SceI': 'TAGGGATAACAGGGTAAT', 'KasI': 'GGCGCC', 'KflI': 'GGGWCCC', 'Kpn2I': 'TCCGGA', 'KpnI': 'GGTACC', 'KspAI': 'GTTAAC', 'LguI': 'GCTCTTC', 'LpnPI': 'CCDG', 'Lsp1109I': 'GCAGC', 'LweI': 'GCATC', 'MauBI': 'CGCGCGCG', 'MbiI': 'GAGCGG', 'MboI': 'GATC', 'MboII': 'GAAGA', 'MfeI': 'CAATTG', 'MlsI': 'TGGCCA', 'MluCI': 'AATT', 'MluI': 'ACGCGT', 'MlyI': 'GAGTC', 'MmeI': 'TCCRAC', 'MnlI': 'CCTC', 'Mph1103I': 'ATGCAT', 'MreI': 'CGCCGGCG', 'MscI': 'TGGCCA', 'MseI': 'TTAA', 'MslI': 'CAYNNNNRTG', 'MspA1I': 'CMGCKG', 'MspI': 'CCGG', 'MssI': 'GTTTAAAC', 'MunI': 'CAATTG', 'Mva1269I': 'GAATGC', 'MvaI': 'CCWGG', 'MwoI': 'GCNNNNNNNGC', 'NaeI': 'GCCGGC', 'NarI': 'GGCGCC', 'NciI': 'CCSGG', 'NcoI': 'CCATGG', 'NdeI': 'CATATG', 'NgoMIV': 'GCCGGC', 'NheI': 'GCTAGC', 'NlaIII': 'CATG', 'NlaIV': 'GGNNCC', 'NmeAIII': 'GCCGAG', 'NmuCI': 'GTSAC', 'NotI': 'GCGGCCGC', 'NruI': 'TCGCGA', 'NsbI': 'TGCGCA', 'NsiI': 'ATGCAT', 'NspI': 'RCATGY', 'OliI': 'CACNNNNGTG', 'PacI': 'TTAATTAA', 'PaeI': 'GCATGC', 'PaeR7I': 'CTCGAG', 'PagI': 'TCATGA', 'PaqCI': 'CACCTGC', 'PasI': 'CCCWGGG', 'PauI': 'GCGCGC', 'PciI': 'ACATGT', 'PdiI': 'GCCGGC', 'PdmI': 'GAANNNNTTC', 'PfeI': 'GAWTC', 'Pfl23II': 'CGTACG', 'PflFI': 'GACNNNGTC', 'PflMI': 'CCANNNNNTGG', 'PfoI': 'TCCNGGA', 'PhoI': 'GGCC', 'PleI': 'GAGTC', 'PluTI': 'GGCGCC', 'PmeI': 'GTTTAAAC', 'PmlI': 'CACGTG', 'PpiI': 'GAACNNNNNCTC', 'Ppu21I': 'YACGTR', 'PpuMI': 'RGGWCCY', 'PscI': 'ACATGT', 'PshAI': 'GACNNNNGTC', 'PsiI': 'TTATAA', 'Psp1406I': 'AACGTT', 'Psp5II': 'RGGWCCY', 'PspFI': 'CCCAGC', 'PspGI': 'CCWGG', 'PspOMI': 'GGGCCC', 'PspXI': 'VCTCGAGB', 'PstI': 'CTGCAG', 'PsuI': 'RGATCY', 'PsyI': 'GACNNNGTC', 'PteI': 'GCGCGC', 'PvuI': 'CGATCG', 'PvuII': 'CAGCTG', 'RruI': 'TCGCGA', 'RsaI': 'GTAC', 'RseI': 'CAYNNNNRTG', 'RsrII': 'CGGWCCG', 'SacI': 'GAGCTC', 'SacII': 'CCGCGG', 'SalI': 'GTCGAC', 'SapI': 'GCTCTTC', 'SaqAI': 'TTAA', 'SatI': 'GCNGC', 'Sau3AI': 'GATC', 'Sau96I': 'GGNCC', 'SbfI': 'CCTGCAGG', 'ScaI': 'AGTACT', 'SchI': 'GAGTC', 'ScrFI': 'CCNGG', 'SdaI': 'CCTGCAGG', 'SduI': 'GDGCHC', 'SexAI': 'ACCWGGT', 'SfaAI': 'GCGATCGC', 'SfaNI': 'GCATC', 'SfcI': 'CTRYAG', 'SfiI': 'GGCCNNNNNGGCC', 'SfoI': 'GGCGCC', 'SgrAI': 'CRCCGGYG', 'SgrDI': 'CGTCGACG', 'SgsI': 'GGCGCGCC', 'SmaI': 'CCCGGG', 'SmiI': 'ATTTAAAT', 'SmlI': 'CTYRAG', 'SmoI': 'CTYRAG', 'SmuI': 'CCCGC', 'SnaBI': 'TACGTA', 'SpeI': 'ACTAGT', 'SphI': 'GCATGC', 'SrfI': 'GCCCGGGC', 'SsiI': 'CCGC', 'SspDI': 'GGCGCC', 'SspI': 'AATATT', 'StuI': 'AGGCCT', 'StyD4I': 'CCNGG', 'StyI': 'CCWWGG', 'SwaI': 'ATTTAAAT', 'TaaI': 'ACNGT', 'TaiI': 'ACGT', 'TaqI': 'TCGA', 'TaqαI': 'TCGA', 'TasI': 'AATT', 'TatI': 'WGTACW', 'TauI': 'GCSGC', 'TfiI': 'GAWTC', 'Tru1I': 'TTAA', 'TscAI': 'CASTG', 'TseI': 'GCWGC', 'TsoI': 'TARCCA', 'Tsp45I': 'GTSAC', 'Tsp509I': 'AATT', 'TspMI': 'CCCGGG', 'TspRI': 'CASTG', 'TstI': 'CACNNNNNNTCC', 'Tth111I': 'GACNNNGTC', 'Van91I': 'CCANNNNNTGG', 'VspI': 'ATTAAT', 'XagI': 'CCTNNNNNAGG', 'XapI': 'RAATTY', 'XbaI': 'TCTAGA', 'XceI': 'RCATGY', 'XcmI': 'CCANNNNNNNNNTGG', 'XhoI': 'CTCGAG', 'XmaI': 'CCCGGG', 'XmaJI': 'CCTAGG', 'XmiI': 'GTMKAC', 'XmnI': 'GAANNNNTTC', 'ZraI': 'GACGTC'}

    # Get hold of any cut sites that were detected and report any sites that were mismatched
    spec_enzymes = {i.lower() for i in enzymes}
    cut_list = {Enz_dict[key] for key in Enz_dict if key.lower() in spec_enzymes}
    missed = spec_enzymes - {key.lower() for key in Enz_dict}
    if len(missed) > 0:
        for item in missed:
            print(item+'not in our enzyme list')

    # provide a list of motifs to avoid in output DNA and combine with enzyme cut sites to avoid
    pattern  = list(['TTTTAGG', 'AAAATCC', 'CCCCCAGG','AATAAA', 'ATTAAA', 'TTTTT', 'TTTTAT', 'AAATAT', 'AAAAA','TATAA', 'ATATAT', 'TATATA','GCCACCATG', 'GCCGCCATG', 'ACCACCATG', 'ACCGCCATG', 'CCCCCC', 'GGGGGG', 'TAAGGAGGT'])
    # collection of banned sequences based on pattern - can expand this
    pattern2  = list(['AATAAA', 'ATTAAA', 'GCCACCATG', 'GCCGCCATG', 'TTTTTAGG', 'TTTTT', 'TATAA', 'ATATAT', 'TATATA', 'AAAAATCC', 'AAAAAA','CCCCCCC', 'GGGGGGG'])
    # less stringent list, in case codon options are too restricted
    for site in cut_list:
        pattern.append(site)
        pattern2.append(site)
        if site != str(Seq(site).reverse_complement()):
            pattern.append(str(Seq(site).reverse_complement()))
            pattern2.append(str(Seq(site).reverse_complement()))

    return pattern, pattern2


"""
Copyright (c) 2024 VECTOR FUTURES LTD
All rights reserved.
This file is part of the Vector Futures Codon Optimization App and may not be copied, distributed, or modified without express written permission.
"""