# main.py

from flask import Flask, render_template, request, redirect, url_for
import os
import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import re
import textwrap
import numpy as np
from io import StringIO, BytesIO
from socket import inet_aton
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from collections import defaultdict
from tqdm import tqdm
from app import codon_optimization, enzyme_utils, progress_utils, sequence_utils, pattern_generator



app = Flask(__name__)

# Define a dictionary of codon tables and their corresponding CSV file paths
codon_tables = {
    'Homo sapiens': 'data/Homo.csv',
    'Mus musculus': 'data/Mus_musculus.csv',
    'Drosophila melanogaster': 'data/Drosophila_melanogaster.csv',
    'Glycine max': 'data/Glycine_max.csv',
    'Staph aureus': 'data/Staph_aureus.csv',
    'Chlamydomonas reinhardtii': 'data/Chlamydomonas_reinhardtii.csv',
    'AAV2': 'data/AAV2.csv',
    'HIV-1': 'data/HIV-1.csv'
}

# Define GC optimization levels in a dictionary with predefined ranges
optimization_levels = {
    'Low (40-50% GC)': (45, 55),
    'Low medium (45-55% GC)': (55, 65),
    'High medium (50-60% GC)': (60, 73),
    'High (55-66% GC)': (65, 77)
}

# Allowed characters for sequence input (Amino acid single-letter codes, *, X)
allowed_characters = set("ARNDCQEGHILKMFPSVTWYX*")

@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        try:
            # Retrieve and sanitize the sequence input
            sequence_input = request.form.get("sequence_input").strip().upper()
            sanitized_sequence = ''.join([char for char in sequence_input if char in allowed_characters])

            # Check if sanitized sequence is valid (not empty and meets length criteria)
            if not sanitized_sequence:
                flash("Error: The submitted sequence contains invalid characters or is empty. Please submit a valid DNA or protein sequence.")
                return render_template("submission_form.html")  # Render the submission page again

            # Get the selected codon table
            selected_codon_table = request.form.get("codon_table")
            codon_df = pd.read_csv(codon_tables[selected_codon_table])

            # Handle whether advanced settings are enabled
            use_advanced_settings = request.form.get("use_advanced_settings") == "on"
            
            if use_advanced_settings:
                # Get GC optimization level as a tuple
                selected_gc_opt = optimization_levels[request.form.get("gc_content")]

                # Get CpG depletion level and codon bias or GC selection
                cpg_depletion_level = request.form.get("cpg_depletion")
                codon_bias_or_gc = request.form.get("codon_bias_or_gc")

                # Get enzymes input (restriction sites) and generate patterns
                enzyme_input = request.form.get("enzymes")
                enzymes = re.split(',| ', enzyme_input) if enzyme_input else []  # Split enzyme input or leave empty

                # Generate patterns using pattern_generator
                pattern, pattern2 = pattern_generator(enzymes)

                # Check if the user wants to expedite processing
                expedition = request.form.get("expedition") == "Yes"
            else:
                # Default settings
                enzymes = []
                selected_gc_opt = (63, 77)
                cpg_depletion_level = "None"
                codon_bias_or_gc = "GC richness"
                pattern, pattern2 = pattern_generator(enzymes)
                expedition = True

            # Call codon optimization function
            optimized_seq = codon_optimization(
                uploaded_seq_file=sanitized_sequence,
                cpg_depletion_level=cpg_depletion_level,
                codon_bias_or_GC=codon_bias_or_gc,
                pattern=pattern,
                pattern2=pattern2,
                codon_df=codon_df,
                selected_GC_opt=selected_gc_opt,
                expedite=expedition
            )
            
            return render_template("result.html", optimized_seq=optimized_seq)


        except Exception as e:
            # Log the error for debugging purposes (optional)
            print(f"Unexpected error occurred: {str(e)}")


    # Render the form on GET request
    return render_template("index.html", codon_tables=codon_tables.keys(), optimization_levels=optimization_levels.keys())

# Result page to display the optimized sequence
@app.route("/result")
def result():
    return render_template("result.html")

if __name__ == "__main__":
    app.run(debug=True)



"""
Copyright (c) 2024 VECTOR FUTURES LTD
All rights reserved.
This file is part of the Vector Futures Codon Optimization App and may not be copied, distributed, or modified without express written permission.
"""