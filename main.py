from flask import Flask, render_template, request, flash
from Bio.Seq import Seq
from itertools import product
from typing import List, Tuple
import re
import random
from app.motif_optimizer import *

app = Flask(__name__)
app.secret_key = 'your_secret_key'  # Required for flashing messages


# Allowed characters for sequence input
allowed_characters = set("ATUGC*")

# Define options
IUPAC_options = {
    'Insert': "Insert",
    'Remove': "Remove"
}

@app.route("/", methods=["GET", "POST"])
def index():
    optimized_seq = None  # Initialize the result as None
    sanitized_sequence = None
    exact_matches = 0
    new_exact_matches = 0
    if request.method == "POST":
        try:
            # Retrieve and sanitize the sequence input
            sequence_input = request.form.get("sequence_input", "").strip().upper()
            modified_sequence_input = sequence_input.replace("U", "T").replace("X", "*")
            sanitized_sequence = ''.join([char for char in modified_sequence_input if char in allowed_characters])

            # Check if sanitized sequence is valid (not empty and meets length criteria)
            if not sanitized_sequence:
                flash("Error: The submitted sequence contains invalid characters or is empty. Please submit a valid DNA sequence.")
                return render_template("index.html", optimized_seq=None, exact_matches=None, new_exact_matches=None)

            # Choose whether to insert or delete motifs
            selected_optimization = IUPAC_options[request.form.get("iupac_option")]
            
            # Get the iupac motif
            iupac_motifs = request.form.get("iupac", "")
            iupac_list = re.split(r',|\s+', iupac_motifs.strip()) if iupac_motifs else []
            iupac_list = list(filter(None, iupac_list))
            print(iupac_list)

            # Call the motif optimization function
            if selected_optimization == "Insert":
                for iupac_motif in iupac_list:
                    x, y, z = insert_motif(
                        dna_sequence=sanitized_sequence,
                        iupac=iupac_motif.upper())
                    sanitized_sequence = x
                    exact_matches += y
                for iupac_motif in iupac_list:
                    c = count_motif(
                        dna_sequence=sanitized_sequence,
                        iupac=iupac_motif.upper())
                    new_exact_matches += c
                    
            elif selected_optimization == "Remove":
                for iupac_motif in iupac_list:
                    x, y, z = destroy_motif(
                        dna_sequence=sanitized_sequence,
                        iupac=iupac_motif.upper())
                    sanitized_sequence = x
                    exact_matches += y
                for iupac_motif in iupac_list:
                    c = count_motif(
                        dna_sequence=sanitized_sequence,
                        iupac=iupac_motif.upper())
                    new_exact_matches += c

        except Exception as e:
            flash(f"Unexpected error occurred: {str(e)}")
            return render_template("index.html", optimized_seq=None, exact_matches=None, new_exact_matches=None)

    # Render the form and (if available) the optimized sequence
    return render_template("index.html", optimized_seq=sanitized_sequence, exact_matches=exact_matches, new_exact_matches=new_exact_matches)

if __name__ == "__main__":
    app.run(debug=True)


"""
Copyright (c) 2024 VECTOR FUTURES LTD
All rights reserved.
This file is part of the Vector Futures Codon Optimization App and may not be copied, distributed, or modified without express written permission.
"""