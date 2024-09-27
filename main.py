from flask import Flask, render_template, request, flash
import pandas as pd
import re
from app import codon_optimization, pattern_generator

app = Flask(__name__)
app.secret_key = 'vector-barcelona-rain-tehran-obvious-shenanigans'  # Required for flashing messages

# Define a dictionary of codon tables and their corresponding CSV file paths
codon_tables = {
    'Homo sapiens': 'data/homo.csv',
    'Mus musculus': 'data/Mus_musculus.csv',
    'Drosophila melanogaster': 'data/Drosophila_melanogaster.csv',
    'Glycine max': 'data/Glycine_max.csv',
    'Staph aureus': 'data/Staph_aureus.csv',
    'Chlamydomonas reinhardtii': 'data/Chlamydomonas_reinhardtii.csv',
    'AAV2': 'data/AAV2.csv',
    'HIV-1': 'data/HIV-1.csv'
}

# Define GC optimization levels
optimization_levels = {
    'Low (40-50% GC)': (45, 55),
    'Low medium (45-55% GC)': (55, 65),
    'High medium (50-60% GC)': (60, 73),
    'High (55-66% GC)': (64, 77),
    'Very high (60-70% GC)': (70, 80)
}

# Allowed characters for sequence input (Amino acid single-letter codes, *, X)
allowed_characters = set("ARNDCQEGHILKMFPSVTWYX*")

@app.route("/", methods=["GET", "POST"])
def index():
    optimized_seq = None  # Initialize the result as None
    if request.method == "POST":
        try:
            # Retrieve and sanitize the sequence input
            sequence_input = request.form.get("sequence_input").strip().upper()
            sanitized_sequence = ''.join([char for char in sequence_input if char in allowed_characters])

            # Check if sanitized sequence is valid (not empty and meets length criteria)
            if not sanitized_sequence:
                flash("Error: The submitted sequence contains invalid characters or is empty. Please submit a valid DNA or protein sequence.")
                return render_template("index.html", codon_tables=codon_tables.keys(), optimization_levels=optimization_levels.keys(), optimized_seq=None)

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
                # Get enzymes input and generate patterns
                enzyme_input = request.form.get("enzymes")
                enzymes = re.split(',| ', enzyme_input) if enzyme_input else []  # Split enzyme input or leave empty
                pattern, pattern2 = pattern_generator(enzymes)
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

        except Exception as e:
            flash(f"Unexpected error occurred: {str(e)}")
            return render_template("index.html", codon_tables=codon_tables.keys(), optimization_levels=optimization_levels.keys(), optimized_seq=None)

    # Render the form and (if available) the optimized sequence
    return render_template("index.html", codon_tables=codon_tables.keys(), optimization_levels=optimization_levels.keys(), optimized_seq=optimized_seq)

if __name__ == "__main__":
    app.run(debug=True)
