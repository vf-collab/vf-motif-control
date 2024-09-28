from flask import Flask, render_template, request, jsonify, Response, stream_with_context
import time
import pandas as pd
from threading import Thread
from app import codon_optimization, pattern_generator
import uuid

app = Flask(__name__)
app.secret_key = 'vector-barcelona-rain-tehran-obvious-eixample'  # Required for flashing messages

# A dictionary to store progress of each task
progress = {}

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
    'Low (40-45% GC)': (45, 55),
    'Low medium (45-50% GC)': (55, 65),
    'High medium (50-55% GC)': (60, 73),
    'High (55-60% GC)': (64, 77),
    'Very high (60-65% GC)': (75, 85)
}

# Allowed characters for sequence input (Amino acid single-letter codes, *, X)
allowed_characters = set("ARNDCQEGHILKMFPSVTWYX*")


def run_optimization(task_id, sanitized_sequence, selected_gc_opt, cpg_depletion_level, codon_bias_or_gc, pattern, pattern2, codon_df, expedition):
    """
    The function that performs codon optimization in a separate thread and updates progress.
    """
    total_steps = 100  # Assuming 100 steps for simplicity, adjust this dynamically based on your real task
    progress[task_id] = 0  # Initialize progress

    try:
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

        # Simulate progress as optimization proceeds
        for step in range(1, total_steps + 1):
            time.sleep(0.1)  # Simulate a step (you should update this with real task steps)
            progress[task_id] = step  # Update global progress
        progress[task_id] = 100  # Mark the task as completed

    except Exception as e:
        progress[task_id] = -1  # Error occurred, mark with -1
        print(f"Error in codon optimization: {e}")


import logging
logging.basicConfig(level=logging.DEBUG)

@app.route("/", methods=["GET", "POST"])
def index():
    optimized_seq = None
    if request.method == "POST":
        try:
            logging.debug('Form submission started')
            # Retrieve and sanitize the sequence input
            form_data = request.get_json()  # This captures JSON data from the request body
            sequence_input = form_data.get("sequence_input").strip().upper()
            sanitized_sequence = ''.join([char for char in sequence_input if char in allowed_characters])

            logging.debug('Sanitized Sequence: %s', sanitized_sequence)

            # Check if sanitized sequence is valid
            if not sanitized_sequence:
                return jsonify({"error": "Invalid sequence"}), 400

            # Get the selected codon table
            selected_codon_table = form_data.get("codon_table")
            codon_df = pd.read_csv(codon_tables[selected_codon_table])

            # Handle advanced settings
            use_advanced_settings = form_data.get("use_advanced_settings") == "on"
            if use_advanced_settings:
                selected_gc_opt = optimization_levels[form_data.get("gc_content")]
                cpg_depletion_level = form_data.get("cpg_depletion")
                codon_bias_or_gc = form_data.get("codon_bias_or_gc")
                enzyme_input = form_data.get("enzymes")
                enzymes = re.split(',| ', enzyme_input) if enzyme_input else []
                pattern, pattern2 = pattern_generator(enzymes)
                expedition = form_data.get("expedition") == "Yes"
            else:
                enzymes = []
                selected_gc_opt = (63, 77)
                cpg_depletion_level = "None"
                codon_bias_or_gc = "GC richness"
                pattern, pattern2 = pattern_generator(enzymes)
                expedition = True

            logging.debug('Calling codon optimization')
            # Trigger the background task
            task = run_codon_optimization.delay(sanitized_sequence, selected_codon_table, use_advanced_settings, {
                'selected_gc_opt': selected_gc_opt,
                'cpg_depletion_level': cpg_depletion_level,
                'codon_bias_or_gc': codon_bias_or_gc,
                'pattern': pattern,
                'pattern2': pattern2,
                'expedition': expedition
            })

            logging.debug('Task started, task id: %s', task.id)
            # Return task ID for progress tracking
            return jsonify({"task_id": task.id})

        except Exception as e:
            logging.exception('Error occurred during form submission')
            return jsonify({"error": str(e)}), 500

    return render_template("index.html", codon_tables=codon_tables.keys(), optimization_levels=optimization_levels.keys())


@app.route("/progress/<task_id>")
def progress_stream(task_id):
    """
    Server-Sent Events (SSE) to provide real-time progress updates.
    """
    def generate():
        while True:
            task_progress = progress.get(task_id, 0)
            yield f"data: {task_progress}\n\n"
            if task_progress >= 100 or task_progress == -1:  # Stop when the task is done or an error occurs
                break
            time.sleep(1)
    return Response(stream_with_context(generate()), mimetype="text/event-stream")


if __name__ == "__main__":
    app.run(debug=True)
