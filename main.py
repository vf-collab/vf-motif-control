from flask import Flask, render_template, request, redirect, url_for, jsonify
from celery import Celery
import pandas as pd
import re
import numpy as np
from app import codon_optimization, enzyme_utils, progress_utils, sequence_utils, pattern_generator

# Initialize Flask app
app = Flask(__name__)

# Celery configuration
app.config.update(
    CELERY_BROKER_URL='redis://209.97.136.147:6379/0',
    CELERY_RESULT_BACKEND='redis://209.97.136.147:6379/0'
)

# Initialize Celery
def make_celery(app):
    celery = Celery(
        app.import_name,
        backend=app.config['CELERY_RESULT_BACKEND'],
        broker=app.config['CELERY_BROKER_URL']
    )
    celery.conf.update(app.config)
    return celery

celery = make_celery(app)

# Codon tables and settings
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

optimization_levels = {
    'Low (40-50% GC)': (45, 55),
    'Low medium (45-55% GC)': (55, 65),
    'High medium (50-60% GC)': (60, 73),
    'High (55-66% GC)': (65, 77),
    'Very high (60-70% GC)': (73, 83)
}

allowed_characters = set("ARNDCQEGHILKMFPSVTWYX*")

# Task for background codon optimization
@celery.task
def run_codon_optimization(sanitized_sequence, selected_codon_table, use_advanced_settings, params):
    # Retrieve codon table data
    codon_df = pd.read_csv(codon_tables[selected_codon_table])

    if use_advanced_settings:
        selected_gc_opt = params['selected_gc_opt']
        cpg_depletion_level = params['cpg_depletion_level']
        codon_bias_or_gc = params['codon_bias_or_gc']
        pattern = params['pattern']
        pattern2 = params['pattern2']
        expedition = params['expedition']
    else:
        selected_gc_opt = (63, 77)
        cpg_depletion_level = "None"
        codon_bias_or_gc = "GC richness"
        pattern = params['pattern']
        pattern2 = params['pattern2']
        expedition = True

    # Call the codon optimization function
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
    return optimized_seq

# Flask route to submit form and process codon optimization
@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        try:
            # Retrieve and sanitize the sequence input
            sequence_input = request.form.get("sequence_input").strip().upper()
            sanitized_sequence = ''.join([char for char in sequence_input if char in allowed_characters])

            if not sanitized_sequence:
                return jsonify({"error": "Invalid sequence"}), 400

            # Get form data
            selected_codon_table = request.form.get("codon_table")
            use_advanced_settings = request.form.get("use_advanced_settings") == "on"

            params = {}
            if use_advanced_settings:
                params['selected_gc_opt'] = optimization_levels[request.form.get("gc_content")]
                params['cpg_depletion_level'] = request.form.get("cpg_depletion")
                params['codon_bias_or_gc'] = request.form.get("codon_bias_or_gc")
                enzyme_input = request.form.get("enzymes")
                params['pattern'], params['pattern2'] = pattern_generator(re.split(',| ', enzyme_input) if enzyme_input else [])
                params['expedition'] = request.form.get("expedition") == "Yes"
            else:
                params['pattern'], params['pattern2'] = pattern_generator([])

            # Start the background task for codon optimization
            task = run_codon_optimization.delay(sanitized_sequence, selected_codon_table, use_advanced_settings, params)

            # Return task ID to allow the client to poll for results
            return jsonify({"task_id": task.id})

        except Exception as e:
            return jsonify({"error": str(e)}), 500

    return render_template("index.html", codon_tables=codon_tables.keys(), optimization_levels=optimization_levels.keys())

# Route to check the task status
@app.route('/task_status/<task_id>')
def task_status(task_id):
    task = run_codon_optimization.AsyncResult(task_id)
    if task.state == 'PENDING':
        response = {'state': task.state, 'progress': 0}
    elif task.state == 'PROGRESS':
        response = {'state': task.state, 'progress': task.info.get('progress', 0)}
    elif task.state == 'SUCCESS':
        response = {'state': task.state, 'result': task.result}
    else:
        response = {'state': task.state, 'result': str(task.info)}
    return jsonify(response)

if __name__ == "__main__":
    app.run(debug=True)