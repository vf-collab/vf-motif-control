Overview
This is a Flask-based web application for optimizing DNA or protein sequences. The app allows users to input sequences and provides various settings for optimizing the sequence, such as GC content, CpG depletion levels, and codon bias.

Features
Sequence Input: Enter DNA or protein sequences as raw text.
Codon Table Selection: Select from a variety of codon tables.
Advanced Settings: Customize optimization settings, including:
Target GC richness.
CpG depletion level.
Codon bias or GC richness optimization.
Enzyme restriction sites exclusion.
Result Display: View the optimized sequence with an option to copy the result to your clipboard.
Requirements
Python 3.x
Flask
Gunicorn (for production)
pandas, BioPython, numpy, tqdm
Installation
Clone the repository:

bash
Copy code
git clone https://github.com/vf-collab/vf-app-base.git
cd vf-app-base
Create and activate a virtual environment:

bash
Copy code
python3 -m venv venv
source venv/bin/activate
Install the dependencies:

bash
Copy code
pip install -r requirements.txt
Run the app locally:

bash
Copy code
python main.py
Deployment
This application can be deployed using various platforms such as DigitalOcean, Heroku, or Vercel. For DigitalOcean, you can use a Droplet or App Platform to deploy.

For production, you can use Gunicorn to run the app.

