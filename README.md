Overview: This is a Flask-based web application for controlling the presence or absence of specific motifs in a synthetic DNA coding sequence. The app allows users to search and destroy IUPAC motifs, or force IUPAC motifs into their DNA sequence, without altering the encoded protein sequence.

Features
Sequence Input: Enter DNA or mRNA sequence as raw text.
Select whether to insert or destroy motifs.
Define the motif to insert or destroy, written in IUPAC notation.
Result Display: View the modified sequence with an option to copy the result to your clipboard.
Requirements
Python 3.x
Flask
Gunicorn (for production)
pandas, BioPython, numpy
Installation
Clone the repository:

bash
Copy code
git clone https://github.com/vf-collab/vf-motif_control.git
cd vf-motif-control
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

