from flask import Flask, render_template, request
import os
import requests
from Bio.PDB import PDBParser

app = Flask(__name__)

PDB_DIR = "static/pdb"
os.makedirs(PDB_DIR, exist_ok=True)

# -----------------------------
# Download protein
# -----------------------------
def download_pdb(pdb_id):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    path = os.path.join(PDB_DIR, f"{pdb_id}.pdb")
    r = requests.get(url)

    if not r.ok:
        return None

    with open(path, "w") as f:
        f.write(r.text)

    return path

# -----------------------------
# Protein quality evaluation
# -----------------------------
def quality_score(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    gaps = 0
    atoms = 0

    for model in structure:
        for chain in model:
            residues = list(chain)
            for i in range(len(residues) - 1):
                if residues[i+1].id[1] - residues[i].id[1] > 1:
                    gaps += 1
            for r in residues:
                atoms += len(r)

    score = 100 - min(30, gaps * 5)
    if atoms < 1000:
        score -= 20

    return score, gaps, atoms

# -----------------------------
# Routes
# -----------------------------
@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        pdb_id = request.form["pdb_id"].lower()

        pdb_path = download_pdb(pdb_id)
        if not pdb_path:
            return render_template("index.html", error="PDB ID not found")

        score, gaps, atoms = quality_score(pdb_path)

        return render_template(
            "result.html",
            pdb_id=pdb_id.upper(),
            score=score,
            gaps=gaps,
            atoms=atoms,
            pdb_file=f"pdb/{pdb_id}.pdb"
        )

    return render_template("index.html")

if __name__ == "__main__":
    app.run(debug=True)

