import os
import requests
import MDAnalysis as mda
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from Bio.PDB import PDBParser
import py3Dmol

# -----------------------------
# 1. Download protein from RCSB
# -----------------------------
def download_protein(pdb_id):
    os.makedirs("protein_structures", exist_ok=True)
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    path = f"protein_structures/{pdb_id}.pdb"

    r = requests.get(url)
    if not r.ok:
        raise ValueError(f"❌ PDB ID {pdb_id} not found in RCSB")

    with open(path, "w") as f:
        f.write(r.text)

    print(f"✅ Downloaded {path}")
    return path

# -----------------------------
# 2. Fix protein with PDBFixer
# -----------------------------
def fix_protein(pdb_path):
    output = pdb_path.replace(".pdb", "_fixed.pdb")

    fixer = PDBFixer(filename=pdb_path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)
    fixer.removeHeterogens(keepWater=False)

    with open(output, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

    print(f"✅ Fixed protein saved: {output}")
    return output

# -----------------------------
# 3. Convert PDB → PDBQT
# -----------------------------
def pdb_to_pdbqt(pdb_path):
    pdbqt = pdb_path.replace(".pdb", ".pdbqt")
    u = mda.Universe(pdb_path)
    u.atoms.write(pdbqt)
    print(f"✅ PDBQT created: {pdbqt}")
    return pdbqt

# -----------------------------
# 4. Protein quality score
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

    return max(score, 0), gaps, atoms

# -----------------------------
# 5. Visualization (py3Dmol)
# -----------------------------
def visualize_structure(pdb_file):
    print("🧬 Opening 3D visualization...")
    with open(pdb_file) as f:
        view = py3Dmol.view(width=800, height=600)
        view.addModel(f.read(), "pdb")
        view.setStyle({"cartoon": {"color": "spectrum"}})
        view.zoomTo()
        view.show()

# -----------------------------
# MAIN
# -----------------------------
if __name__ == "__main__":

    pdb_id = input("🔬 Enter PDB ID (example: 5k5x): ").strip().lower()

    try:
        print("\n🔽 Downloading protein...")
        raw_pdb = download_protein(pdb_id)

        print("🛠 Fixing protein...")
        fixed_pdb = fix_protein(raw_pdb)

        print("🔁 Converting to PDBQT...")
        pdbqt = pdb_to_pdbqt(fixed_pdb)

        score, gaps, atoms = quality_score(fixed_pdb)

        print("\n=== PROTEIN QUALITY REPORT ===")
        print(f"PDB ID          : {pdb_id.upper()}")
        print(f"Quality score  : {score}/100")
        print(f"Missing gaps   : {gaps}")
        print(f"Atom count     : {atoms}")
        print("==============================")

        if score < 70:
            print("⚠️ Protein NOT recommended for docking")
        else:
            print("✅ Protein READY for docking")

        # ---- Visualization ----
        visualize_structure(fixed_pdb)

    except Exception as e:
        print(f"\n❌ ERROR: {e}")

