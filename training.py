import os
import numpy as np
from Bio.PDB import PDBParser
from scipy.spatial.distance import pdist, squareform
from pathlib import Path

# -------------------------------------------------------------
# Configuration
# -------------------------------------------------------------
RNA_NT = {"A", "U", "G", "C"}
TARGET_ATOM = "C3'"
NUM_BINS = 20
BIN_WIDTH = 1.0
MAX_DIST = 20.0
MIN_SEQ_SEP = 4  # i and i+4

# Update these paths
BASE_PATH = input(Path(" copy your folder path Input"))
pdb=input(str(" copy your pdb file name Input"))
PDB_FILE = BASE_PATH /pdb
ENERGY_FOLDER = BASE_PATH /os.makedirs("training_sets")


# -------------------------------------------------------------
# 1. Load Energy Profiles (Symmetric)
# -------------------------------------------------------------
def load_energy_profiles(folder):
    """
    Loads energy profiles into a dictionary. 
    Ensures keys are sorted tuples so (A, U) is the same as (U, A).
    """
    energy = {}
    folder = Path(folder)
    
    # We loop through unique combinations to avoid double checking
    for n1 in RNA_NT:
        for n2 in RNA_NT:
            # Try to find file for n1n2
            filename = folder / f"{n1}{n2}.txt"
            
            if filename.exists():
                with open(filename, 'r') as f:
                    # Convert to numpy array for fast access
                    values = np.array([float(line.strip()) for line in f.readlines()])
                    
                # Store symmetrically: sorted tuple key ensures (A,G) == (G,A)
                key = tuple(sorted((n1, n2)))
                energy[key] = values
                
    print(f"Loaded {len(energy)} energy profiles.")
    return energy

# -------------------------------------------------------------
# 2. Extract Data (Coordinates & Sequence)
# -------------------------------------------------------------
def extract_structure_data(pdb_file, chain_id="A"):
    """
    Parses PDB and returns:
    - coords: Numpy array (N, 3)
    - sequence: List of nucleotide names corresponding to coords
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("RNA", pdb_file)
    
    try:
        chain = structure[0][chain_id]
    except KeyError:
        raise ValueError(f"Chain {chain_id} not found in {pdb_file}")

    coords = []
    sequence = []

    for residue in chain.get_residues():
        res_name = residue.get_resname().strip()
        
        # Skip non-RNA residues
        if res_name not in RNA_NT:
            continue

        # Find C3' atom
        if TARGET_ATOM in residue:
            atom = residue[TARGET_ATOM]
            coords.append(atom.get_coord())
            sequence.append(res_name)
            
    return np.array(coords), sequence

# -------------------------------------------------------------
# 3. Optimized Scoring Calculation
# -------------------------------------------------------------
def interpolated_score_fast(dist, profile):
    """
    Vectorized or Scalar linear interpolation.
    """
    # Calculate lower bin index
    bin_low = np.floor(dist / BIN_WIDTH).astype(int)
    
    # Check bounds (if dist is exactly MAX_DIST, handle carefully)
    if bin_low >= NUM_BINS - 1:
        return profile[-1] # Return last bin value
        
    r_low = bin_low * BIN_WIDTH
    
    # Calculate weight
    w = (dist - r_low) / BIN_WIDTH
    
    # Interpolate: (1-w)*low + w*high
    return (1 - w) * profile[bin_low] + w * profile[bin_low + 1]


def score_structure(pdb_file, energy_profiles, chain="A"):
    # 1. Get Data
    coords, sequence = extract_structure_data(pdb_file, chain)
    n_atoms = len(sequence)
    
    if n_atoms < 2:
        return 0.0

    # 2. Compute Distance Matrix (Fastest method)
    # pdist returns a condensed distance matrix (1D array)
    dists_condensed = pdist(coords)
    
    # Convert to square matrix (N x N) to easily map to sequence indices
    dist_matrix = squareform(dists_condensed)

    total_score = 0.0

    # 3. Iterate through upper triangle of the matrix
    # We use range loops here because we need to look up the nucleotide pair types
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            
            # Filter 1: Sequence Separation (using array index)
            if (j - i) < MIN_SEQ_SEP:
                continue

            dist = dist_matrix[i, j]

            # Filter 2: Max Distance cut-off
            if dist >= MAX_DIST:
                continue

            # Identify Pair Type
            nt1 = sequence[i]
            nt2 = sequence[j]
            
            # Sort pair to match our loading logic
            pair_key = tuple(sorted((nt1, nt2)))

            if pair_key in energy_profiles:
                profile = energy_profiles[pair_key]
                score = interpolated_score_fast(dist, profile)
                total_score += score

    return total_score

# -------------------------------------------------------------
# MAIN
# -------------------------------------------------------------
if __name__ == "__main__":
    if not PDB_FILE.exists():
        print(f"Error: PDB file not found at {PDB_FILE}")
    elif not ENERGY_FOLDER.exists():
        print(f"Error: Energy folder not found at {ENERGY_FOLDER}")
    else:
        # Load profiles once
        profiles = load_energy_profiles(ENERGY_FOLDER)
        
        # Run Score
        try:
            gibbs_score = score_structure(PDB_FILE, profiles, chain="A")
            print("-" * 30)
            print(f"Structure: {PDB_FILE.name}")
            print(f"Total Interpolated Score: {gibbs_score:.4f}")
            print("-" * 30)
        except Exception as e:
            print(f"Calculation failed: {e}")