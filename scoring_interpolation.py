import itertools
import math
from collections import defaultdict
from Bio.PDB import PDBParser
import os
import sys

# --- CONFIGURATION ---
# Update this path to where your files actually are
PATH = input(str(" copy your scoring  Input path:  "))
PDB_FILENAME = input(str("copy your pdb file name  : ")) # Ensure this file exists in PATH

RNA_NT = ["A", "U", "G", "C"]
TARGET_ATOM = "C3'"
NUM_BINS = 20
BIN_WIDTH = 1.0  # Angstrom
MAX_DIST = NUM_BINS * BIN_WIDTH

def get_unique_pair(r1, r2):
    """
    Sorts residues to ensure symmetry.
    (G, A) and (A, G) both become ('A', 'G').
    """
    return tuple(sorted((r1, r2)))

def distance(a, b):
    diff = a.coord - b.coord
    return float((diff * diff).sum() ** 0.5)

def get_bin_index(dist):
    """Return distance bin index (0–19). Returns None if out of range."""
    if dist < 0 or dist >= MAX_DIST:
        return None
    return int(dist // BIN_WIDTH)

def process_structure(filepath, chain_id="A"):
    """
    Parses PDB and returns a list of (res1, res2, distance).
    Only for C3' atoms separated by sequence >= 3.
    """
    if not os.path.exists(filepath):
        print(f"Error: File {filepath} not found.")
        return []

    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("RNA", filepath)
    except Exception as e:
        print(f"Error parsing PDB: {e}")
        return []

    # Find the chain
    atoms = []
    found_chain = False
    
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                found_chain = True
                for residue in chain.get_residues():
                    resname = residue.get_resname().strip()
                    if resname not in RNA_NT:
                        continue
                    if TARGET_ATOM in residue:
                        atoms.append((residue[TARGET_ATOM], resname, residue.id[1]))
                break # Stop after finding chain
        if found_chain: break

    if not atoms:
        print(f"No atoms found for Chain {chain_id}")
        return []

    # Compute distances
    pair_data = []
    for (a1, name1, idx1), (a2, name2, idx2) in itertools.combinations(atoms, 2):
        if abs(idx1 - idx2) < 3:
            continue

        d = distance(a1, a2)
        
        # Ensure we treat A-G same as G-A immediately
        unique_pair = get_unique_pair(name1, name2)
        pair_data.append((unique_pair[0], unique_pair[1], d))

    return pair_data

def compute_potentials(pair_data):
    """
    Computes the discrete energy table (u) for 20 bins.
    Returns: dictionary u where u[(R1, R2)] = [20 floats]
    """
    
    # Initialize all 10 unique pairs with 0 counts to ensure file creation
    unique_pairs = list(itertools.combinations_with_replacement(sorted(RNA_NT), 2))
    
    N_obs = {pair: [0] * NUM_BINS for pair in unique_pairs}
    N_ref = [0] * NUM_BINS
    N_obs_total = defaultdict(int)
    N_ref_total = 0

    for r1, r2, dist in pair_data:
        bin_idx = get_bin_index(dist)
        if bin_idx is None:
            continue
        
        pair = (r1, r2) # Already sorted in process_structure
        
        N_obs[pair][bin_idx] += 1
        N_obs_total[pair] += 1
        
        N_ref[bin_idx] += 1
        N_ref_total += 1

    # Frequencies and Energy
    f_ref = [0.0] * NUM_BINS
    if N_ref_total > 0:
        f_ref = [x / N_ref_total for x in N_ref]

    u = {}
    for pair in unique_pairs:
        u[pair] = []
        pair_count = N_obs_total[pair]
        
        for r in range(NUM_BINS):
            val = 10.0 # Default high energy
            if pair_count > 0:
                f_obs = N_obs[pair][r] / pair_count
                if f_obs > 0 and f_ref[r] > 0:
                    val = -math.log(f_obs / f_ref[r])
            
            u[pair].append(val)
            
    return u

def get_interpolated_score(dist, energy_bins):
    """
    Linear Interpolation.
    energy_bins: list of 20 values.
    dist: float distance.
    """
    # 1. Map distance to float index
    # If bin width is 1.0, 5.5A is index 5.5
    float_idx = dist / BIN_WIDTH
    
    # 2. Determine lower and upper integer indices
    idx_floor = int(math.floor(float_idx))
    idx_ceil = idx_floor + 1
    
    # 3. Handle out of bounds
    if idx_floor < 0: return 10.0
    if idx_floor >= NUM_BINS: return 0.0 # Distance too large usually means 0 interaction or ignore
    if idx_ceil >= NUM_BINS: return energy_bins[idx_floor] # Cap at last bin

    # 4. Get values
    val_lower = energy_bins[idx_floor]
    val_upper = energy_bins[idx_ceil]
    
    # 5. Calculate fraction
    fraction = float_idx - idx_floor
    
    # 6. Interpolate: y = y1 + (y2 - y1) * fraction
    interpolated_val = val_lower + (val_upper - val_lower) * fraction
    
    return interpolated_val

def calculate_total_score(pair_data, potential_table):
    """
    Sums up interpolated scores for all pairs in the structure.
    """
    total_score = 0.0
    count = 0
    
    for r1, r2, dist in pair_data:
        pair = (r1, r2)
        if pair in potential_table and dist < MAX_DIST:
            score = get_interpolated_score(dist, potential_table[pair])
            total_score += score
            count += 1
            
    return total_score, count

def save_energies(u):
    """Save the 10 files to PATH."""
    if not os.path.exists(PATH):
        try:
            os.makedirs(PATH)
        except OSError:
            print(f"Warning: Could not create directory {PATH}")

    for (r1, r2), scores in u.items():
        filename = os.path.join(PATH, f"{r1}{r2}.txt")
        with open(filename, "w") as f:
            for val in scores:
                f.write(f"{val:.4f}\n")
    print(f"Saved 10 energy files to {PATH}")

# --- MAIN SCRIPT ---
if __name__ == "__main__":
    full_pdb_path = os.path.join(PATH, PDB_FILENAME)
    
    print(f"Processing {full_pdb_path}...")

    # 1. Extract distances (list of tuples: r1, r2, dist)
    pair_data = process_structure(full_pdb_path, chain_id="A")
    print(f"Extracted {len(pair_data)} valid pairs.")

    if len(pair_data) > 0:
        # 2. Compute the discrete potentials (Training step)
        u_table = compute_potentials(pair_data)

        # 3. Save the discrete energies to .txt files
        save_energies(u_table)

        # 4. Calculate Total Score using Linear Interpolation (Scoring step)
        total_energy, num_scored = calculate_total_score(pair_data, u_table)
        
        print("-" * 30)
        print(f"Total Gibs free energy scores : {total_energy:.4f}")
        print(f"Number of pairs scored:   {num_scored}")
        print("-" * 30)
    else:
        print("No pairs found. Check PDB path, Chain ID, or atom names.")