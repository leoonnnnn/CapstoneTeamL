########################################################################################
### This script generates data for testing your TSP/VRP approach from real PDB file. ###
### The only input you need is a pdb file: https://www.rcsb.org/                     ###
########################################################################################

# The idea is that we can simplify sugars and phosphates into single nodes (a sugar node and a phosphate node),
# instead of dealing with all the atoms that compose sugars and phosphates. This makes our TSP/VRP problem easier.
# To find each node, we'll just find the centers of sugars and phosphates.
# Finding the center of phosphates is fairly easy; just take the P atom.
# Finding the center of sugars is a little harder, but can be estimated by finding the center of the sugar ring.
# The ring part of the sugar is composed of O4', C1', C2', C3', and C4'.
# Then we can find the center (technically the centroid) by averaging the positions of those 5 atoms.


# For ease of access,
# Specify input/output locations HERE:
#####################################
# -=-REQUIRED-=-
input_pdb_file = r"C:\Users\noelu\Downloads\6gtg.pdb"   # r for raw string, helps deal with spaces

# -=-Optional-=-
# Leaving this blank as r"" is OK. It'll default to the current folder. 
# If you do specify a folder, make sure it exists locally first, since I didn't set permissions in the code to make new directories.
output_folder = r"C:\Users\noelu\Downloads\DELETEME"   
#####################################


# Helpful resources:
#   PDB file format (go to page 2): https://www.biostat.jhsph.edu/~iruczins/teaching/260.655/links/pdbformat.pdf
#   nucleotide atom numbering: https://www.biosyn.com/tew/numbering-convention-for-nucleotides.aspx
#   Biopython documentation: https://biopython.org/


from Bio import PDB
import numpy as np
import os


# This is a helper function for find_centers().
# This function calculates the centroid (i.e. geometric center) of a list of atom coordinates.
def calculate_centroid(atoms):
    coords = [atom.get_coord() for atom in atoms]
    centroid = np.mean(coords, axis=0)  # Calculate mean of x, y, z
    return centroid


# This function finds the center of sugar rings and phosphates for each RNA residue (nucleotide) in a PDB file.
def find_centers(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    
    results = []

    for model in structure:
        for chain in model:
            for residue in chain:
                sugar_atoms = []
                phosphate_center = None

                # Collect relevant atoms for sugar and phosphate
                for atom in residue:
                    if atom.get_name() in {'O4\'', 'C1\'', 'C2\'', 'C3\'', 'C4\''}:
                        sugar_atoms.append(atom)
                    elif atom.get_name() == 'P':
                        phosphate_center = atom.get_coord()

                # Calculate sugar centroid
                sugar_centroid = calculate_centroid(sugar_atoms) if len(sugar_atoms) == 5 else None

                # Store the results
                results.append({
                    "residue_seq": residue.get_id()[1],
                    "sugar_centroid": sugar_centroid,
                    "phosphate_center": phosphate_center
                })

    return results


# This is a utility function for writing results to an output file in a human readable format.
def write_to_human_readable_file(results, output_file):
    """Write sugar centroids and phosphate centers to an output file."""
    with open(output_file, 'w') as f:
        for center in results:
            residue_seq = center["residue_seq"]
            sugar_centroid = center["sugar_centroid"]
            phosphate_center = center["phosphate_center"]
            
            # Format the output for each residue
            sugar_str = (
                f"[{sugar_centroid[0]:.3f}, {sugar_centroid[1]:.3f}, {sugar_centroid[2]:.3f}]"
                if sugar_centroid is not None else "Not found"
            )
            phosphate_str = (
                f"[{phosphate_center[0]:.3f}, {phosphate_center[1]:.3f}, {phosphate_center[2]:.3f}]"
                if phosphate_center is not None else "Not found"
            )
            
            f.write(f"Residue {residue_seq:>4}:\n")
            f.write(f"  Sugar centroid: {sugar_str}\n")
            f.write(f"  Phosphate center: {phosphate_str}\n")
            f.write("\n")  # Add a blank line between residues


# This is a utility function for writing results to pdb format.
# It gives you 2 pdb files, one for the sugar nodes and one for the phosphate nodes.
def write_to_pdb(results, sugar_file, phosphate_file):
    # Create PDBIO object for writing
    io = PDB.PDBIO()
    sugar_structure = PDB.Structure.Structure('SugarCentroids')
    phosphate_structure = PDB.Structure.Structure('PhosphateCenters')

    # Create models
    sugar_model = PDB.Model.Model(0)
    phosphate_model = PDB.Model.Model(0)

    # Assign valid chain IDs (use A for sugar and B for phosphate)
    sugar_chain = PDB.Chain.Chain('A')
    phosphate_chain = PDB.Chain.Chain('B')

    # Unique residue counters
    sugar_residue_counter = 1
    phosphate_residue_counter = 1

    for result in results:
        # Create sugar centroid atom
        if result["sugar_centroid"] is not None:
            x, y, z = result["sugar_centroid"]
            sugar_residue = PDB.Residue.Residue((' ', sugar_residue_counter, ' '), 'CEN', '')
            sugar_atom = PDB.Atom.Atom('SUG', [x, y, z], 1.0, 0.0, ' ', 'SUG', sugar_residue_counter)
            sugar_residue.add(sugar_atom)
            sugar_chain.add(sugar_residue)
            sugar_residue_counter += 1

        # Create phosphate center atom
        if result["phosphate_center"] is not None:
            x, y, z = result["phosphate_center"]
            phosphate_residue = PDB.Residue.Residue((' ', phosphate_residue_counter, ' '), 'PHO', '')
            phosphate_atom = PDB.Atom.Atom('PHO', [x, y, z], 1.0, 0.0, ' ', 'PHO', phosphate_residue_counter)
            phosphate_residue.add(phosphate_atom)
            phosphate_chain.add(phosphate_residue)
            phosphate_residue_counter += 1

    # Add chains to models
    sugar_model.add(sugar_chain)
    phosphate_model.add(phosphate_chain)

    # Add models to structures
    sugar_structure.add(sugar_model)
    phosphate_structure.add(phosphate_model)

    # Write structures to PDB files
    io.set_structure(sugar_structure)
    io.save(sugar_file)

    io.set_structure(phosphate_structure)
    io.save(phosphate_file)


# Generates 5 output files: 1 txt file that's just for easier readability, 2 pdb files, and 2 xyz files.
# The information is the same, just different formats to avoid needing to convert btwn formats later.
def main():
    # Feel free to set the file paths here instead of at the top
    # Doing it here gives more customizability with output file paths.
    pdb_file = input_pdb_file

    output_dir = os.getcwd()     # default output directory
    if output_folder != "":
        output_dir = output_folder
    filename = os.path.basename(pdb_file)
    filename, _ = os.path.splitext(filename)    # strip off the file extension
    output_file = os.path.join(output_dir, filename + "_spbbcenters.txt")     # spbb == sugar phosphate backbone
    sugar_file = os.path.join(output_dir, filename + "_sugar_centers.pdb")
    phosphate_file = os.path.join(output_dir, filename + "_phosphate_centers.pdb")
    
    # # debugging, ignore this
    # print(pdb_file)
    # print(filename)
    # print(output_file)
    # print(sugar_file)
    # print(phosphate_file)
    # print("Does raw empty string == empty string?", r"" == "")
    # print(os.path.splitext(sugar_file)[0] + ".xyz")
    # print(os.path.splitext(phosphate_file)[0] + ".xyz")
    # assert(False)      # an assert statement that evalutes to false stops the program

    # Find centers
    centers = find_centers(pdb_file)
    
    # # Print results
    # for center in centers:
    #     print(f"Residue {center['residue_seq']:>4}:")
    #     if center['sugar_centroid'] is not None:
    #         print(f"  Sugar centroid: {center['sugar_centroid']}")
    #     else:
    #         print("  Sugar centroid: Not found")
        
    #     if center['phosphate_center'] is not None:
    #         print(f"  Phosphate center: {center['phosphate_center']}")
    #     else:
    #         print("  Phosphate center: Not found")

    # Write results to file(s)      If you only want certain files, comment out the rest. But for now,
    # you need the pdb files to get the xyz files since I'm just using a function I already had for converting pdb to xyz. 
    # But rewriting the write_to_pdb() function into a write_to_xyz() function should be pretty straightforward.
    # Conversely, extracting coordinate data out of a pdb file is pretty straightforward too. Biopython should have something.

    # All info in 1 txt file
    write_to_human_readable_file(centers, output_file)
    print(f"Results (human readable) written to: \n{output_file}")

    # Saving as 2 pdb files
    write_to_pdb(centers, sugar_file, phosphate_file)
    print(f"Results (pdb) written to: \n{sugar_file} and {phosphate_file}")

    # Saving as 2 xyz files.
    from pdb2xyz import pdb_to_xyz   # this is a local file btw
    pdb_to_xyz(sugar_file, os.path.splitext(sugar_file)[0] + ".xyz")
    pdb_to_xyz(phosphate_file, os.path.splitext(phosphate_file)[0] + ".xyz")


if __name__ == "__main__":
    main()

