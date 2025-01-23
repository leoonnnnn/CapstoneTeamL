#################################################################################
### This script reads the PDB file and extracts the sugar-phosphate backbone. ###
#################################################################################

# PDB file format (go to page 2): https://www.biostat.jhsph.edu/~iruczins/teaching/260.655/links/pdbformat.pdf


file_path = r"C:\Users\noelu\Capstone 2024-25 (177B as a helper bruh)\code sandbox\bio_structures\real\6eri-pdb-bundle1.pdb"
output_path = r"C:\Users\noelu\Capstone 2024-25 (177B as a helper bruh)\code sandbox\bio_structures\artificial\6eri-pdb-bundle1_spbackbone.pdb"

# Define backbone atom names for RNA
backbone_atoms = {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'", "O2'"}   # O2' can be removed (hydroxyl group on 2' carbon)

# Open the input file and filter lines for the backbone atoms
with open(file_path, "r") as pdb_file, open(output_path, "w") as output_file:
    for line in pdb_file:
        # if line.startswith("ATOM") or line.startswith("HETATM"):
        if line.startswith("ATOM"):
            atom_name = line[13:16].strip()    # change to 12:16 if including HETATM
            # print(atom_name)
            if atom_name in backbone_atoms:
                output_file.write(line)
        elif line.startswith("TER") or line.startswith("END"):  # Keep termination lines
            output_file.write(line)

print(f"output file location:\n{output_path}")
