############################################################################
### This is a utility function for converting a pdb file to an xyz file. ###
### You can run this file as a script or import it as a module.          ###
############################################################################

def pdb_to_xyz(pdb_file, xyz_file):
    """
    Convert a PDB file to an XYZ file.

    Args:
        pdb_file (str): Path to the input PDB file.
        xyz_file (str): Path to the output XYZ file.
    """
    atoms = []
    with open(pdb_file, 'r') as pdb:
        for line in pdb:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Extract atom name and coordinates
                atom = line[76:78].strip()  # Atom type (e.g., C, O, N)
                x = float(line[30:38].strip())  # X-coordinate
                y = float(line[38:46].strip())  # Y-coordinate
                z = float(line[46:54].strip())  # Z-coordinate
                atoms.append((atom, x, y, z))

    with open(xyz_file, 'w') as xyz:
        xyz.write(f"{len(atoms)}\n")
        xyz.write("Converted from PDB\n")
        for atom, x, y, z in atoms:
            xyz.write(f"{atom} {x:.3f} {y:.3f} {z:.3f}\n")

# Example usage:
if __name__ == "__main__":
    # Put your files in the
    pdb_file = r""  # Path to your PDB file
    xyz_file = r""  # Path to save the XYZ file
    pdb_to_xyz(pdb_file, xyz_file)