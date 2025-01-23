TLDR: This is data generated from a real pdb file, and it corresponds to the data already in your colab notebook. 6gtg_phosphate_centers.xyz is like the "ground truth" version (not exactly true) of the data in your colab notebook. Ideally, your code should perform really well on this, since it's (hopefully) not missing nodes.




0. data_from_your_colabnotebook.xyz - These are the datapoints that were already in your colab notebook. I think they come from a cryoem image (EMD-0065) fed into an ML model, so some nodes are missing.

1. 6gtg.pdb is real data from https://www.rcsb.org/structure/6GTG. This is the fitted model for EMD-0065.
The following files are generated from this file.

2. ..._spbackbone.pdb is the same pdb but with only the lines that correspond to the sugar-phosphate backbone.
This file can help you visualize the backbone.

You can open both these files in ChimeraX to view them; might be helpful for visualizing some stuff. 
They might also be helpful for checking bond lengths/angles. I think you there's some ruler tool in ChimeraX.


But for your TSP/VRP code, you should use the following files because they simplify each phosphate and sugar into single points/nodes.
This was done by taking the center of each phosphate and the center of each sugar.
For phosphates, just take the P atom.
For sugars, take the center (centroid technically) of the O4' C1' C2' C3' C4' ring, so average the positions of those 5 atoms.

Note that these sugar nodes are therefore not real atoms in the real structure (see the .JPG). 
The phosphate nodes however are real P atoms.
3. ..._phosphate_centers.pdb and ..._sugar_centers.pdb
4. ..._phosphate_centers.xyz and ..._sugar_centers.xyz

The coordinate information is the same for both file types. Feel free to use either.
EDIT: slight bug (data is still fine to use): the atom names show up as X in the .xyz files because theyre generated from their pdb counterparts, but I used made-up atom names in the pdb files, and my pdb2xyz conversion script puts X for atom names that don't exist in real life. The coordinates are fine tho. As long as you can distinguish between sugar and phosphate nodes, these atom names won't matter.

Note that the sugar and phosphate nodes are separated into different files, so you'll need to read in two files. In your code, you can store them in separate arrays or in one array or a pandas dataframe or whatever else.