# merge_pdb
Merges single model pdb files, changing atom order if needed

Assumes that there's no dual assignment and that there's the same number of atoms.

To merge 1.pdb and 2.pdb in traj.pdb:
./merge_pdb.py 1.pdb 2.pdb > traj.pdb
