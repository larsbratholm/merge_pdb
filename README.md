# merge_pdb
Merges single model pdb files, changing atom order if needed

Assumes that the chain column is blank, that there's no dual assignment and that there's the same number of atoms.
Atom names are sorted, which might break certain readers.

To merge 1.pdb and 2.pdb in traj.pdb:
./merge_pdb.py 1.pdb 2.pdb > traj.pdb
