#!/usr/bin/env python

import sys
import numpy as np

filenames = sys.argv[1:]

ordering = {
            "ALA": dict(list(enumerate(["N","CA","C","O","CB"]))),
            "ARG": dict(list(enumerate(["N","CA","C","O","CB","CG","CD","NE","CZ","NH1","NH2"]))),
            "ASN": dict(list(enumerate(["N","CA","C","O","CB","CG","OD1","ND2"]))),
            "ASP": dict(list(enumerate(["N","CA","C","O","CB","CG","OD1","OD2"]))),
            "CYS": dict(list(enumerate(["N","CA","C","O","CB","SG"]))),
            "GLN": dict(list(enumerate(["N","CA","C","O","CB","CG","CD","OE1","NE2"]))),
            "GLU": dict(list(enumerate(["N","CA","C","O","CB","CG","CD","OE1","OE2"]))),
            "GLY": dict(list(enumerate(["N","CA","C","O"]))),
            "HIS": dict(list(enumerate(["N","CA","C","O","CB","CG","ND1","CD2","CE1","NE2"]))),
            "ILE": dict(list(enumerate(["N","CA","C","O","CB","CG1","CG2","CD1"]))),
            "LEU": dict(list(enumerate(["N","CA","C","O","CB","CG","CD1","CD2"]))),
            "LYS": dict(list(enumerate(["N","CA","C","O","CB","CG","CD","CE","NZ"]))),
            "MET": dict(list(enumerate(["N","CA","C","O","CB","CG","SD","CE"]))),
            "MSE": dict(list(enumerate(["N","CA","C","O","CB","CG","SE","CE"]))),
            "PHE": dict(list(enumerate(["N","CA","C","O","CB","CG","CD1","CD2","CE1","CE2","CZ"]))),
            "PRO": dict(list(enumerate(["N","CA","C","O","CB","CG","CD"]))),
            "SER": dict(list(enumerate(["N","CA","C","O","CB","OG"]))),
            "THR": dict(list(enumerate(["N","CA","C","O","CB","OG1","CG2"]))),
            "TRP": dict(list(enumerate(["N","CA","C","O","CB","CG","CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2"]))),
            "TYR": dict(list(enumerate(["N","CA","C","O","CB","CG","CD1","CD2","CE1","CE2","CZ","OH"]))),
            "VAL": dict(list(enumerate(["N","CA","C","O","CB","CG1","CG2"])))
            }

c = 0
for filename in filenames:
    c += 1
    atomlines = []
    last_names = None
    oxt = None
    with open(filename) as f:
        for line in f.readlines():
            line = line[:56]
            if line.startswith("ATOM"):
                # to make sure oxt is last
                line.replace("OC1","OXT")
                line.replace("OC2","O")

                tokens = [line[:6], line[6:11], line[12:16], line[17:20], line[21:22], line[22:26], line[30:38], line[38:46], line[46:54]]

                if "OXT" in line:
                    oxt = tokens[:]
                else:
                    atomlines.append(tokens[:])

    atomlines = np.asarray(atomlines)
    if oxt != None:
        oxt = np.asarray(oxt)

    last_resid = atomlines[0,5]
    id_counter = 0
    for i, line in enumerate(atomlines):
        resid = line[5]
        aa = line[4]
        atom = line[2]
        id_counter += 1
        if last_resid != resid:
            size = len(ordering[aa])
            last_resid = resid
            sorted_order = np.argsort(order)
            atomlines[i-id_counter:i] = atomlines[i-id_counter_i][sorted_order]
            order = []
            id_counter = 0

        if atom in ordering[aa]:
            order.append(ordering[aa][atom])
        else: # hydrogens
            order.append(99)


    atomlines = atomlines[order]
    atomlines[:,1] = (np.arange(atomlines.shape[0]) + 1).astype(str)
    if oxt != None:
        oxt[1] = str(atomlines.shape[0] + 1)
        atomlines = np.concatenate([atomlines, oxt[None,:]])
    if last_names != None:
        assert(atomlines[:,1] == last_names)
    else:
        last_names = atomlines[:,1]
    print "MODEL",c
    for line in atomlines:
        print "{:<6}{:>5} {:<4}{:<3}{:1}{:>4} {:>8}{:>8}{:>8}".format(*line)
    print "TER"
    print "ENDMDL"

