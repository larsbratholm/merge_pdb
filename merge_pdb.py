#!/usr/bin/env python

import sys
import numpy as np

filenames = sys.argv[1:]

ordering = {
            "ALA": dict([(j,i) for i,j in list(enumerate(["N","CA","C","O","CB"]))]),
            "ARG": dict([(j,i) for i,j in list(enumerate(["N","CA","C","O","CB","CG","CD","NE","CZ","NH1","NH2"]))]),
            "ASN": dict([(j,i) for i,j in list(enumerate(["N","CA","C","O","CB","CG","OD1","ND2"]))]),
            "ASP": dict([(j,i) for i,j in list(enumerate(["N","CA","C","O","CB","CG","OD1","OD2"]))]),
            "CYS": dict([(j,i) for i,j in list(enumerate(["N","CA","C","O","CB","SG"]))]),
            "GLN": dict([(j,i) for i,j in list(enumerate(["N","CA","C","O","CB","CG","CD","OE1","NE2"]))]),
            "GLU": dict([(j,i) for i,j in list(enumerate(["N","CA","C","O","CB","CG","CD","OE1","OE2"]))]),
            "GLY": dict([(j,i) for i,j in list(enumerate(["N","CA","C","O"]))]),
            "HIS": dict([(j,i) for i,j in list(enumerate(["N","CA","C","O","CB","CG","ND1","CD2","CE1","NE2"]))]),
            "ILE": dict([(j,i) for i,j in list(enumerate(["N","CA","C","O","CB","CG1","CG2","CD1"]))]),
            "LEU": dict([(j,i) for i,j in list(enumerate(["N","CA","C","O","CB","CG","CD1","CD2"]))]),
            "LYS": dict([(j,i) for i,j in list(enumerate(["N","CA","C","O","CB","CG","CD","CE","NZ"]))]),
            "MET": dict([(j,i) for i,j in list(enumerate(["N","CA","C","O","CB","CG","SD","CE"]))]),
            "MSE": dict([(j,i) for i,j in list(enumerate(["N","CA","C","O","CB","CG","SE","CE"]))]),
            "PHE": dict([(j,i) for i,j in list(enumerate(["N","CA","C","O","CB","CG","CD1","CD2","CE1","CE2","CZ"]))]),
            "PRO": dict([(j,i) for i,j in list(enumerate(["N","CA","C","O","CB","CG","CD"]))]),
            "SER": dict([(j,i) for i,j in list(enumerate(["N","CA","C","O","CB","OG"]))]),
            "THR": dict([(j,i) for i,j in list(enumerate(["N","CA","C","O","CB","OG1","CG2"]))]),
            "TRP": dict([(j,i) for i,j in list(enumerate(["N","CA","C","O","CB","CG","CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2"]))]),
            "TYR": dict([(j,i) for i,j in list(enumerate(["N","CA","C","O","CB","CG","CD1","CD2","CE1","CE2","CZ","OH"]))]),
            "VAL": dict([(j,i) for i,j in list(enumerate(["N","CA","C","O","CB","CG1","CG2"]))])
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
                line.replace("OC2","O  ")

                tokens = [line[:6].strip(), line[6:11].strip(), line[12:16].strip(), line[17:20].strip(), line[21:22].strip(), line[22:26].strip(), line[30:38].strip(), line[38:46].strip(), line[46:54].strip()]
                if len(tokens[2]) != 4:
                    tokens[2] = " " + tokens[2]

                if "OXT" in line:
                    oxt = tokens[:]
                else:
                    atomlines.append(tokens[:])

    atomlines = np.asarray(atomlines)
    if oxt != None:
        oxt = np.asarray(oxt)

    last_resid = -99999999
    id_counter = 0
    for i, line in enumerate(atomlines):
        resid = line[5]
        aa = line[3]
        atom = line[2]
        if i == atomlines.shape[0]-1:
            if atom in ordering[aa]:
                order.append(ordering[aa][atom])
            else: # hydrogens
                order.append(99)
            size = len(ordering[aa])
            sorted_order = np.argsort(order)
            atomlines[i-id_counter:] = atomlines[i-id_counter:][sorted_order]
        elif last_resid == resid:
            id_counter += 1
            if atom in ordering[aa]:
                order.append(ordering[aa][atom])
            else: # hydrogens
                order.append(99)
        else:
            size = len(ordering[aa])
            last_resid = resid
            if i != 0:
                sorted_order = np.argsort(order)
                atomlines[i-id_counter:i] = atomlines[i-id_counter:i][sorted_order]
            order = []
            id_counter = 0


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
        print "{:<6}{:>5}{:<4}{:<3}{:1}{:>4}{:>8}{:>8}{:>8}".format(*line)
    print "TER"
    print "ENDMDL"

