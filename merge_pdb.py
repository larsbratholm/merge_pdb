#!/usr/bin/env python

import sys
import numpy as np

filenames = sys.argv[1:]

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

    order = np.lexsort((atomlines[:,2],atomlines[:,5].astype(int)))

    atomlines = atomlines[order]
    atomlines[:,numindex] = (np.arange(atomlines.shape[0]) + 1).astype(str)
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

