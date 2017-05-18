#!/usr/bin/env python

import sys
import numpy as np

filenames = sys.argv[1:]

c = 0
for filename in filenames:
    c += 1
    atomlines = []
    with open(filename) as f:
        for line in f.readlines():
            if line.startswith("ATOM"):
                # to make sure oxt is last
                if "OXT" in line:
                    oxt = line.split()
                else:
                    atomlines.append(line.split())

    atomlines = np.asarray(atomlines)
    oxt = np.asarray(oxt)
    atomindex = np.where(atomlines == "CA")[1][0]
    numindex = atomindex - 1
    # check if there's chain identifier
    try:
        resindex = int(atomindex + 2)
        has_chain = False
    except:
        resindex = atomindex + 3
        has_chain = True

    order = np.lexsort((atomlines[:,atomindex],atomlines[:,resindex].astype(int)))
    atomlines = atomlines[order]
    atomlines[:,numindex] = (np.arange(atomlines.shape[0]) + 1).astype(str)
    oxt[numindex] = str(atomlines.shape[0] + 1)
    atomlines = np.concatenate([atomlines, oxt[None,:]])
    print "MODEL",c
    for line in atomlines:
        if has_chain:
            print "{:<6}{:>5} {:<4}{:<3}{:1}{:>4} {:>8}{:>8}{:>8}".format(*line)
        else:
            print "{:<6}{:>5} {:<4}{:<3} {:>4} {:>8}{:>8}{:>8}".format(*line)
    print "TER"
    print "ENDMDL"

