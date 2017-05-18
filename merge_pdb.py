#!/usr/bin/env python

import sys
import numpy as np

resindex = 4
numindex = 1
atomindex = 2

filenames = sys.argv[1:]

c = 0
for filename in filenames:
    c += 1
    atomlines = []
    with open(filename) as f:
        for line in f.readlines():
            if line.startswith("ATOM"):
                atomlines.append(line.split())

    atomlines = np.asarray(atomlines)
    order = np.lexsort((atomlines[:,atomindex],atomlines[:,resindex]))
    atomlines = atomlines[order]
    atomlines[:,numindex] = (np.arange(atomlines.shape[0]) + 1).astype(str)
    print "MODEL",c
    for line in atomlines:
        print "{:<6s}{:>5s} {:<4s}{:<3} {:>4} {:>8}{:>8}{:>8}".format(*line)
    print "TER"
    print "ENDMDL"

