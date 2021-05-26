#!/usr/bin/env python2

import sys
import os

hits = {}
for line in sys.stdin:
    cols = line.strip("\n").split("\t")
    ref = cols[0]
    mag = cols[1]
    dist = float(cols[2])
    pvalue = float(cols[3])
    if os.path.basename(ref) != os.path.basename(mag):
        if mag not in hits.keys():
            hits[mag] = [ref,dist,pvalue]
        else:
            if hits[mag][1] > dist:
                hits[mag] = [ref,dist,pvalue]
            elif hits[mag][1] == dist:
                if hits[mag][2] > pvalue:
                    hits[mag] = [ref,dist,pvalue]

for ele in hits.keys():
    print "%s\t%s\t%.4f\t%.4f" % (ele, hits[ele][0], hits[ele][1], hits[ele][2])
    
