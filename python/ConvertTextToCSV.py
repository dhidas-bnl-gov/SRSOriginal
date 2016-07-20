#!/usr/bin/env python

import sys
import numpy as np


def ConvertTextToCSV (fileName) :
    f = open(fileName, 'r')

    x = []
    y = []
    z = []


    with open(fileName, 'r') as f:
        for line in f:
            v = line.split()
            if (len(v) > 0 and v[0].startswith('#')):
                continue
    
            if (len(v) < 3) :
                continue
    
            x.append(float(v[0]))
            y.append(float(v[2]))
            z.append(float(v[3]))


    x2 = list(set(x))
    y2 = list(set(y))

    nx = len(x2)
    ny = len(y2)

    x2.sort()
    y2.sort()


    print ' ,', ', '.join(map(str, x2))
    for i in xrange(ny):
        print str(y2[i]) + ',',
        for j in xrange(nx):
            if j < nx - 1:
                print str(z[i * nx + j]) + ',',
            else:
                print z[i * nx + j]





if __name__ == "__main__":
    if len(sys.argv) is not 2:
        print 'Usage: ', sys.argv[1], ' [InFile]'
        exit(0)

    ConvertTextToCSV(sys.argv[1])
