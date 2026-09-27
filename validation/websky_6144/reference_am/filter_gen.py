#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import sys

bx   = float(sys.argv[1])
n    = int(sys.argv[2])
rmax = float(sys.argv[3])

rmincell = 1.65   # smallest smoothing radius in cellsize units
#rmincell = 2.   # smallest smoothing radius in cellsize units

spacing = 1.15 #best spacing to use for filters
#spacing = 1.2 #best spacing to use for filters
cellsize = bx / n

rmin = rmincell * cellsize

# get number of filter scales
r=rmin
i=1
rlast = rmax*0.07
while (r<rmax-rlast):    
    i += 1
    r = r*spacing
if(r<rmax-rlast): #If last filter basically at rmax don't add another at rmax
    i+=1 #add filter at rmax
print("{0:<2}".format(i)) # Prints the number of filter scales, read in as "nic" in hpkvd.f90

#make filters and save to file
r=rmin
i=1
while (r<rmax):    
    print("{0:<2}".format(i),1.686,"{:6.4e}".format(r),0.300)
    i += 1
    r = r*spacing
i = i-1
r = r/spacing

if(r<rmax-rlast):
    i=i+1
    print("{0:<2}".format(i),1.686,"{:6.4e}".format(rmax),0.300)
