# Script for computing the position of the point source relative to the galaxy.
# in glafic scripts we always define the galaxy to be at (0, 0)

import os, sys
import numpy as np

# this is ACS image
pixscale = 0.05

# galaxy
xg = 49.8566
yg = 50.4680
ar = 0.4733
pa = 68.8965

# qso
xq = 37.9065
yq = 51.0470

# compute the relative positions of the quasar. Set xg, yg = 0, 0

xp = (xq-xg)*pixscale
yp = (yq-yg)*pixscale

print(xp, yp, ar, pa)
