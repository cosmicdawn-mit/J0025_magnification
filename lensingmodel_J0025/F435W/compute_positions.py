import os, sys
import numpy as np

# this is ACS image
pixscale = 0.05

# galaxy
xg = 50.8655
yg = 50.0253
ar = 0.4751
pa = 69.7385

# qso
xq = 37.9065
yq = 51.0470

# compute the relative positions of the quasar. Set xg, yg = 0, 0

xp = (xq-xg)*pixscale
yp = (yq-yg)*pixscale

print(xp, yp, ar, pa)
