"""This will make a two panel plot with upper panel being
the total intensity and lower panel being the polarization
intensity for comparison. It actually does no computing but
to load saved maps computed separately using plot_multifreq
or plot_multifreq2 (for polarization). 

"""

import argparse, os, os.path as op
import numpy as np
import matplotlib.pyplot as plt
from pixell import enmap, enplot, utils as u
from matplotlib import colors

from common import *
import lib
import plotstyle

# parser defined in common
parser.add_argument("-T", default=None)
parser.add_argument("-P", default=None)
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)

# load two maps
tmap = np.load(args.T)
pmap = np.load(args.P)

# start plotting
popts = {
    'origin': 'lower',
}

# plots:
# -> upper panel: temperature
fig, axes = plt.subplots(2, 1, figsize=(8,4))
axes[0].imshow(tmap, **popts)
axes[1].imshow(pmap, **popts)
for ax in axes.flat: ax.axis('off')
plt.tight_layout(h_pad=0.5)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
