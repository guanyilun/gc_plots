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
parser.add_argument("--axis", action='store_true')
parser.add_argument("--figsize", default="(8,4)")
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)
if args.figsize: figsize=eval(args.figsize)
else: figsize=None
# define box of interests
box = boxes[args.area]
# load a map for wcs only
imap = load_map(filedb['f090']['coadd'], box=box, fcode='f090')
# load two maps
tmap = np.load(args.T)
pmap = np.load(args.P)

# start plotting
popts = {
    'origin': 'lower',
    # 'extent': box2extent(box)/np.pi * 180
}

# plots:
# -> upper panel: temperature
fig, axes = plt.subplots(2, 1, figsize=figsize, subplot_kw={'projection':imap.wcs})
axes[0].imshow(tmap, **popts)
axes[1].imshow(pmap, **popts)
if not args.axis:
    for ax in axes.flat: ax.axis('off')
    plt.tight_layout(h_pad=0.5)
else:
    for ax in axes:
        ax.tick_params(axis='x', colors='white', which='both', labelcolor='black')
        ax.tick_params(axis='y', colors='white', which='both', labelcolor='black')
        ax.set_aspect('equal')
        for side in ['left','right','top','bottom']:
            ax.spines[side].set_visible(True)    
            ax.spines[side].set_color('white')    
    plotstyle.setup_axis(axes[0])
    plotstyle.setup_axis(axes[1])
    axes[0].set_ylabel(" ")
    axes[1].set_ylabel(" ")
    axes[0].set_xticklabels([])
    axes[0].set_yticklabels([])
    axes[1].set_xlabel('l')
    axes[1].set_ylabel('b')
    plt.tight_layout(h_pad=0.1)

ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
