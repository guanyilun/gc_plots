"""This script plots planck and act polarization map side by side as a 1x2 grid
for one freq only"""

import argparse, os, os.path as op
import numpy as np, glob
import matplotlib.pyplot as plt
from pixell import enmap, colorize

# setup common plotstyle
import plotstyle
# load common variables and utility functions
from common import *

def process_map(imap):
    return np.sum(imap[1:]**2, axis=0)**0.5 / 1e9

# parser defined in common
parser.add_argument("--freq", default="f090")
parser.add_argument("--figsize", default=None)
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)
if args.figsize: figsize = eval(args.figsize)
else: figsize = None

# Build figure
fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)
box = boxes[args.area]

plot_opts = {
    'origin': 'lower',
    'cmap': args.cmap,
    'vmin': args.min,
    'vmax': args.max
}
# plot planck
imap = load_map(filedb[args.freq]['planck'], box, fcode=args.freq)
imap = process_map(imap)
axes[0].imshow(imap, **plot_opts)
imap = load_map(filedb[args.freq]['coadd'], box, fcode=args.freq)
imap = process_map(imap)
im = axes[1].imshow(imap, **plot_opts)
fig.subplots_adjust(right=0.9, wspace=0)
cax = fig.add_axes([0.91, 0.12, 0.02, 0.76])
fig.colorbar(im, cax=cax).set_label('Polarization intensity [MJy/sr]')
# turn-off axis
for ax in axes.flat: ax.axis('off')

# setup labels: Planck, ACT+Planck
axes[0].text(0.4, 1.05, 'Planck',
             verticalalignment='bottom', horizontalalignment='left',
             transform=axes[0].transAxes, fontsize=12)
axes[1].text(0.35, 1.05, 'ACT+Planck',
             verticalalignment='bottom', horizontalalignment='left',
             transform=axes[1].transAxes, fontsize=12)

# setup labels: f090, f150, f220
axes[0].text(-0.05, 0.40, args.freq, rotation=90,
             verticalalignment='bottom', horizontalalignment='left',
             transform=axes[0].transAxes, fontsize=12)
# plt.tight_layout(w_pad=0.2)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
