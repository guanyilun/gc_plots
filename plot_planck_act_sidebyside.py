"""This script plots planck and act temperature maps side by side as a 3x2 grid"""

import argparse, os, os.path as op
import numpy as np, glob
import matplotlib.pyplot as plt
from pixell import enmap, colorize

# setup common plotstyle
import plotstyle
# load common variables and utility functions
from common import *

def process_map(imap):
    return np.log10(imap[0])

parser = argparse.ArgumentParser()
parser.add_argument("-o","--odir", default='plots')
parser.add_argument("--oname", default="test.pdf")
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)

# Build figure
fig, axes = plt.subplots(3, 2, gridspec_kw={'hspace': 0.05, 'wspace': 0.05}, figsize=(9,7))
box = np.array([[-1,2],[1,-2]]) / 180*np.pi

plot_opts = {
    'origin': 'lower',
    'cmap': 'planck',
    'vmin': 1.5,
    'vmax': 5
}
# row 0: f090
imap = load_map(filedb['f090']['planck'], box)
imap = process_map(imap)
axes[0,0].imshow(imap, **plot_opts)
imap = load_map(filedb['f090']['coadd'], box)
imap = process_map(imap)
axes[0,1].imshow(imap, **plot_opts)

# row 1: f150
plot_opts.update({'vmax': 4.5})
imap = load_map(filedb['f150']['planck'], box)
imap = process_map(imap)
axes[1,0].imshow(imap, **plot_opts)
imap = load_map(filedb['f150']['coadd'], box)
imap = process_map(imap)
axes[1,1].imshow(imap, **plot_opts)

# row 2: f220
plot_opts.update({'vmax': 5})
imap = load_map(filedb['f220']['planck'], box)
imap = process_map(imap)
axes[2,0].imshow(imap, **plot_opts)
imap = load_map(filedb['f220']['coadd'], box)
imap = process_map(imap)
axes[2,1].imshow(imap, **plot_opts)

# turn-off axis
for i in range(3):
    for j in range(2):
        axes[i,j].axis('off')

# setup labels: Planck, ACT+Planck
axes[0,0].text(0.4, 1.05, 'Planck',
               verticalalignment='bottom', horizontalalignment='left',
               transform=axes[0,0].transAxes, fontsize=12)
axes[0,1].text(0.35, 1.05, 'ACT+Planck',
               verticalalignment='bottom', horizontalalignment='left',
               transform=axes[0,1].transAxes, fontsize=12)

# setup labels: f090, f150, f220
axes[0,0].text(-0.05, 0.35, 'f090', rotation=90,
               verticalalignment='bottom', horizontalalignment='left',
               transform=axes[0,0].transAxes, fontsize=12)
axes[1,0].text(-0.05, 0.35, 'f150', rotation=90,
               verticalalignment='bottom', horizontalalignment='left',
               transform=axes[1,0].transAxes, fontsize=12)
axes[2,0].text(-0.05, 0.35, 'f220', rotation=90,
               verticalalignment='bottom', horizontalalignment='left',
               transform=axes[2,0].transAxes, fontsize=12)

ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
