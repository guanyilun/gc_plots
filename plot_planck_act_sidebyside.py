"""This script plots planck and act temperature maps side by side as a 3x2 grid"""

import argparse, os, os.path as op
import numpy as np, glob
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors
from pixell import enmap, colorize

# setup common plotstyle
import plotstyle
# load common variables and utility functions
from common import *

def process_map(imap):
    return np.log10(imap[0])

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    """build a colormap by truncating a known colormap"""
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

parser = argparse.ArgumentParser()
parser.add_argument("-o","--odir", default='plots')
parser.add_argument("--oname", default="test.pdf")
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)

# Build figure
fig, axes = plt.subplots(3, 2, figsize=(9,7))
box = np.array([[-1,2],[1,-2]]) / 180*np.pi

planck_half = truncate_colormap(plt.get_cmap('planck'), minval=0.5)

plot_opts = {
    'origin': 'lower',
    'cmap': planck_half,
    'vmin': 8.5,
    'vmax': 10.5
}
# row 0: f090
imap = load_map(filedb['f090']['planck'], box, fcode='f090')
imap = process_map(imap)
axes[0,0].imshow(imap, **plot_opts)
imap = load_map(filedb['f090']['coadd'], box, fcode='f090')
imap = process_map(imap)
axes[0,1].imshow(imap, **plot_opts)

# row 1: f150
imap = load_map(filedb['f150']['planck'], box, fcode='f150')
imap = process_map(imap)
axes[1,0].imshow(imap, **plot_opts)
imap = load_map(filedb['f150']['coadd'], box, fcode='f150')
imap = process_map(imap)
axes[1,1].imshow(imap, **plot_opts)

# row 2: f220
plot_opts.update({'vmin':9.5, 'vmax': 11})
imap = load_map(filedb['f220']['planck'], box, fcode='f220')
imap = process_map(imap)
axes[2,0].imshow(imap, **plot_opts)
imap = load_map(filedb['f220']['coadd'], box, fcode='f220')
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
axes[0,0].text(-0.05, 0.40, 'f090', rotation=90,
               verticalalignment='bottom', horizontalalignment='left',
               transform=axes[0,0].transAxes, fontsize=12)
axes[1,0].text(-0.05, 0.40, 'f150', rotation=90,
               verticalalignment='bottom', horizontalalignment='left',
               transform=axes[1,0].transAxes, fontsize=12)
axes[2,0].text(-0.05, 0.40, 'f220', rotation=90,
               verticalalignment='bottom', horizontalalignment='left',
               transform=axes[2,0].transAxes, fontsize=12)
plt.tight_layout(w_pad=0.2, h_pad=0.2)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
