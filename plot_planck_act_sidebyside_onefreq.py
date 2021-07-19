"""This script plots planck and act temperature maps side by side as a 3x2 grid

2021 Mar 26
- add colorbar and axis coordinates
"""

import argparse, os, os.path as op
import numpy as np, glob
import matplotlib as mpl
mpl.rcParams['xtick.direction']='in'
mpl.rcParams['ytick.direction']='in'
mpl.rcParams["xtick.major.size"]=2
mpl.rcParams["ytick.major.size"]=2
mpl.rcParams["xtick.minor.size"]=1
mpl.rcParams["ytick.minor.size"]=1

import matplotlib.pyplot as plt
from matplotlib import colors
from pixell import enmap, colorize

# setup common plotstyle
import plotstyle
# load common variables and utility functions
from common import *

parser = argparse.ArgumentParser()
parser.add_argument("-o","--odir", default='plots')
parser.add_argument("--oname", default="test.pdf")
parser.add_argument("--use", default='coadd')
parser.add_argument("--area", default='half')
parser.add_argument("--figsize", default="(9,6.5)")
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)
if args.figsize: figsize = eval(args.figsize)
else: figsize = None
box = boxes[args.area]

# Build figure
# load a temperorary file to get wcs right
imap = load_map(filedb['f090']['planck'], box, fcode='f090')/1e9
fig, axes = plt.subplots(1, 2, figsize=figsize, gridspec_kw={'width_ratios':[0.93,1]},
                         subplot_kw={'projection': imap.wcs})

plot_opts = {
    'origin': 'lower',
    'cmap': 'planck_half',
}

for freq in ['f090']:
    for j in range(2):
        # plotstyle
        if   freq == 'f090': norm = colors.LogNorm(vmin=10**-0.5, vmax=10**1.5)
        elif freq == 'f150': norm = colors.LogNorm(vmin=10**-0.5, vmax=10**1.5)
        else:        norm = colors.LogNorm(vmin=10**0.5,  vmax=10**2)
        # setup axes projection
        ax = axes[j]
        if j == 1:
            xticks, yticks = True, False
        else:
            xticks, yticks = True, True
        if j == 0:
            ax.set_ylabel('$b$')
        ax.set_xlabel('$l$')
        plotstyle.setup_axis(ax, fmt="d.d", xticks=xticks, yticks=yticks, nticks=[10,5])
        # load data
        if j == 0: use = 'act'
        if j == 1: use = 'coadd'
        imap = load_map(filedb[freq][use], box, fcode=freq)/1e9
        imap[0,imap[0]<10**-0.5] = 10**-0.5
        im = ax.imshow(imap[0], norm=norm, **plot_opts)
        cb = plotstyle.add_colorbar(fig, axes[1])
        fig.colorbar(im, cax=cb).set_label(texify('Total Intensity [MJy/sr]'), fontsize=14)


# setup labels: Planck, ACT+Planck
axes[0].text(0.4, 1.05, texify('ACT'),
               verticalalignment='bottom', horizontalalignment='left',
               transform=axes[0].transAxes, fontsize=14)
axes[1].text(0.35, 1.05, texify('ACT+Planck'),
               verticalalignment='bottom', horizontalalignment='left',
               transform=axes[1].transAxes, fontsize=14)

plt.tight_layout(w_pad=0.01, h_pad=0)
props = dict(alpha=1, facecolor='white')
plt.text(0.44, 0.75,  texify('f090'), transform=fig.transFigure, fontsize=12,
         usetex=True, horizontalalignment='left', bbox=props)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
