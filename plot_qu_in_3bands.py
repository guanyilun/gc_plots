"""This script plots Q/U maps in 3 bands"""

import argparse, os, os.path as op
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from pixell import utils as u
import plotstyle
from common import *

parser = argparse.ArgumentParser()
parser.add_argument("-o","--odir",default="plots")
parser.add_argument('--smooth', type=float, default=0)
parser.add_argument("--oname",default="QU_3bands.pdf")
parser.add_argument("--min",type=float, default=-1e8)
parser.add_argument("--max",type=float, default=1e8)
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)

# load datafile
def process_map(m, component='Q'):
    if   component == 'Q': i = 1
    elif component == 'U': i = 2
    else: raise ValueError("Only Q/U are allowed")
    if args.smooth > 0:
        omap = enmap.smooth_gauss(m[i], args.smooth*u.fwhm*u.arcmin)
    else: omap = m[i]
    return omap
# box of interests in format: [[fromy, fromx],[toy, tox]]
box = np.array([[-1,2],[1,-2]]) / 180*np.pi

# initialize plot
fig, axes = plt.subplots(3,2,figsize=(9,7))
plot_opts = {
    'origin': 'lower',
    'cmap': 'planck',
    'vmin': args.min,
    'vmax': args.max,
    'extent': box2extent(box)/np.pi*180
}
# row 0: f090
imap = load_map(filedb['f090']['coadd'], box, fcode='f090')
qmap = process_map(imap, 'Q')
im = axes[0,0].imshow(qmap, **plot_opts)
umap = process_map(imap, 'U')
axes[0,1].imshow(umap, **plot_opts)
ofile = op.join(args.odir, args.oname)
# row 1: f150
imap = load_map(filedb['f150']['coadd'], box, fcode='f150')
qmap = process_map(imap, 'Q')
axes[1,0].imshow(qmap, **plot_opts)
umap = process_map(imap, 'U')
axes[1,1].imshow(umap, **plot_opts)
ofile = op.join(args.odir, args.oname)
# row 2: f220
imap = load_map(filedb['f220']['coadd'], box, fcode='f220')
qmap = process_map(imap, 'Q')
axes[2,0].imshow(qmap, **plot_opts)
umap = process_map(imap, 'U')
axes[2,1].imshow(umap, **plot_opts)
ofile = op.join(args.odir, args.oname)

# axes settings
for i in range(3):
    for j in range(2):
        ax = axes[i,j]
        ax.set_axis_off()
        # if j>0:
        #     ax.axes.yaxis.set_ticklabels([])
        # if i<2:
        #     ax.axes.xaxis.set_ticklabels([])
        # if i==2:
        #     ax.xaxis.set_major_locator(ticker.FixedLocator([1.5,1,0.5,0,-0.5,-1,-1.5]))
        #     # ax.set_xlabel('l [deg]')
        # if j==0:
        #     ax.yaxis.set_major_locator(ticker.FixedLocator([-0.5,0,0.5]))
        #     # ax.set_ylabel('b [deg]')

# labels
# setup labels: Planck, ACT+Planck
axes[0,0].text(0.5, 1.05, 'Q',
               verticalalignment='bottom', horizontalalignment='left',
               transform=axes[0,0].transAxes, fontsize=12)
axes[0,1].text(0.5, 1.05, 'U',
               verticalalignment='bottom', horizontalalignment='left',
               transform=axes[0,1].transAxes, fontsize=12)

# setup labels: f090, f150, f220
for i, label in zip(range(3), ['f090','f150','f220']):
    axes[i,0].text(-0.1, 0.35, label, rotation=90,
                   verticalalignment='bottom', horizontalalignment='left',
                   transform=axes[i,0].transAxes, fontsize=12)

# add colorbar
# fig.subplots_adjust(right=0.8,wspace=0.02, hspace=0.02)
# cb_ax = fig.add_axes([0.83, 0.1, 0.02, 0.8])
fig.colorbar(im, ax=axes, shrink=0.8).set_label()

# plt.tight_layout()        
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
