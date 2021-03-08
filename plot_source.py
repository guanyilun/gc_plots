"""A more source tailored version of plot_region"""

import argparse, os, os.path
import numpy as np

import lib
from common import *
import plotstyle
from pixell import enplot, utils as u, enmap
from matplotlib import pyplot as plt
from scipy.optimize import leastsq

parser = argparse.ArgumentParser()
parser.add_argument("-o","--odir", default="plots")
parser.add_argument("--oname", default="mouse.pdf")
parser.add_argument("--smooth", type=float,default=2)
parser.add_argument("--tonly", action='store_true', default=False)
parser.add_argument("--dust-removal", action='store_true', default=False)
parser.add_argument("--dust-factor-f090", type=float, default=1)
parser.add_argument("--dust-factor-f150", type=float, default=1)
parser.add_argument("--auto-adj", action='store_true')
parser.add_argument("-l", help="in deg", type=float, default=None)
parser.add_argument("-b", help="in deg", type=float, default=None)
parser.add_argument("--margin", help="margin in deg", type=float, default=0.1)
parser.add_argument("--cmap", default='planck')
parser.add_argument("--title", default=None)
parser.add_argument("--tmin", default="0,0,0")
parser.add_argument("--tmax", default="1,1,1")
parser.add_argument("--pmin", default="0,0,0")
parser.add_argument("--pmax", default="1,1,1")
parser.add_argument("--figsize", default="(8,7)")

args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)
tmin = args.tmin.split(",")
tmax = args.tmax.split(",")
pmin = args.pmin.split(",")
pmax = args.pmax.split(",")

# define box
xmin, xmax = args.l-args.margin, args.l+args.margin
ymin, ymax = args.b-args.margin, args.b+args.margin
box = np.array([[ymin, xmax], [ymax, xmin]]) / 180*np.pi
# set figure size if needed
if args.figsize is not None: figsize=eval(args.figsize)
else: figsize=None
# load maps
imap_f090 = load_map(filedb['f090']['coadd'], fcode='f090')/1e9
imap_f150 = load_map(filedb['f150']['coadd'], fcode='f150')/1e9
imap_f220 = load_map(filedb['f220']['coadd'], fcode='f220')/1e9

# form dust template at each frequency
beta = 1.53
s_f090 = (100/217)**(beta+2)
s_f150 = (143/217)**(beta+2)
dust_f090 = lib.beam_match(imap_f220, 'f090', 'f220')*s_f090*args.dust_factor_f090
dust_f150 = lib.beam_match(imap_f220, 'f150', 'f220')*s_f150*args.dust_factor_f150

def preprocess_T(imap, dust_removal=False, dust_template=None, box=None, auto_adj=False):
    omap = enmap.submap(imap[0], box=box)
    if dust_removal:
        dust = enmap.submap(dust_template[0], box=box)
        if auto_adj:
            par, _ = leastsq(lambda x: np.ravel(omap-x*dust), x0=1)
            factor = par[0]
        else: factor = 1
        print("factor:", factor)
        omap -= dust * factor
    return omap

def preprocess_P(imap, dust_removal=None, dust_template=None, smooth=None, box=None):
    omap = np.sum(imap[1:]**2,axis=0)**0.5
    if dust_removal:
        omap -= np.sum(dust_template[1:]**2,axis=0)**0.5
    if smooth is not None:
        omap = enmap.smooth_gauss(omap, smooth*u.fwhm*u.arcmin)
    return enmap.submap(omap, box=box)

# preprocess
I_f090 = preprocess_T(imap_f090, dust_removal=args.dust_removal,
                      dust_template=dust_f090, box=box, auto_adj=args.auto_adj)
I_f150 = preprocess_T(imap_f150, dust_removal=args.dust_removal,
                      dust_template=dust_f150, box=box, auto_adj=args.auto_adj)
I_f220 = preprocess_T(imap_f220, box=box)

P_f090 = preprocess_P(imap_f090, dust_removal=args.dust_removal,
                      dust_template=dust_f090,
                      smooth=args.smooth, box=box)
P_f150 = preprocess_P(imap_f150, dust_removal=args.dust_removal,
                      dust_template=dust_f150,
                      smooth=args.smooth, box=box)
P_f220 = preprocess_P(imap_f220, smooth=args.smooth, box=box)

popts = {
    'origin': 'lower',
    'cmap': args.cmap,
    'extent': box2extent(box)/np.pi*180
}
if args.dust_removal: names = ['f090 - f220','f150 - f220']
else: names = ['f090', 'f150', 'f220']
maps = ['T','P']
if args.tonly:
    nrow, ncol = 1, len(names)
else:
    nrow, ncol = 2, len(names)

fig, axes = plt.subplots(nrow,ncol, figsize=figsize)
data = [[I_f090, I_f150, I_f220],[P_f090, P_f150, P_f220]]
mins = [tmin, pmin]
maxs = [tmax, pmax]
for i in range(nrow):
    for j in range(ncol):
        axes[i,j].imshow(data[i][j], vmin=mins[i][j], vmax=maxs[i][j], **popts)
for i in range(nrow):
    for j in range(ncol):
        ax = axes[i,j]
        if i == 0:
            ax.text(0.3, 1.02, names[j],
                   verticalalignment='bottom', horizontalalignment='left',
                   transform=ax.transAxes, fontsize=14)
        if j == 0:
            ax.text(-0.2, 0.5, maps[i], rotation=90,
                   verticalalignment='bottom', horizontalalignment='left',
                   transform=ax.transAxes, fontsize=14)
        # axes[i,j].axis('off')
if args.title is not None: plt.suptitle(args.title, y = 1.02, fontsize=16)
plt.tight_layout(h_pad=0.1, w_pad=0.1)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
