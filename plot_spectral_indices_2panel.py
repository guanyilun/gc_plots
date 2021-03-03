"""This script aims to plot a map of spectral indices
in terms of its difference to dust and synchrotron spectral
indices. Two maps are plotted as a two-panel image with
upper panel being f090/f150, and lower pannel being f150/f220

"""

import argparse, os, os.path as op
import numpy as np

import lib
from common import *
import plotstyle
from matplotlib import pyplot as plt

def calc_beta(map1, map2, f1, f2):
    return np.log(map1[0]/map2[0])/np.log(f1/f2) - 2

def beam_match(imap, f1, f2):
    """f1, f2 are fcodes instead of the actual frequency centers. It
    assumes that the first one (f1) has larger beam, so f2 will be matched
    to it.
    """
    l = imap.modlmap()
    bmap_f1 = np.exp(-0.5*l**2*(fwhms[f1]*u.fwhm*u.arcmin)**2)
    bmap_f2 = np.exp(-0.5*l**2*(fwhms[f2]*u.fwhm*u.arcmin)**2)
    rmap = enmap.ifft(enmap.fft(imap) * (bmap_f1 / np.maximum(bmap_f2, 1e-3))).real
    return rmap

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--odir', default='plots')
parser.add_argument('-v', action='store_true', default=False)
parser.add_argument("--oname", default="spec.pdf")
parser.add_argument("--area", default='half')
parser.add_argument("--cmap", default='planck')
parser.add_argument("--figsize",default=None)
parser.add_argument("--f1", default='f090')
parser.add_argument("--f2", default='f150')
parser.add_argument("--f3", default='f220')
parser.add_argument("--min1", type=float, default=None)
parser.add_argument("--max1", type=float, default=None)
parser.add_argument("--min2", type=float, default=None)
parser.add_argument("--max2", type=float, default=None)
parser.add_argument("--beam-match", action='store_true')
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)
if args.figsize: figsize = eval(args.figsize)
else: figsize = None

# first load maps from different frequencies
# define a box of interests
box = boxes[args.area]

# make two panel figure
fig, axes = plt.subplots(2,1, figsize=figsize, sharex=True)
opts = {
    'origin': 'lower',
    'cmap': args.cmap, 
    'extent': box2extent(box)/np.pi * 180,
    'vmin': args.min1,
    'vmax': args.max1
}
label = r"$\beta$"

# upper pannel:
# calculate the ratio between f1 and f2
map_1 = load_map(filedb[args.f1]['coadd'], fcode=args.f1, box=box)
map_2 = load_map(filedb[args.f2]['coadd'], fcode=args.f2, box=box)
if args.beam_match: map_2 = beam_match(map_2, args.f1, args.f2)
f1    = fcenters[args.f1]
f2    = fcenters[args.f2]
bmap  = calc_beta(map_1, map_2, f1, f2)
# plotting
im1   = axes[0].imshow(bmap, **opts)
axes[0].set_ylabel('b [deg]')
# axes[0].set_title(f'{f1}-{f2}')
fig.colorbar(im1, ax=axes[0], orientation='vertical', shrink=0.9).set_label(label)

opts = {
    'origin': 'lower',
    'cmap': args.cmap, 
    'extent': box2extent(box)/np.pi * 180,
    'vmin': args.min2,
    'vmax': args.max2
}
# lower pannel f2/f3
# map_2 loaded previously
map_3 = load_map(filedb[args.f3]['coadd'], fcode=args.f3, box=box)
if args.beam_match: map_3 = beam_match(map_3, args.f2, args.f3)
f3    = fcenters[args.f3]
bmap  = calc_beta(map_2, map_3, f2, f3)
# plotting
im2   = axes[1].imshow(bmap, **opts)
fig.colorbar(im2, ax=axes[1], orientation='vertical', shrink=0.9).set_label(label)
axes[1].set_ylabel('b [deg]')
axes[1].set_xlabel('l [deg]')
# axes[1].set_title(f'{f1}-{f2}')

plt.tight_layout(h_pad=0.2, w_pad=0.0)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
