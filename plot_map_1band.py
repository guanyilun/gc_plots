"""This script plots only one coadd map from a given frequency.
This is to prepare for the object identification work

"""

import argparse, os, os.path as op
import numpy as np
import matplotlib.pyplot as plt
from pixell import enmap, utils as u

# setup common plotstyle
import plotstyle
# load common variables and utility functions
from common import *

parser = argparse.ArgumentParser()
parser.add_argument("-o","--odir", default='plots')
parser.add_argument("--oname", default="coadd_map.pdf")
parser.add_argument("--freq", default="f150")
parser.add_argument("--min", type=float, default=1.5)
parser.add_argument("--max", type=float, default=5)
parser.add_argument("--deconv", action="store_true", default=False)
parser.add_argument("--comp", type=int, default=0)
parser.add_argument("--smooth", type=float, default=None)
parser.add_argument("--cmap", default="planck")
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)

fwhms = {
    'f090': 2.05,
    'f150': 1.40,
    'f220': 0.98,
}

def process_map(imap, box=None, deconvolve=False, freq=None, smooth=None):
    if deconvolve:
        fwhm = fwhms[freq]
        l = imap.modlmap()
        bmap = np.exp(-0.5*l**2*(fwhm*u.fwhm*u.arcmin)**2)
        bmap = np.maximum(bmap, 1e-3)
        fmap = enmap.fft(imap)
        omap = enmap.ifft(fmap/bmap).real
    else:
        omap = imap
    if smooth is not None:
        omap = enmap.smooth_gauss(omap, smooth*u.fwhm*u.arcmin)
    if box is not None: omap = omap.submap(box)
    if args.comp == 0: return np.log10(omap[0])
    if args.comp == 12: return np.sum(omap[1:]**2, axis=0)**0.5
    else: return omap[args.comp]

# initialize figure
fig = plt.figure(figsize=(12,6))
# box = np.array([[-1,2],[1,-2]]) / 180*np.pi
box = None

plot_opts = {
    'origin': 'lower',
    'cmap': args.cmap,
    'vmin': args.min,
    'vmax': args.max,
    # 'extent': [2,-2,-1,1],
}

imap = process_map(load_map(filedb[args.freq]['coadd']),
                   box=box, freq=args.freq, deconvolve=args.deconv,
                   smooth=args.smooth)
plt.imshow(imap, **plot_opts)
# plt.grid(c='k', alpha=0.5)
plt.axis('off')
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
