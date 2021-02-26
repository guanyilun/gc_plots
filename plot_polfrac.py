"""This script aims to produce a polarization fraction plot for a
given frequency

"""

import argparse, os, os.path as op
from matplotlib import pyplot as plt
import numpy as np
import plotstyle
from common import *
import lib
from pixell import utils as u

parser = argparse.ArgumentParser()
parser.add_argument("-o","--odir", default="plots")
parser.add_argument("--freq",default="f090")
parser.add_argument("--oname", default="polfrac.pdf")
parser.add_argument('--smooth', help="smooth in arcmin", type=float, default=0)
parser.add_argument("--min", type=float, default=-3)
parser.add_argument("--max", type=float, default=0)
parser.add_argument('--log', action='store_true')
parser.add_argument('--cmap', default='planck')
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)

# load map and ivar from the specified
box = np.array([[-1,2],[1,-2]]) / 180*np.pi
imap = load_map(filedb[args.freq]['coadd'], box, fcode=args.freq)
# ivar = load_ivar(filedb[args.freq]['coadd_ivar'], box, fcode=args.freq)

# compute polarization fraction
P = np.sum(imap[1:]**2, axis=0)**0.5
p = P/imap[0]

# optionally smooth the map
if args.smooth>0:
    p = enmap.smooth_gauss(p, args.smooth*u.fwhm*u.arcmin)
# import ipdb; ipdb.set_trace()
plot_opts = {
    'origin': 'lower',
    'cmap': args.cmap,
    'vmin': args.min,
    'vmax': args.max,
    'extent': [2,-2,-1,1]
}

fig = plt.figure()
if args.log: plt.imshow(np.log10(p), **plot_opts)
else: plt.imshow(p, **plot_opts)
plt.xlabel('l [deg]')
plt.ylabel('b [deg]')
if args.log: plt.colorbar(shrink=0.5).set_label(r"$\log_{10}$(P/I)")
else: plt.colorbar(shrink=0.5).set_label(r"P/I")
plt.tight_layout()
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
