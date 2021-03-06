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
parser.add_argument("--downgrade", type=float, default=1)
parser.add_argument("--cmap", default='planck')
parser.add_argument("--min", type=float, default=0)
parser.add_argument("--max", type=float, default=20)
parser.add_argument("--smooth", type=float, default=0)
parser.add_argument("--area", default='half')
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)

# load map and ivar from the specified
box = boxes[args.area]
# this doesn't need to be done in mjy coordinate as the results
# are dimensionless
imap = load_map(filedb[args.freq]['coadd'], box=box, fcode=args.freq)
ivar = load_ivar(filedb[args.freq]['coadd_ivar'], box=box, fcode=args.freq)

# optionally smooth the map
if args.smooth > 0:
    imap  = enmap.smooth_gauss(imap, args.smooth*u.arcmin*u.fwhm)
    ivar  = enmap.smooth_gauss(ivar, args.smooth*u.arcmin*u.fwhm)
    # weight down noise level by smoothing
    ivar *= sfactor(args.freq, args.smooth)**2
# optionally downgrade the map
if args.downgrade > 1:
    # smooth with a fwhm equal to the pixel size after downgrade
    imap = imap.downgrade(args.downgrade)
    ivar = ivar.downgrade(args.downgrade)*args.downgrade**2
# compute polarization angle error
P_err = lib.Pangle_error(imap, ivar, deg=True)

plot_opts = {
    'origin': 'lower',
    'cmap': args.cmap,
    'vmin': args.min,
    'vmax': args.max,
    'extent': box2extent(box)/np.pi*180
}

fig = plt.figure()
plt.imshow(P_err, **plot_opts)
plt.xlabel('l [deg]')
plt.ylabel('b [deg]')
plt.colorbar(shrink=0.5).set_label(r"$\delta \psi$ [deg]")
plt.tight_layout()
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
