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
parser.add_argument("--area", default='half')
parser.add_argument('--cmap', default='planck')
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)

# load map and ivar from the specified
# box = np.array([[-1,2],[1,-2]]) / 180*np.pi
print(args.area)
box = boxes[args.area]
imap = load_map(filedb[args.freq]['coadd'], box=box, fcode=args.freq)/1e9
ivar = load_ivar(filedb[args.freq]['coadd_ivar'], box, fcode=args.freq)*1e18

# optionally smooth the map
if args.smooth>0:
    imap = enmap.smooth_gauss(imap, args.smooth*u.fwhm*u.arcmin)
    ivar = enmap.smooth_gauss(ivar, args.smooth*u.fwhm*u.arcmin)
    s    = sfactor(args.freq, args.smooth)
else:
    s    = 1

# compute polarization fraction
P = np.sum(imap[1:]**2, axis=0)**0.5
p = P/imap[0]
# calculate p error
P_err = lib.P_error(imap, ivar*s**2)
I_err = ivar[0]**-0.5 / s
p_err = (imap[0] * P_err - P * I_err) / imap[0]**2
mask = imap[0] > 1.6
mask_max = p == np.max(p[mask])
p[~mask] = 1e-6
print(f'max: {p[mask_max]} +- {p_err[mask_max]}')

# import ipdb; ipdb.set_trace()
plot_opts = {
    'origin': 'lower',
    'cmap': args.cmap,
    'vmin': args.min,
    'vmax': args.max,
    'extent': box2extent(box)/np.pi*180
}

fig = plt.figure()
if args.log: plt.imshow(np.log10(p), **plot_opts)
else: plt.imshow(p, **plot_opts)
# p[mask_max] = 1
# plt.imshow(p, **plot_opts)
plt.xlabel('l [deg]')
plt.ylabel('b [deg]')
plt.gca().set_aspect('equal')
if args.log: plt.colorbar(shrink=0.5).set_label(r"$\log_{10}$(P/I)")
else: plt.colorbar(shrink=0.5).set_label(r"P/I")
plt.tight_layout()
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
