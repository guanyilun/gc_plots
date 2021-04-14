import argparse, os, os.path as op
from common import *
import lib
from matplotlib import pyplot as plt
from pixell import enmap
import plotstyle

# parser defined in common
parser.add_argument("--freq", default='f090')
parser.add_argument("--smooth", type=float, default=None)
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)

box = boxes[args.area]
imap = load_map(filedb[args.freq]['coadd'], box=box, fcode=args.freq) / 1e9
if args.smooth:
    imap = enmap.smooth_gauss(imap, args.smooth*u.fwhm*u.arcmin)  # u defined in common
# calculate polarization angle
Pangle = 90-lib.Bangle(imap[1], imap[2], toIAU=True) / np.pi * 180
print(f"Pangle = {np.median(Pangle)} +- {np.std(Pangle, ddof=1)}")
plt.hist(np.ravel(Pangle), bins=100)
plt.xlabel('Tilt w.r.t Galactic plane [deg]')
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
