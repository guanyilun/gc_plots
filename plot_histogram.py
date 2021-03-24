"""This script supports histogram of various quantities. The primary purpose is
the histogram of the polarization angle, used to test the agreement between
histogram spread and the uncertainty calculation"""

import argparse, os, os.path as op
import plotstyle
from common import *
import lib
from matplotlib import pyplot as plt

# parser defined in common
parser.add_argument("--freq", default="f090")
parser.add_argument("--use", default="coadd")
parser.add_argument("-Q", "--quantity", help='quantity of interests', default='bangle')
parser.add_argument("--figsize", default=None)
parser.add_argument("--bins", type=int, default=50)
parser.add_argument("--ofile")
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)
if args.figsize: figsize = eval(args.figsize)
else: figsize=None
# load data
box = boxes[args.area]
imap = load_map(filedb[args.freq][args.use], box=box, fcode=args.freq)

if args.quantity == 'bangle':
    # histogram of polarization angle
    quantity = np.ravel(lib.Bangle(imap[1], imap[2], toIAU=True))/np.pi*180
    xlabel   = r"$\psi_B$ [deg]"
else: raise ValueError("Specify a quantity")
# make a figure
fig = plt.figure(figsize=None)
ax = fig.add_subplot(111)
ax.hist(quantity, bins=args.bins)
ax.set_xlabel(xlabel)
ax.set_ylabel('Count')
ofile = op.join(args.odir, args.ofile)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
