"""This script aims to produce a polarization fraction plot for a
given frequency

"""

import argparse, os, os.path as op
from matplotlib import pyplot as plt
import numpy as np
import plotstyle
from common import *
import lib

parser = argparse.ArgumentParser()
parser.add_argument("-o","--odir", default="plots")
parser.add_argument("--freq",default="f090")
parser.add_argument("--oname", default="polfrac.pdf")
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)

# load map and ivar from the specified
box = np.array([[-1,2],[1,-2]]) / 180*np.pi
imap = load_map(filedb[args.freq]['coadd'], box)
ivar = load_ivar(filedb[args.freq]['coadd_ivar'], box)

# compute polarization fraction
P = np.sum(imap[1:]**2, axis=0)**0.5
p = P/imap[0]
# import ipdb; ipdb.set_trace()
plot_opts = {
    'origin': 'lower',
    'cmap': 'planck',
    'vmin': -3,
    'vmax': 0,
    'extent': [2,-2,-1,1]
}

fig = plt.figure()
plt.imshow(np.log10(p), **plot_opts)
plt.xlabel('l [deg]')
plt.ylabel('b [deg]')
plt.colorbar(shrink=0.5).set_label(r"$\log_{10}$(P/I)")
plt.tight_layout()
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
