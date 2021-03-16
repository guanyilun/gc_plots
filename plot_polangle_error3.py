"""This script aims to produce a polarization fraction plot for a
given frequency

update: this script tests the use of full covariance matrix from ACT
only maps result: no noticable difference. This can also be used to
estimate the error from ignoring covariance between Q and U, so the
script is worth keeping.

"""

import argparse, os, os.path as op
from matplotlib import pyplot as plt
import numpy as np
import plotstyle
from common import *
import lib
from pixell import utils as u
from mpl_toolkits.axes_grid1 import make_axes_locatable

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
parser.add_argument("--use", default='coadd')
parser.add_argument("--method", type=int, default=1)
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)

# not pretty
divs = {
    'f090': enmap.read_map("/home/yilung/work/mapdata/gc_8x4_201018/gc_f090_night/sky_div.fits"),
    'f150': enmap.read_map("/home/yilung/work/mapdata/gc_8x4_201018/gc_f150_night/sky_div.fits"),
    'f220': enmap.read_map("/home/yilung/work/mapdata/gc_8x4_201018/gc_f220_night/sky_div.fits")
}

def Pangle_error(imap, div, deg=True):
    """modofied from the lib version, uses covariance matrix instead"""
    Q, U = imap[1], imap[2]
    P = np.sqrt(Q**2+U**2)
    cov = enmap.multi_pow(div, -1)  # might be slow
    QQ = cov[1,1]
    UU = cov[2,2]
    QU = cov[1,2]
    # QU = 0
    # import ipdb; ipdb.set_trace()
    if deg: factor = 28.65  # 0.5 (in radian) converted to deg
    else:   factor = 0.5    # in radian
    err = factor*np.sqrt(U**2*QQ + Q**2*UU - 2*Q*U*QU) / (P**2)
    return err    
    
# define area of interests
box = boxes[args.area]
# define levels to show
levels = [0,5,10,15,30,45]
fig, axes = plt.subplots(3, 1, sharex=True, figsize=(7,7))
for fi, freq in enumerate(['f090','f150','f220']):
    # load map and ivar from the specified
    imap = load_map(filedb[freq]['coadd'], box=box, fcode=freq, mJy=False)
    # ivar = load_ivar(filedb[freq]['coadd_ivar'], box=box, fcode=freq, mJy=False)
    div = divs[freq].submap(box)
    # optionally smooth the map
    if args.smooth > 0:
        imap = enmap.smooth_gauss(imap, args.smooth*u.arcmin*u.fwhm)
        div  = enmap.smooth_gauss(div, args.smooth*u.arcmin*u.fwhm)
        # weight down noise level by smoothing
        div *= sfactor(freq, args.smooth)**2
    # compute polarization angle error
    P_err = Pangle_error(imap, div, deg=True)
    # get meshgrid
    y, x = imap.posmap()/np.pi*180
    ax = axes[fi]
    cf = ax.contourf(x, y, P_err, cmap=args.cmap, levels=levels)
    ax.set_aspect(1)
    ax.set_xlim([x.max(), x.min()])
# plt.tight_layout(h_pad=-0.5)
axes[-1].set_xlabel('l [deg]')
axes[-2].set_ylabel('b [deg]')
fig.subplots_adjust(right=0.9, wspace=0.02, hspace=0.02)
cax = fig.add_axes([0.93, 0.11, 0.02, 0.77])
fig.colorbar(cf, cax=cax, orientation='vertical').set_label(r"$\delta \psi$ [deg]")
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
