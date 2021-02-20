"""This script plots the magnetic field orientation on top of some
other plots

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
# box = np.array([[-1,2],[1,-2]]) / 180*np.pi
box = np.array([[-0.5,1],[0.5,-1]]) / 180*np.pi
imap = load_map(filedb[args.freq]['coadd'], box)
ivar = load_ivar(filedb[args.freq]['coadd_ivar'], box)

P = np.sum(imap[1:]**2,axis=0)**0.5
p = P / imap[0]
# compute polarization fraction
# P_err = lib.Pangle_error(imap, ivar, deg=True)
# P_err = enmap.smooth_gauss(P_err, 1.4*arcmin/2.355)
# test downgrade
# P_err = enmap.smooth_gauss(P_err, 5*arcmin / 2.355)/10
# P_err = P_err.downgrade(10)/10

plot_opts = {
    'origin': 'lower',
    'cmap': 'planck',
    'vmin': -3,
    'vmax': -0.5,
    # 'extent': [2,-2,-1,1]
    'extent': [1,-1,-0.5,0.5]
}

fig = plt.figure()
plt.imshow(np.log10(p), **plot_opts)
plt.xlabel('l [deg]')
plt.ylabel('b [deg]')
plt.colorbar(shrink=0.5).set_label(r"$\log_{10}$P/I")
# add magnetic field
ax = fig.add_subplot(1,1,1)
imap_ds = enmap.smooth_gauss(imap, 2*arcmin/2.355).downgrade(4)
Y, X = imap_ds.posmap() / np.pi * 180
# import ipdb; ipdb.set_trace()
theta = lib.Bangle(imap_ds[1], imap_ds[2], toIAU=True)
theta += (np.pi/2.)
# x- and y-components of magnetic field
Bx = np.cos(theta)
By = np.sin(theta)
ax.quiver(X,Y,Bx,By,pivot='middle', headlength=0, headaxislength=0)
plt.tight_layout()
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
