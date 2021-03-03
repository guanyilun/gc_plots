"""This script plots the magnetic field orientation on top of some
other plots

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
parser.add_argument("--box", default=None)
parser.add_argument("--underlay", default='T')
parser.add_argument("--cmap", default='planck')
parser.add_argument("--min", type=float, default=None)
parser.add_argument("--max", type=float, default=None)
parser.add_argument("--downgrade", help="magnetic field downgrade", type=int, default=None)
parser.add_argument("--area", default='quat')
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)

if args.box is None:
    box = boxes[args.area]
else:
    box = np.array(eval(args.box)) / 180*np.pi

imap = load_map(filedb[args.freq]['coadd'], fcode=args.freq, box=box)

if args.underlay == 'T':
    back  = imap[0]
    label = 'Temperature'
elif args.underlay == 'P':
    P     = np.sum(imap[1:]**2,axis=0)**0.5
    back  = enmap.smooth_gauss(P, 1*u.arcmin*u.fwhm)
    label = 'Polarization Intensity'
elif args.underlay == 'plog':
    P     = np.sum(imap[1:]**2,axis=0)**0.5
    p     = P / imap[0]
    back  = np.log10(p)
    label = r"$\log_{10}$P/I"
elif args.underlay == 'plin':
    P     = np.sum(imap[1:]**2,axis=0)**0.5
    p     = P / imap[0]
    back  = p
    label = r"P/I"
if args.downgrade is None:
    imap_ds = imap
else:
    # smooth before downgrade
    imap_ds = enmap.smooth_gauss(
        imap, args.downgrade*0.5*arcmin/2.355).downgrade(args.downgrade)
    
Y, X = imap_ds.posmap()
plot_opts = {
    'origin': 'lower',
    'cmap': args.cmap,
    'vmin': args.min,
    'vmax': args.max,
    'extent': [X.max()+0.25*u.arcmin, X.min()-0.25*u.arcmin,
               Y.min()-0.25*u.arcmin, Y.max()+0.25*u.arcmin]
}
fig = plt.figure()
plt.imshow(back, **plot_opts)
plt.xlabel('l [deg]')
plt.ylabel('b [deg]')
plt.colorbar(shrink=0.48).set_label(label)
# add magnetic field
ax = fig.add_subplot(1,1,1)

theta = lib.Bangle(imap_ds[1], imap_ds[2], toIAU=True)
theta += (np.pi/2.)
# x- and y-components of magnetic field
Bx = np.cos(theta)
By = np.sin(theta)
ax.quiver(X,Y,Bx,By,pivot='middle', headlength=0, headaxislength=0)
plt.tight_layout()
plt.axis('off')
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
