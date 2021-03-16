"""This script plots the magnetic field orientation on top of some
other plots

update: plot three freqs together
"""
import argparse, os, os.path as op
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
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
parser.add_argument("--downgrade", help="magnetic field downgrade", type=int, default=1)
parser.add_argument("--area", default='quat')
parser.add_argument("--contour", action='store_true')
parser.add_argument("--smooth", type=float, default=None)
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)
norm  = Normalize(vmin=args.min, vmax=args.max)

# define box of interests
if args.box is None:
    box = boxes[args.area]
else:
    box = np.array(eval(args.box)) / 180*np.pi

# define figure
fig, axes = plt.subplots(3,1,sharex=True, figsize=(6,9))
freqs = ['f090','f150','f220']
for i, freq in enumerate(freqs):
    # load map
    imap = load_map(filedb[freq]['coadd'], fcode=freq, box=box) 
    # define underlay to plot under magnetic field orientation
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
        # for contour plot
        level = np.array([0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2,0.3])
        color = plt.get_cmap(args.cmap)(norm(level**0.5))
    elif args.underlay == 'Perr':
        ivar    = load_ivar(filedb[freq]['coadd_ivar'], fcode=freq, box=box)
        # not accurate but close
        imap_   = enmap.smooth_gauss(imap, args.downgrade*0.5*u.fwhm*u.arcmin)
        ivar_   = enmap.smooth_gauss(ivar, args.downgrade*0.5*u.fwhm*u.arcmin)
        s       = sfactor(freq, 0.5*args.downgrade)
        P_err   = lib.Pangle_error(imap_, ivar_*s**2, deg=True)
        back    = P_err
        label   = r'$\delta \psi$ [deg]'
        level   = [0,5,10,15,30,45,60]
        color   = plt.get_cmap(args.cmap)(norm(level))
    # smooth if necessary
    if args.smooth:
        imap_ds = enmap.smooth_gauss(imap, args.smooth*u.arcmin*u.fwhm)
    else:
        imap_ds = imap
    # doengrade if necessary
    if args.downgrade > 1:
        imap_ds = imap_ds.downgrade(args.downgrade)
    # get meshgrid to plot
    Y, X = imap_ds.posmap()/np.pi*180
    if not args.contour:
        plot_opts = {
            'origin': 'lower',
            'cmap': args.cmap,
            'vmin': args.min,
            'vmax': args.max,
            'extent': box2extent(box)/np.pi*180
        }
        im = axes[i].imshow(back, **plot_opts)
    else:
        Y_, X_ = imap.posmap()/np.pi*180
        im = axes[i].contourf(X_, Y_, back, colors=color, levels=level)
        xmin, xmax = box2extent(box)[:2]/np.pi*180
        axes[i].set_xlim([xmin, xmax])  # revert x axis
    theta = lib.Bangle(imap_ds[1], imap_ds[2], toIAU=True)
    theta += (np.pi/2.)
    # x- and y-components of magnetic field
    Bx = np.cos(theta)
    By = np.sin(theta)
    if args.cmap == 'magma': color = 'w'
    else: color='k'
    axes[i].quiver(X,Y,Bx,By,pivot='middle', headlength=0, headaxislength=0, color=color, alpha=0.7)
axes[-1].set_xlabel('l [deg]')
axes[-2].set_ylabel('b [deg]')
fig.subplots_adjust(right=0.9, hspace=0.02, wspace=0.02)
cax = fig.add_axes([0.93,0.11, 0.02,0.77])
fig.colorbar(im, cax=cax, orientation='vertical').set_label(label)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
