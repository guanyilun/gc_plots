"""This script plots the magnetic field orientation on top of some
other plots

update: plot three freqs together
"""
import argparse, os, os.path as op
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import plotstyle
from common import *
import lib
from pixell import utils as u

parser = argparse.ArgumentParser()
parser.add_argument("-o","--odir", default="plots")
parser.add_argument("--freq",default=None)
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
parser.add_argument("--figsize", default=None)
parser.add_argument("--mask", help='mask Pangle err > x degree', type=float, default=None)
parser.add_argument("--color", default='white')
parser.add_argument("--alpha", type=float, default=1)
parser.add_argument("--largefont", action="store_true")
parser.add_argument("--scale", type=int, default=None)
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)
if args.figsize: figsize=eval(args.figsize)
else: figsize=None
if args.largefont: mpl.rcParams['font.size'] = 18
norm = Normalize(vmin=args.min, vmax=args.max)

# define box of interests
if args.box is None:
    box = boxes[args.area]
else:
    box = np.array(eval(args.box)) / 180*np.pi
# frequency of interests
if args.freq: freqs = [args.freq]
else: freqs = ['f090','f150','f220']
# define figure
fig, axes = plt.subplots(len(freqs),1,sharex=True, figsize=figsize)
for i, freq in enumerate(freqs):
    # load map and ivar
    imap = load_map(filedb[freq]['coadd'], fcode=freq, box=box)
    ivar = load_ivar(filedb[freq]['coadd_ivar'], fcode=freq, box=box)
    if len(freqs) == 1: ax = axes
    else: ax = axes[i]
    # define underlay to plot under magnetic field orientation
    if args.underlay == 'T':
        back  = imap[0]
        label = 'Temperature'
    elif args.underlay == 'P':
        P     = np.sum(imap[1:]**2,axis=0)**0.5/1e9
        back  = enmap.smooth_gauss(P, 1*u.arcmin*u.fwhm)
        label = 'Polarization Intensity [MJy/sr]'
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
        imap_ds = enmap.smooth_gauss(imap, args.smooth*u.fwhm*u.arcmin)
        ivar_ds = enmap.smooth_gauss(ivar, args.smooth*u.fwhm*u.arcmin)
        s       = sfactor(freq, args.smooth)
    else:
        imap_ds = imap
        ivar_ds = ivar
        s       = 1
    # downgrade if necessary
    if args.downgrade > 1:
        imap_ds = imap_ds.downgrade(args.downgrade)
        ivar_ds = ivar_ds.downgrade(args.downgrade)
    # get meshgrid to plot
    Y, X = imap_ds.posmap()/np.pi*180
    if not args.contour:
        plot_opts = {
            'origin': 'lower',
            'cmap': args.cmap,
            'vmin': args.min,
            'vmax': args.max,
            'extent': [X.max()+0.25/60, X.min()-0.25/60,
                       Y.min()-0.25/60, Y.max()+0.25/60]
        }
        im = ax.imshow(back, **plot_opts)
    else:
        Y_, X_ = imap.posmap()/np.pi*180
        im = ax.contourf(X_, Y_, back, colors=color, levels=level)
        xmin, xmax = box2extent(box)[:2]/np.pi*180
        ax.set_xlim([xmin, xmax])  # revert x axis
    theta = lib.Bangle(imap_ds[1], imap_ds[2], toIAU=True)
    theta += (np.pi/2.)
    # x- and y-components of magnetic field
    Bx = np.cos(theta)
    By = np.sin(theta)
    # also compute the angle uncertainty
    P_err = lib.Pangle_error(imap_ds, ivar_ds*s**2, deg=True)
    # optionally mask regions above certain angular uncertainty level
    # and high the segments through alpha
    if args.mask:
        mask = P_err > args.mask
        cmap_ = plt.get_cmap('gray')  # actual cmap doesn't matter
        color = cmap_(np.ones_like(X))
        color[mask,-1] = 0
        color[~mask,-1] = args.alpha
        color=color.reshape(color.shape[0]*color.shape[1],4)
        ax.quiver(X,Y,Bx,By,pivot='middle', headlength=0, headaxislength=0, color=color, scale=args.scale)
    else:
        color = args.color
        ax.quiver(X,Y,Bx,By,pivot='middle', headlength=0,
                  headaxislength=0, color=color, alpha=args.alpha, scale=args.scale)
    ax.set_aspect('equal')
if len(freqs) == 1:
    axes.set_xlabel('l [deg]')
    axes.set_ylabel('b [deg]')
else:
    axes[-1].set_xlabel('l [deg]')
    axes[-2].set_ylabel('b [deg]')
fig.subplots_adjust(right=0.9, hspace=0.02, wspace=0.02)
cax = fig.add_axes([0.93,0.11, 0.04,0.77])
fig.colorbar(im, cax=cax, orientation='vertical').set_label(label)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
