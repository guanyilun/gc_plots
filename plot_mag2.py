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

def parse_expr(expr):
    if expr:
        return [float(f) for f in expr.split(',')]
    else:
        return [None]*3

parser = argparse.ArgumentParser()
parser.add_argument("-o","--odir", default="plots")
parser.add_argument("--freq",default=None)
parser.add_argument("--oname", default="polfrac.pdf")
parser.add_argument("--box", default=None)
parser.add_argument("--underlay", default='T')
parser.add_argument("--cmap", default='planck')
parser.add_argument("--min", type=str, default=None)
parser.add_argument("--max", type=str, default=None)
parser.add_argument("--downgrade", help="magnetic field downgrade", type=int, default=1)
parser.add_argument("--area", default='quat')
parser.add_argument("--contour", action='store_true')
parser.add_argument("--smooth", type=float, default=None)
parser.add_argument("--figsize", default=None)
parser.add_argument("--mask", help='mask Pangle err > x degree', type=float, default=None)
parser.add_argument("--color", default='white')
parser.add_argument("--alpha", type=float, default=1)
parser.add_argument("--largefont", type=int, default=None)
parser.add_argument("--scale", type=int, default=None)
parser.add_argument("--title", default=None)
parser.add_argument("--transpose", action='store_true')
parser.add_argument("--show-freq", action='store_true')
parser.add_argument("--sep-colorbar", action='store_true')
parser.add_argument("--axdir", default='in')
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)
if args.figsize: figsize=eval(args.figsize)
else: figsize=None
if args.largefont: mpl.rcParams['font.size'] = args.largefont
if args.axdir == 'out':
    mpl.rcParams['xtick.direction'] = 'out'
    mpl.rcParams['ytick.direction'] = 'out'
    mpl.rcParams['xtick.top']       = False
    mpl.rcParams['ytick.right']     = False   
    
# define box of interests
if args.box is None:
    box = boxes[args.area]
else:
    box = np.array(eval(args.box)) / 180*np.pi
# frequency of interests
if args.freq: freqs = [args.freq]
else: freqs = ['f090','f150','f220']

maxs = parse_expr(args.max)
mins = parse_expr(args.min)

# define figure
if not args.transpose:
    nrow, ncol = len(freqs), 1
else:
    nrow, ncol = 1, len(freqs)
# temp loading to get wcs
imap = load_map(filedb['f090']['coadd'], fcode='f090', box=box)
fig, axes = plt.subplots(nrow, ncol, figsize=figsize, subplot_kw={'projection': imap.wcs})
for i, freq in enumerate(freqs):
    norm = Normalize(vmin=mins[i], vmax=maxs[i])
    # load map and ivar
    imap = load_map(filedb[freq]['coadd'], fcode=freq, box=box)
    ivar = load_ivar(filedb[freq]['coadd_ivar'], fcode=freq, box=box)
    if len(freqs) == 1: ax = axes
    else: ax = axes[i]
    if len(freqs) == 3:
        if i != 2:
            plotstyle.setup_axis(ax, nticks=[10,5], yticks=True, xticks=False)
        else:
            plotstyle.setup_axis(ax, nticks=[10,5], yticks=True, xticks=True)
        if args.axdir == 'out':
            ax.coords[0].set_ticks_position('b')
            ax.coords[1].set_ticks_position('l')
    else:
        plotstyle.setup_axis(ax, nticks=[5,5], fmt=None)
    # define underlay to plot under magnetic field orientation
    if args.underlay == 'T':
        back  = imap[0]/1e9
        label = 'Total Intensity [MJy/sr]'
    elif args.underlay == 'P':
        P     = np.sum(imap[1:]**2,axis=0)**0.5/1e9
        back  = enmap.smooth_gauss(P, 1*u.arcmin*u.fwhm)
        label = texify('Polarized Intensity [MJy/sr]')
    elif args.underlay == 'plog':
        P     = np.sum(imap[1:]**2,axis=0)**0.5
        p     = P / imap[0]
        back  = np.log10(p)
        label = r"$\log_{10}$P/I"
    elif args.underlay == 'plin':
        P     = np.sum(imap[1:]**2,axis=0)**0.5
        p     = P / imap[0]
        back  = p
        label = texify("Polarization Fraction")
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
    # Y, X = imap_ds.posmap()
    if not args.contour:
        plot_opts = {
            'origin': 'lower',
            'cmap': args.cmap,
            'vmin': mins[i],
            'vmax': maxs[i],
        }
        im = ax.imshow(back, **plot_opts)
    else:
        im = ax.contourf(back, colors=color, levels=level)
    theta = lib.Bangle(imap_ds[1], imap_ds[2], toIAU=True)
    theta += (np.pi/2.)
    # x- and y-components of magnetic field
    Bx = np.cos(theta)
    By = np.sin(theta)
    # also compute the angle uncertainty
    P_err = lib.Pangle_error(imap_ds, ivar_ds*s**2, deg=True)
    # optionally mask regions above certain angular uncertainty level
    # and high the segments through alpha
    # get meshgrid to plot
    Y, X = imap_ds.posmap()/np.pi*180
    if args.mask:
        mask = P_err > args.mask
        cmap_ = plt.get_cmap('gray')  # actual cmap doesn't matter
        color = cmap_(np.ones_like(X))
        color[mask,-1] = 0
        color[~mask,-1] = args.alpha
        color=color.reshape(color.shape[0]*color.shape[1],4)
        ax.quiver(X,Y,Bx,By,pivot='middle', headlength=0,
                  headaxislength=0, color=color, scale=args.scale,
                  transform=ax.get_transform('world'))
    else:
        color = args.color
        ax.quiver(X,Y,Bx,By,pivot='middle', headlength=0,
                  transform=ax.get_transform('world'), headaxislength=0, color=color,
                  alpha=args.alpha, scale=args.scale)
    ax.set_aspect('equal')
    if args.sep_colorbar:
        cbar = plotstyle.add_colorbar(fig, ax, size="5%")
        if i == len(freqs)-1: fig.colorbar(im, cax=cbar).set_label(label, fontsize=14)
        else: fig.colorbar(im, cax=cbar)
    if args.show_freq: ax.set_title(texify(freqs[i]), fontsize=14)
    
if len(freqs) == 1:
    axes.set_xlabel('$l$')
    axes.set_ylabel('$b$')
    # axes.set_title(texify(args.freq))
else:
    if not args.transpose:
        axes[-1].set_xlabel('$l$')
        axes[-1].set_ylabel('$b$')
        axes[0].set_ylabel('$b$')
        axes[1].set_ylabel('$b$')
    else:
        axes[0].set_ylabel('$b$')
        axes[0].set_xlabel('$l$')
        axes[-1].set_xticklabels([])
        axes[-2].set_xticklabels([])
        axes[-1].set_yticklabels([])
        axes[-2].set_yticklabels([])
    for ax in axes:
        # for side in ['left','right','top','bottom']:
        #     ax.spines[side].set_visible(True)
        #     ax.spines[side].set_color('red')
        ax.tick_params(axis='x', colors='black', which='both', labelcolor='black')
        ax.tick_params(axis='y', colors='black', which='both', labelcolor='black')
if not args.sep_colorbar:
    fig.subplots_adjust(right=0.9, hspace=0.03, wspace=0.02)
    # cax = fig.add_axes([0.93, 0.11, 0.04, 0.77])
    # cax = fig.add_axes([0.92, 0.12, 0.02, 0.74])
    cax = fig.add_axes([0.92, 0.3, 0.02, 0.4])
    fig.colorbar(im, cax=cax, orientation='vertical').set_label(label)    
else:
    fig.subplots_adjust(right=0.9, hspace=0.02)

if args.title: plt.suptitle(texify(args.title), fontsize=16)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
