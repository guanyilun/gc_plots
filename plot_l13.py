"""This script aims to plot l=1.3 molecular cloud complex using external data
as the contour.

"""


from common import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import lib
import plotstyle
from matplotlib import colors

mpl.rcParams['axes.formatter.useoffset'] = False
# parser defined in common
parser.add_argument("--back", default="external/HIGAL_PLW1299n-00_500_RM.fits")
parser.add_argument("--title", default="l=1.3 Complex")
parser.add_argument("--figsize", default=None)
parser.add_argument("--smooth", type=float, default=None)
parser.add_argument("--freq", default='f220')
parser.add_argument("--mask", action='store_true')
parser.add_argument("--min2", type=float, default=None)
parser.add_argument("--max2", type=float, default=None)
parser.add_argument("--scale", type=float, default=100)
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)
box = boxes[args.area]

# load coadd map
imap = load_map(filedb[args.freq]['coadd'], fcode=args.freq, box=box) / 1e9
ivar = load_ivar(filedb[args.freq]['coadd_ivar'], fcode=args.freq, box=box) * 1e18
# smooth if necessary
if args.smooth:
    imap = enmap.smooth_gauss(imap, args.smooth*u.fwhm*u.arcmin)
    ivar = enmap.smooth_gauss(ivar, args.smooth*u.fwhm*u.arcmin)
    s    = sfactor(args.freq, args.smooth)
else:
    imap = imap
    ivar = ivar
    s    = 1

# load external data
irmap = enmap.read_map(args.back)
irmap = enmap.submap(irmap, box=box)

# setup two panel plot
if args.figsize: figsize=eval(args.figsize)
else: figsize=None
Y, X = imap.posmap()/np.pi*180
opts = {
    'cmap': args.cmap,
    'vmin': args.min,
    'vmax': args.max,
}
fig = plt.figure(figsize=figsize)

# left panel: coadd + contour from herschel
ax = plt.subplot(121, projection=imap.wcs)
plotstyle.setup_axis(ax, nticks=[5,5], fmt=None)
opts.update({
    # 'vmin': args.min2,
    # 'vmax': args.max2,
    'vmin': None,
    'vmax': None,
    'norm': colors.LogNorm(vmin=30, vmax=80),
})
im = ax.imshow(imap[0], **opts)
# P  = np.sum(imap[1:]**2, axis=0)**0.5
# im = ax.imshow(P, **opts)
Y_, X_ = irmap.posmap()/np.pi*180
# contour plot
# l_ = np.percentile(np.ravel(irmap), [70,80,90])
l_ = np.percentile(np.ravel(irmap), [50,70,90])
levels = [irmap.min(), l_[0], l_[1], l_[2], irmap.max()]
ax.contour(X_, Y_, irmap, levels=levels, cmap='gray', transform=ax.get_transform('world'))
ax.set_xlabel('$l$')
ax.set_ylabel('$b$')
# add colorbar
# cax = plotstyle.add_colorbar(fig, ax)
# cbar = fig.colorbar(im, cax=cax)
cax = plotstyle.add_colorbar_hpad(ax, hpad="50%", size="2%")
cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
cbar.set_label(texify("I [MJy/sr]"), fontsize=10)
cbar.set_ticks([30, 40, 60, 80])
cbar.set_ticklabels([texify("30"), texify("40"), texify("60"), texify("80")])
cax.xaxis.set_ticks_position('top')
cax.xaxis.set_label_position('top')
ax.text(0.12, 1.02, texify("f220"), transform=ax.transAxes, fontsize=14)

# right panel
# ax = plt.subplot(122, projection=irmap.wcs)
ax = plt.subplot(122, projection=irmap.wcs)
plotstyle.setup_axis(ax, yticks=False, nticks=[5,5], fmt=None)
opts.update({
    # 'vmin': args.min2,
    # 'vmax': args.max2,
    'vmin': None,
    'vmax': None,
    'norm': colors.LogNorm(vmin=5e2, vmax=3e3),
})
im = ax.imshow(irmap, **opts)
# polarization angle plot
# reload imap to get the original resolution
theta = lib.Bangle(imap[1], imap[2], toIAU=True)
theta += (np.pi/2)  # this gets the B-field angle corrected
# x- and y-components of magnetic field
Bx = np.cos(theta)
By = np.sin(theta)

# mask by polarization intensity
if args.mask:
    P     = np.sum(imap[1:]**2,axis=0)**0.5
    P_err = lib.P_error(imap, ivar*s**2)
    Psnr  = P / P_err
    mask  = Psnr < 3
    # Pangle_err = lib.Pangle_error(imap, ivar*s**2, deg=True)    
    # mask = Pangle_err > 10 
    cmap_ = plt.get_cmap('binary')  # actual cmap doesn't matter
    color = cmap_(np.ones_like(X))
    color[ mask,-1] = 0.3
    color[~mask,-1] = 1
    color=color.reshape(color.shape[0]*color.shape[1],4)
else:
    color='k'
ax.quiver(X,Y,Bx,By,pivot='middle', headlength=0, headaxislength=0,
          color=color, transform=ax.get_transform('world'), scale=args.scale)
ax.set_xlabel('$l$')
# colorbar
# cax = plotstyle.add_colorbar(fig, ax)
# fig.colorbar(im, cax=cax).set_label(texify("Total Intensity [MJy/sr]"), fontsize=12)
cax = plotstyle.add_colorbar_hpad(ax, hpad="50%", size="2%")
cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
cbar.set_label(texify("I [MJy/sr]"), fontsize=10)
cax.xaxis.set_ticks([600, 1000, 2000, 3000])
cax.xaxis.set_ticklabels([texify("600"), texify("1000"), texify("2000"), texify("3000")], fontsize=10)
cax.xaxis.set_ticks_position('top')
cax.xaxis.set_label_position('top')
ax.text(0.10, 1.02, texify("Herschel")+" $500 \mu$"+texify("m"), transform=ax.transAxes, fontsize=10)

plt.subplots_adjust(hspace=0, wspace=0.1)
if args.title: plt.suptitle(r"$l=1.3$ {\rm Complex}", fontsize=16)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
