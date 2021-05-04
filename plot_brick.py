"""This script aims to plot brick molecular clouds only using external data
as the contour.

"""


from common import *
import matplotlib.pyplot as plt
from matplotlib import ticker
import lib
import plotstyle

# parser defined in common
parser.add_argument("--back", default="external/HIGAL_PLW0252p001_500_RM.fits")
parser.add_argument("--title", default="Brick")
parser.add_argument("--figsize", default=None)
parser.add_argument("--smooth", type=float, default=None)
parser.add_argument("--freq", default='f220')
parser.add_argument("--mask", action='store_true')
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
    # 'origin': 'lower',
    # 'extent': [X.max()+0.25/60, X.min()-0.25/60,
    #            Y.min()-0.25/60, Y.max()+0.25/60],
    'cmap': args.cmap,
    'vmin': args.min,
    'vmax': args.max,
    # 'interpolation': 'nearest'
}
fig = plt.figure(figsize=figsize)

# left panel: coadd + contour from herschel
ax = plt.subplot(121, projection=imap.wcs)
plotstyle.setup_axis(ax, nticks=[5,5], fmt=None)
im = ax.imshow(imap[0], **opts)
Y_, X_ = irmap.posmap()/np.pi*180
# contour plot
l_ = np.percentile(np.ravel(irmap), [50,70,90])
levels = [irmap.min(), l_[0], l_[1], l_[2], irmap.max()]
ax.contour(X_, Y_, irmap, levels=levels, cmap='gray', transform=ax.get_transform('world'))
ax.set_xlabel('$l$')
ax.set_ylabel('$b$')
# add colorbar
# cax = plotstyle.add_colorbar(fig, ax)
locator = ticker.MaxNLocator(nbins=5)
cax = plotstyle.add_colorbar_hpad(ax, hpad='50%', loc='top')
cb  = fig.colorbar(im, cax=cax, orientation='horizontal', ticks=locator,
                   ticklocation='top')
cb.set_label(texify("I [MJy/sr]"), fontsize=12)
ax.text(0.15, 1.03, texify("f220"), transform=ax.transAxes, fontsize=14)

# right panel
ax = plt.subplot(122, projection=irmap.wcs)
plotstyle.setup_axis(ax, yticks=False, nticks=[5,5], fmt=None)
im = ax.imshow(irmap, **opts)
# polarization angle plot
# reload imap to get the original resolution
theta = lib.Bangle(imap[1], imap[2], toIAU=True)
theta+= np.pi/2  # this gets the B-field angle corrected
# x- and y-components of magnetic field
Bx = np.cos(theta)
By = np.sin(theta)

# mask by polarization intensity
if args.mask:
    P     = np.sum(imap[1:]**2,axis=0)**0.5
    P_err = lib.P_error(imap, ivar*s**2)
    Psnr  = P / P_err
    mask  = Psnr < 3
    cmap_ = plt.get_cmap('binary')  # actual cmap doesn't matter
    color = cmap_(np.ones_like(X))
    # color[ mask,-1] = 0.3
    # color[~mask,-1] = 1
    val   = np.min([Psnr, np.ones_like(Psnr)*3], axis=0)
    val  /= 3
    color[...,-1] = val
    color=color.reshape(color.shape[0]*color.shape[1],4)
else:
    color='k'
ax.quiver(X,Y,Bx,By,pivot='middle', headlength=0, headaxislength=0,
          color=color, transform=ax.get_transform('world'))
ax.set_xlabel('$l$')
# colorbar
# cax = plotstyle.add_colorbar(fig, ax)
locator = ticker.MaxNLocator(nbins=5)
cax = plotstyle.add_colorbar_hpad(ax, hpad='50%', loc='top')
fig.colorbar(im, cax=cax, ticks=locator,
             orientation='horizontal').set_label(texify("I [MJy/sr]"),
                                                 fontsize=12)
cax.xaxis.set_label_position('top')
cax.xaxis.set_ticks_position('top')
ax.text(0.04, 1.09, texify("Herschel")+" $500 \mu$"+texify("m"), transform=ax.transAxes, fontsize=12)
ax.text(0.04, 1.02, texify("B")+"-"+texify("fields")+": "+texify("f220"), transform=ax.transAxes, fontsize=12)
# ax.text(0.06, 1.03, texify("Herschel")+" $500 \mu$"+texify("m"), transform=ax.transAxes, fontsize=12)
plt.subplots_adjust(hspace=0, wspace=0.1)
if args.title: plt.suptitle(texify(args.title), fontsize=16)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
