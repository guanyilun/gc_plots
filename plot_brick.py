"""This script aims to plot brick molecular clouds only using external data
as the contour.

"""

import plotstyle
from common import *
import matplotlib.pyplot as plt
import lib

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
imap = load_map(filedb[args.freq]['coadd'], fcode=args.freq, box=box)
ivar = load_ivar(filedb[args.freq]['coadd_ivar'], fcode=args.freq, box=box)
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
    'origin': 'lower',
    'extent': [X.max()+0.25/60, X.min()-0.25/60,
               Y.min()-0.25/60, Y.max()+0.25/60],
    'cmap': args.cmap,
    'vmin': args.min,
    'vmax': args.max,
    'interpolation': 'nearest'
}
fig = plt.figure(figsize=figsize)

# left panel: coadd + contour from herschel
ax = plt.subplot(121)
im = ax.imshow(imap[0], **opts)
Y_, X_ = irmap.posmap()/np.pi*180
# contour plot
l_ = np.percentile(np.ravel(irmap), [50,70,90])
levels = [irmap.min(), l_[0], l_[1], l_[2], irmap.max()]
ax.contour(X_, Y_, irmap, levels=levels, cmap='gray')
ax.set_xlabel('l [deg]')
ax.set_ylabel('b [deg]')
# add colorbar
cax = plotstyle.add_colorbar(fig, ax)
fig.colorbar(im, cax=cax).set_label("Total Intensity [MJy/sr]")

# right panel
ax = plt.subplot(122)
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
    color[ mask,-1] = 0.3
    color[~mask,-1] = 1
    color=color.reshape(color.shape[0]*color.shape[1],4)
else:
    color='k'
ax.quiver(X,Y,Bx,By,pivot='middle', headlength=0, headaxislength=0, color=color)
ax.set_yticklabels([])
# colorbar
cax = plotstyle.add_colorbar(fig, ax)
fig.colorbar(im, cax=cax).set_label("Total Intensity [MJy/sr]")

plt.subplots_adjust(hspace=0)
if args.title: plt.suptitle(args.title, fontsize=14)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile)
