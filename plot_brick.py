"""This script aims to plot brick molecular clouds only using external data
as the contour.

"""

import plotstyle
from common import *
import matplotlib.pyplot as plt
# parser defined in common
parser.add_argument("--back", default="external/HIGAL_PLW0252p001_500_RM.fits")
parser.add_argument("--title", default="Brick")
parser.add_argument("--figsize", default=None)
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)
box = boxes[args.area]

# load coadd map
imap  = load_map(filedb['f220']['coadd'], fcode='f220')/1e9

# load external data
irmap = enmap.read_map(args.back)
irmap = enmap.submap(irmap, box=box)
# irmap = irmap.project(imap.shape, imap.wcs)
# setup figure and axis
imap = imap.project(irmap.shape, irmap.wcs, order=0)

# two panel plot
if args.figsize: figsize=eval(args.figsize)
else: figsize=None
fig = plt.figure(figsize=figsize)
# gridspec = fig.add_gridspec(ncols=4, width_ratios=[20,1,20,1], hspace=0) # two cols for colorbar
# gridspec = fig.add_gridspec(ncols=4, width_ratios=[20,1,20,1], hspace=0) # two cols for colorbar

# left panel: coadd + contour from herschel
# ax = plt.subplot(gridspec[0], projection=irmap.wcs)
ax = plt.subplot(121)
opts = {
    'origin': 'lower',
    'extent': box2extent(box)/np.pi*180,
    'cmap': args.cmap,
    'vmin': args.min,
    'vmax': args.max,
    'interpolation': 'nearest'
}
im = ax.imshow(imap[0], **opts)
# contour plot
ax.contour(irmap, levels=3, cmap='gray')
ax.set_xlabel('l [deg]')
ax.set_ylabel('b [deg]')
# add colorbar
cax = plotstyle.add_colorbar(fig, ax)
# cax = plt.subplot(gridspec[1], projection=None)
fig.colorbar(im, cax=cax).set_label("Total Intensity [MJy/sr]")
# fig.colorbar(im, cax=cax).set_label("Total Intensity [MJy/sr]")

# right panel
# ax = plt.subplot(gridspec[2], projection=irmap.wcs)
ax = plt.subplot(122)
im = ax.imshow(irmap, **opts)
# polarization angle plot
Y, X = imap.posmap()/np.pi*180
theta = lib.Bangle(imap[1], imap[2], toIAU=True)
# theta += (np.pi/2.)
# x- and y-components of magnetic field
Bx = np.cos(theta)
By = np.sin(theta)
ax.quiver(X,Y,Bx,By,pivot='middle', headlength=0, headaxislength=0, color='k')
# colorbar
cax = plotstyle.add_colorbar(fig, ax)
fig.colorbar(im, cax=cax).set_label("Total Intensity [MJy/sr]")
# fig.colorbar(im, cax=cax)# .set_label("Total Intensity [MJy/sr]")
ax.set_xlabel('$l$')

if args.title: plt.suptitle(args.title, fontsize=14)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile)
