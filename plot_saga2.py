"""Plot mouse PWN: v2

2 panel plot: f090 only -- brainstorm below

f090 total intensity, with magnetic field orientation plotted on top
f090 polarization fraction, with contours showing polarization fraction
 
"""

from common import *
import matplotlib.pyplot as plt
import plotstyle, lib
from matplotlib import colors, ticker

# parser defined in common
parser.add_argument("--tmin", type=float, default=None)
parser.add_argument("--tmax", type=float, default=None)
parser.add_argument("--pmin", type=float, default=None)
parser.add_argument("--pmax", type=float, default=None)
parser.add_argument("--smooth", type=float, default=None)
parser.add_argument("--radio", default='external/meerkat/MeerKAT_radio_bubbles.fits')
parser.add_argument("--freq", default='f090')
parser.add_argument("--title", default="GCRA and Sgr A")
parser.add_argument("--mask", action='store_true')
parser.add_argument("--figsize", default=None)
parser.add_argument("--scale", type=float, default=100)
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(argsn.odir)
if args.figsize: figsize = eval(args.figsize)
else: figsize = None
# define box of interests
box = boxes[args.area]

# load map
imap = load_map(filedb[args.freq]['coadd'], fcode=args.freq, box=box) / 1e9 # MJy/sr
ivar = load_ivar(filedb[args.freq]['coadd_ivar'], fcode=args.freq, box=box) * 1e18 # MJy/sr
# optionally smooth the map
if args.smooth:
    imap_sm = enmap.smooth_gauss(imap, args.smooth*u.fwhm*u.arcmin)
    ivar_sm = enmap.smooth_gauss(ivar, args.smooth*u.fwhm*u.arcmin)
    s       = sfactor(args.freq, args.smooth)
else:
    imap_sm = imap
    ivar_sm = ivar
    s       = 1

# initialize figure
fig = plt.figure(figsize=figsize)

##############
# left panel #
##############
ax   = plt.subplot(131, projection=imap.wcs)
ax   = plotstyle.setup_axis(ax, nticks=[5,5], fmt=None)
opts = {
    'cmap': 'planck_half',
    'interpolation': 'nearest',
    'vmin': 0,
    'vmax': 20,
    # 'norm': colors.LogNorm(vmin=5, vmax=50),    
}
# background: total intensity
im = ax.imshow(imap_sm[0], **opts)
cax = plotstyle.add_colorbar_hpad(ax, pad="1%", hpad="50%")
locator = ticker.MaxNLocator(nbins=4)
fig.colorbar(im, cax=cax, orientation='horizontal', ticks=locator).set_label(texify("I [MJy/sr]"), fontsize=10)
cax.xaxis.set_label_position('top')
cax.xaxis.set_ticks_position('top')
ax.text(0.12, 1.03, texify("f090"), transform=ax.transAxes, fontsize=12)
ax.set_xlabel(r"$l$")
ax.set_ylabel(r"$b$")
# foreground: magnetic field orientation

theta = lib.Bangle(imap_sm[1], imap_sm[2], toIAU=True)
theta+= np.pi/2  # this gets the B-field angle corrected
# x- and y-components of magnetic field
Bx = np.cos(theta)
By = np.sin(theta)
# calculate polarization error as masking
P     = np.sum(imap[1:]**2,axis=0)**0.5
P_err = lib.P_error(imap, ivar*s**2)
P_snr = P / P_err
mask  = P_snr < 3
cmap_ = plt.get_cmap('binary')  # actual cmap doesn't matter
color = cmap_(np.ones_like(P))
val   = np.min([P_snr, np.ones_like(P_snr)*3], axis=0)
val  /= 3
color[...,-1] = val
color = color.reshape(color.shape[0]*color.shape[1],4)
ax.quiver(Bx,By,pivot='middle', headlength=0, headaxislength=0, color=color, scale=args.scale)

################
# middle panel #
################
ax = plt.subplot(132, projection=imap.wcs)
ax = plotstyle.setup_axis(ax, nticks=[5,5], xticks=True, yticks=False, fmt=None)
ax.tick_params(axis='x', colors='white', which='both', labelcolor='black')
ax.tick_params(axis='y', colors='white', which='both', labelcolor='black')
for side in ['left','right','top','bottom']:
    ax.spines[side].set_visible(True)
    ax.spines[side].set_color('white')
ax.set_xlabel('$l$')
opts = {
    'cmap': 'magma',
    'norm': colors.LogNorm(vmin=1e-2, vmax=3),
    'interpolation': 'nearest'
}
P   = np.sum(imap_sm[1:]**2, axis=0)**0.5
im  = ax.imshow(P, **opts)
cax = plotstyle.add_colorbar_hpad(ax, pad="1%", hpad="50%")
locator = ticker.MaxNLocator(nbins=4)
fig.colorbar(im, cax=cax, orientation='horizontal', ticks=locator).set_label(texify("P [MJy/sr]"), fontsize=10)
cax.xaxis.set_label_position('top')
cax.xaxis.set_ticks_position('top')
ax.text(0.12, 1.03, texify("f090"), transform=ax.transAxes, fontsize=12)

###############
# right panel #
###############
# load radio data
irmap = enmap.read_map(args.radio, box=box)
ax = plt.subplot(133, projection=irmap.wcs)
ax = plotstyle.setup_axis(ax, nticks=[5,5], xticks=True, yticks=False, fmt=None)
opts = {
    'cmap': 'planck_half',
    'norm': colors.LogNorm(vmin=1e-5, vmax=1e-2),    
}
# fixing irmaps negative
irmap[irmap<0] = 1e-6
im = ax.imshow(irmap, **opts)
theta = lib.Bangle(imap_sm[1], imap_sm[2], toIAU=True)
theta += (np.pi/2)  # this gets the B-field angle corrected
# x- and y-components of magnetic field
Bx = np.cos(theta)
By = np.sin(theta)
Y, X = imap.posmap() / np.pi * 180
if args.mask:
    P     = np.sum(imap[1:]**2,axis=0)**0.5
    P_err = lib.P_error(imap, ivar*s**2)
    Psnr  = P / P_err
    mask  = Psnr < 3
    # Pangle_err = lib.Pangle_error(imap, ivar*s**2, deg=True)    
    # mask = Pangle_err > 10 
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
          scale=args.scale, color=color, transform=ax.get_transform('world'))
ax.set_xlabel('$l$')
ax.set_ylabel('$b$')
ax.set_aspect('equal')
# cax = plotstyle.add_colorbar(fig, ax)
# fig.colorbar(im, cax=cax).set_label(texify("Total Intensity [MJy/sr]"), fontsize=10)
cax = plotstyle.add_colorbar_hpad(ax, pad="1%", hpad="50%")
fig.colorbar(im, cax=cax, orientation='horizontal').set_label(texify("I [MJy/sr]"), fontsize=10)
cax.xaxis.set_label_position('top')
cax.xaxis.set_ticks_position('top')
ax.text(0.04, 1.07, texify("I")+": "+ texify("1.28 GHz"), transform=ax.transAxes, fontsize=8)
# ax.text(0.04, 1.07, texify("1.28 GHz"), transform=ax.transAxes, fontsize=10)
ax.text(0.04, 1.02, texify("B")+"-"+texify("fields")+": "+texify("f090"),
        transform=ax.transAxes, fontsize=8)
if args.title: plt.suptitle(texify(args.title), fontsize=16)
fig.subplots_adjust(hspace=0, wspace=0.1)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
