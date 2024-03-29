"""This script aims to plot 3lp molecular clouds using external data
as the contour.

"""


from common import *
import matplotlib.pyplot as plt
import lib
import plotstyle
from matplotlib import colors, ticker

# parser defined in common
parser.add_argument("--back", default="external/meerkat/MeerKAT_radio_bubbles.fits")
parser.add_argument("--title", default="GCRA and Sgr A")
parser.add_argument("--figsize", default=None)
parser.add_argument("--smooth", type=float, default=None)
parser.add_argument("--freq", default='f090')
parser.add_argument("--mask", action='store_true')
parser.add_argument("--min2", type=float, default=None)
parser.add_argument("--max2", type=float, default=None)
parser.add_argument("--scale", type=int, default=100)
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
fig = plt.figure(figsize=figsize)

##############
# left panel #
##############
ax = plt.subplot(121, projection=imap.wcs)
# ax = plt.subplot(121, projection=irmap.wcs)
opts = {
    # 'cmap': 'magma',
    'norm': colors.LogNorm(vmin=1e-3, vmax=1e-1),
    # 'cmap': 'planck_half',
    'cmap': 'plasma',
    # 'norm': colors.LogNorm(vmin=1e-5, vmax=1e-2),

}
plotstyle.setup_axis(ax, nticks=[5,5], fmt=None)
P  = np.sum(imap[1:]**2,axis=0)**0.5
im = ax.imshow(P/imap[0], **opts)
# im = ax.imshow(irmap, **opts)
ax.set_xlabel('$l$')
ax.set_ylabel('$b$')
cax = plotstyle.add_colorbar_hpad(ax, pad="1%", hpad="50%")
fig.colorbar(im, cax=cax, orientation='horizontal',
             shrink='50%').set_label(texify("P/I "), fontsize=12)
cax.xaxis.set_label_position('top')
cax.xaxis.set_ticks_position('top')
ax.text(0.15, 1.03, texify("f090"), transform=ax.transAxes, fontsize=14)
ax.tick_params(axis='x', colors='white', which='both', labelcolor='black')
ax.tick_params(axis='y', colors='white', which='both', labelcolor='black')
ax.set_aspect('equal')
for side in ['left','right','top','bottom']:
    ax.spines[side].set_visible(True)
    ax.spines[side].set_color('white')

# add contours
# levels = np.logspace(np.log10(5),np.log10(25),5)
# levels = [5,9,13,17,21,25]
levels = np.arange(1,31,2)
# print(levels)
# ax.contour(imap[0], levels=10, cmap='gray', norm=colors.LogNorm(vmin=5, vmax=50))
# ax.contour(imap[0], levels=levels, cmap='Greys', norm=colors.Normalize(vmin=5, vmax=50), alpha=0.6)
ax.contour(imap[0], levels=levels, colors='w', alpha=0.5)
# ax.contour(X, Y, imap[0], levels=levels, colors='k', alpha=0.3, transform=ax.get_transform('world'))
sagax, sagay = [359.94423568, -0.04616002]
ax.scatter([sagax], [sagay], s=15, marker='x', color='w',
           transform=ax.get_transform('world'), alpha=0.7)

# calculate polarization fraction
deg  = 180/np.pi
rmap = imap.modrmap(ref=[sagay/deg, (sagax-360)/deg])*deg
mask = rmap < 1/60
Pfrac= P/imap[0]
Perr = lib.P_error(imap, ivar)
Ierr = ivar[0]**-0.5
# Ierr = 0.01*imap[0]
# dfrac I + dI frac = d P ==> dfrac = (dP - dI frac)/I = I dP / I^2 - P dI / I^2
Pfrac_err = (Perr - Ierr * Pfrac)/imap[0]
print(f"mean polfrac={np.mean(Pfrac[mask]):.4f} +- {np.mean(Pfrac_err[mask])/np.sqrt(np.sum(mask)):.4f}")

###############
# right panel #
###############
ax = plt.subplot(122, projection=irmap.wcs)
opts = {
    'cmap': 'planck_half',
    'norm': colors.LogNorm(vmin=1e-5, vmax=1e-2),
}
plotstyle.setup_axis(ax, nticks=[5,5], yticks=False, fmt=None)
irmap[irmap<0] = 1e-6
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
    # color[ mask,-1] = 0.2
    # color[~mask,-1] = 1
    val   = np.min([Psnr, np.ones_like(Psnr)*3], axis=0)
    val  /= 3
    color[...,-1] = val
    color=color.reshape(color.shape[0]*color.shape[1],4)
else:
    color='k'
sagax, sagay = [359.94423568, -0.04616002]
ax.scatter([sagax], [sagay], s=15, marker='x', color='w',
           transform=ax.get_transform('world'), alpha=0.7)
ax.quiver(X,Y,Bx,By,pivot='middle', headlength=0, headaxislength=0, scale=args.scale,
          color=color, transform=ax.get_transform('world'))
ax.set_xlabel('$l$')
ax.set_ylabel('$b$')
ax.set_aspect('equal')

# colorbar
# fig.colorbar(im, cax=cax).set_label(texify("Total Intensity [MJy/sr]"), fontsize=12)
# new colorbar
cax = plotstyle.add_colorbar_hpad(ax, pad="1%", hpad="50%")
cb = fig.colorbar(im, cax=cax, orientation='horizontal',
                  ticks=[1e-5,1e-4,1e-3,1e-2], ticklocation='top')
cb.set_label(texify("I [MJy/sr]"), fontsize=12)
x_minor = ticker.LogLocator(base=10.0, subs=np.arange(1.0, 10.0)*0.1, numticks=10)
cb.ax.xaxis.set_minor_locator(x_minor)
cb.update_ticks()

ax.text(0.06, 1.07, texify("I")+": "+ texify("1.28 GHz"), transform=ax.transAxes, fontsize=12)
# ax.text(0.04, 1.07, texify("1.28 GHz"), transform=ax.transAxes, fontsize=10)
ax.text(0.06, 1.02, texify("B")+"-"+texify("fields")+": "+texify("f090"),
        transform=ax.transAxes, fontsize=12)
# ax.text(0.1, 1.03, texify("1.28 GHz"), transform=ax.transAxes, fontsize=14)
plt.subplots_adjust(hspace=0, wspace=0.1)

if args.title: plt.suptitle(texify(args.title), fontsize=16)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
