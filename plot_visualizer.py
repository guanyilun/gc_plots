"""This script aims to be a general one panel plotter aimed for
producing plots for the web visualizer.

"""

from common import *
import matplotlib.pyplot as plt
import lib
import plotstyle
from matplotlib import colors

# parser defined in common
parser.add_argument("--back", default='')
parser.add_argument("--front", default='')
parser.add_argument("--title", default=None)
parser.add_argument("--figsize", default=None)
parser.add_argument("--smooth", type=float, default=None)
parser.add_argument("--colorbar", action='store_true')
parser.add_argument("--smooth2", type=float, default=None)
parser.add_argument("--dg2", type=float, default=None)
parser.add_argument("--scale", type=float, default=100)
parser.add_argument("--log", action='store_true')
parser.add_argument("--axis", action='store_true')
parser.add_argument("--mask", action='store_true')
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)
box = boxes[args.area]

# setup one panel plot
if args.figsize: figsize=eval(args.figsize)
else: figsize=None
fig = plt.figure(figsize=figsize)
if args.log:
    opts = {
        'cmap': args.cmap,
        'norm': colors.LogNorm(vmin=args.min, vmax=args.max)
    }
else:
    opts = {
        'cmap': args.cmap,
        'vmin': args.min,
        'vmax': args.max
    }
front = back = None
########
# back #
########
if args.back == 'radio':
    imap = enmap.read_map("external/meerkat/MeerKAT_radio_bubbles.fits")
    imap = enmap.submap(imap, box=box)
    imap[imap<0] = 1e-6
    back = imap
if args.back == 'herschel':
    imap = enmap.read_map("external/HIGAL_PLW1299n-00_500_RM.fits")
    imap = enmap.submap(imap, box=box)
    back = imap
if 'tmap' in args.back:
    fcode = args.back.split('_')[-1]
    if 'planck' in args.back:
        imap = load_map(filedb[fcode]['planck'], fcode=fcode, box=box)/1e9        
    else:
        imap = load_map(filedb[fcode]['coadd'], fcode=fcode, box=box)/1e9
    back = imap[0]
    if fcode != 'f220':
        opts['norm'] = colors.LogNorm(vmin=10**-0.5, vmax=10**1.5)
    else:
        opts['norm'] = colors.LogNorm(vmin=10**0.5,  vmax=10**2)
    label = texify("Total Intensity [MJy/sr]")
if 'pmap' in args.back:
    fcode = args.back.split('_')[-1]
    if 'planck' in args.back:
        imap = load_map(filedb[fcode]['planck'], fcode=fcode, box=box)/1e9        
    else:
        imap = load_map(filedb[fcode]['coadd'], fcode=fcode, box=box)/1e9
    if args.smooth:
        imap = enmap.smooth_gauss(args.smooth*u.fwhm*u.arcmin)
    back = np.sum(imap[1:]**2, axis=0)**0.5
    label = texify("Polarized Intensity [MJy/sr]")
# plotting
if back is not None:
    ax = plt.subplot(111, projection=back.wcs)
    plotstyle.setup_axis(ax, nticks=[10,5])
    im_back = ax.imshow(back, **opts)
#########
# front #
#########
if 'bvec' in args.front:
    fcode = args.front.split('_')[-1]
    imap  = load_map(filedb[fcode]['coadd'], fcode=fcode, box=box)
    ivar  = load_ivar(filedb[fcode]['coadd_ivar'], fcode=fcode, box=box)
    if args.smooth2:
        imap = enmap.smooth_gauss(imap, args.smooth2*u.fwhm*u.arcmin)
        ivar = enmap.smooth_gauss(ivar, args.smooth2*u.fwhm*u.arcmin)
        s    = sfactor(fcode, args.smooth2)
    else:
        s    = 1
    if args.dg2:
        imap = imap.downgrade(args.dg2)
        ivar = ivar.downgrade(args.dg2)
    theta  = lib.Bangle(imap[1], imap[2], toIAU=True)
    theta += (np.pi/2)  # this gets the B-field angle corrected
    Y, X   = imap.posmap()/np.pi*180
    # x- and y-components of magnetic field
    Bx, By = np.cos(theta), np.sin(theta)
    # color of vectors
    if args.mask:
        P     = np.sum(imap[1:]**2,axis=0)**0.5
        P_err = lib.P_error(imap, ivar*s**2)
        Psnr  = P / P_err
        mask  = Psnr < 3
        cmap_ = plt.get_cmap('binary')  # actual cmap doesn't matter
        color = cmap_(np.ones_like(X))
        val   = np.min([Psnr, np.ones_like(X)*3], axis=0)/3
        color[...,-1] = val
        color = color.reshape(color.shape[0]*color.shape[1],4)
    else:
        color = 'k'
    q = ax.quiver(X,Y,Bx,By,pivot='middle', headlength=0, headaxislength=0, width=0.001,
                  color=color, transform=ax.get_transform('world'), scale=args.scale)
# plotting    
if front: im_front = ax.imshow(back, **opts)
########
# axis #
########
if args.axis:
    ax.set_xlabel('$l$')
    ax.set_ylabel('$b$')
    ax.set_aspect('equal')
else:
    ax.axis('off')
plt.tight_layout()
# colorbar
if args.colorbar:
    cax = plotstyle.add_colorbar(fig, ax, size="2%")
    fig.colorbar(im_back, cax=cax).set_label(label, fontsize=12)
if args.title: plt.suptitle(texify(args.title), fontsize=16)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
