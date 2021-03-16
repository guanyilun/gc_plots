"""This script aims to calculate angular dispersion defined as

S^2 = \sum_N (psi(x) - psi(x+d))^2  where psi_i are taken in a fixed annulus d
"""

import matplotlib.pyplot as plt
from common import *
import plotstyle

# parser defined in common
parser.add_argument("--freq", default='f090')
parser.add_argument("--delta", help='delta in annulus in arcmin', type=float, default=2)
parser.add_argument("--ddelta", help='d delta of in arcmin', type=float, default=0.5)
parser.add_argument("--smooth", type=float, default=None)
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)

# load map
box = boxes[args.area]
imap = load_map(filedb[args.freq]['coadd'], fcode=args.freq, box=box)
if args.smooth:
    imap = enmap.smooth_gauss(imap, args.smooth*u.fwhm*u.arcmin)
# build an annulus map
rmap = imap.modrmap()
delta, ddelta = args.delta*u.arcmin, args.ddelta*u.arcmin
annulus = np.abs(rmap-delta) <= ddelta/2
print(np.sum(annulus))
annulus = annulus.astype(float)

# S^2 = <psi_i^2>+<psi>^2-2 psi_i psi
psi    = lib.Bangle(imap[1], imap[2], toIAU=True)/np.pi*180
psi2   = psi**2
# psi_i  = enmap.ifft(enmap.fft(enmap.fftshift(annulus))*enmap.fft(psi)).real
# psi_i2 = enmap.ifft(enmap.fft(enmap.fftshift(annulus))*enmap.fft(psi2)).real
# calculate dispersion
# n = np.sum(annulus)
# S2 = (psi_i2+psi2-2*psi_i*psi)/np.sum(annulus)
# S2 = (psi_i2/n-(psi_i/n)**2)

# method 2: manual convolve
k = np.array([
    [0,0,0,0,0],
    [0,1,1,1,0],
    [0,1,1,1,0],
    [0,1,1,1,0],
    [0,0,0,0,0],
])
from scipy.signal import convolve2d
psi_i  = convolve2d(psi, k, mode='same')
psi_i2 = convolve2d(psi**2, k, mode='same')
n = np.sum(k)
S2 = (psi_i2/n-(psi_i/n)**2)
opts = {
    'origin': 'lower',
    'extent': box2extent(box)/np.pi*180,
    'cmap': args.cmap,
    'vmin': args.min,
    'vmax': args.max
}
plt.imshow(np.sqrt(S2), **opts)
# plt.imshow(enmap.ifft(enmap.fft(enmap.fftshift(annulus))*enmap.fft(psi**2)).real, **opts)
plt.xlabel('l [deg]')
plt.ylabel('b [deg]')
plt.colorbar(shrink=0.55).set_label(r'$S$ [deg]')
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
