import os.path as op
import matplotlib.pyplot as plt
from common import *
import plotstyle

# Note: parser defined in common
args = parser.parse_args()
box = boxes[args.area]

opts = {
    'origin': 'lower',
    'cmap': args.cmap,
    'vmin': args.min,
    'vmax': args.max
}

for fcode in ['f090','f150','f220']:
    ivar = load_ivar(filedb[fcode]['planck_ivar'], fcode=fcode, mJy=False, box=box)
    plt.imshow(ivar[0], **opts)
    plt.colorbar(shrink=0.55).set_label('ivar')
    plt.title(fcode)
    ofile = op.join(args.odir, f'ivar_{fcode}_I.png')
    print("Writing:", ofile)
    plt.savefig(ofile, bbox_inches='tight')
    plt.clf()
