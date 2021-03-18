"""Everything about plotting"""
import matplotlib as mpl
from matplotlib import colors, cm
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

mpl.rcParams['font.size']=10
mpl.rcParams['figure.dpi']=180
mpl.rcParams['xtick.direction']='in'
mpl.rcParams['ytick.direction']='in'
mpl.rcParams['xtick.top']=True
mpl.rcParams['ytick.right']=True
mpl.rcParams['text.usetex']=True

from pixell import colorize
colorize.mpl_register('planck')

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    """build a colormap by truncating a known colormap"""
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

planck_half = truncate_colormap(plt.get_cmap('planck'), minval=0.5)
cm.register_cmap(name='planck_half', cmap=planck_half)

def add_colorbar(fig, ax, size="5%", pad=0.1, **kwargs):
    divider = make_axes_locatable(ax)
    cb = divider.append_axes('right', size=size, pad=pad, **kwargs)
    return cb
