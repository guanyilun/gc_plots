"""Everything about plotting"""
import matplotlib as mpl
from matplotlib import colors, cm
import matplotlib.pyplot as plt
import numpy as np

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
