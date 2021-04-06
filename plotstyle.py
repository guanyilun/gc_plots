"""Everything about plotting"""
import matplotlib as mpl
from matplotlib import colors, cm
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.axes import Axes

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

def add_colorbar(fig, ax, size="5%", pad=0.1, axes_class=Axes, **kwargs):
    divider = make_axes_locatable(ax)
    cb = divider.append_axes('right', size=size, pad=pad, axes_class=axes_class, **kwargs)
    return cb

def setup_axis(ax, format='d.d', minor=True, nminor=[5,5], nticks=[10,10], xticks=True, yticks=True):
    ax.set_aspect('equal')
    ax.coords[0].set_format_unit('deg')
    ax.coords[0].set_ticks(number=nticks[0])
    ax.coords[1].set_format_unit('deg')
    ax.coords[1].set_ticks(number=nticks[1])
    if minor:
        ax.coords[0].display_minor_ticks(True)
        ax.coords[1].display_minor_ticks(True)
        ax.coords[0].set_minor_frequency(nminor[0])
        ax.coords[1].set_minor_frequency(nminor[1])
    if not xticks:
        ax.coords[0].set_axislabel(" ")
        ax.coords[0].set_ticklabel_visible(False)        
    if not yticks:
        ax.coords[1].set_axislabel(" ")
        ax.coords[1].set_ticklabel_visible(False)        
    # ax.coords[0].set_major_formatter(format)
    # ax.coords[1].set_major_formatter(format)
    return ax

def setup_colorbar(ax):
    ax.coords[0].set_ticks_visible(False)
    ax.coords[1].set_ticks_visible(False)
    ax.yaxis.tick_right()
    ax.coords[0].set_ticklabel_visible(False)
    # ax.coords[1].set_ticklabel_visible(False)
    ax.grid('off')
    return ax
