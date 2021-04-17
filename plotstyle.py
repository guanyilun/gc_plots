"""Everything about plotting"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import colors, cm
from matplotlib.axes import Axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.axes_divider import AxesDivider
import mpl_toolkits.axes_grid1.axes_size as Size


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

def add_colorbar(fig, ax, size="5%", pad=0.1, axes_class=Axes, loc='right', **kwargs):
    divider = make_axes_locatable(ax)
    cb = divider.append_axes(loc, size=size, pad=pad, axes_class=axes_class, **kwargs)
    return cb

def add_colorbar_hpad(fig, ax, size="5%", axes_class=Axes, **kwargs):
    divider = MyAxesDivider(ax)
    locator = divider.new_locator(nx=0, ny=0)
    ax.set_axes_locator(locator)
    cb = divider.new_vertical_pad(size, axes_class=axes_class, **kwargs)
    fig.add_axes(cb)
    return cb

def setup_axis(ax, fmt='d.d', minor=True, nminor=[5,5], nticks=[10,10], xticks=True, yticks=True):
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
    if fmt is not None:
        ax.coords[0].set_major_formatter(fmt)
        ax.coords[1].set_major_formatter(fmt)
    return ax


class MyAxesDivider(AxesDivider):
    """a more flexible AxesDivider that supports grid"""
    def __init__(self, *args, **kwargs):
        AxesDivider.__init__(self, *args, **kwargs)

    def new_vertical_pad(self, size, pad=None, hpad=None,
                         hpad_left=True, pack_start=False, **kwargs):
        """similar to new_vertical but adds a horizontal pad"""
        if pad:
            if not isinstance(pad, Size._Base):
                pad = Size.from_any(pad, fraction_ref=self._yref)
            if pack_start:
                self._vertical.insert(0, pad)
                self._yrefindex += 1
            else:
                self._vertical.append(pad)
        # new horizontal pad
        if hpad:
            if not isinstance(hpad, Size._Base):
                pct = float(hpad[:-1])
                hpad = Size.from_any(hpad, fraction_ref=self._xref)
            if hpad_left:
                self._horizontal = [hpad, Size.from_any(f"{100-pct}%", fraction_ref=self._xref)]
            else:
                self._horizontal = [Size.from_any(f"{100-pct}%", fraction_ref=self._xref), hpad]
        if not isinstance(size, Size._Base):
            size = Size.from_any(size, fraction_ref=self._yref)
        # when pack_start is True, the colorbar is added in the start
        # of the list since the origin is at bottom left, this means
        # to add the colorbar at the bottom of the figure, whereas
        # when pack_start is False, the colorbar will be added at the
        # top of the figure.
        #
        # Note that when hpad is used, the horizontal and vertical
        # spliting will look like: (assuming pack_start is False,
        # hpad=50%, size=25%)
        #
        # xxxxx|xxxxx  --> top row: hpad, colorbar
        # -----|-----
        # xxxxx|xxxxx  --> bottom row: image
        # xxxxx|xxxxx
        # xxxxx|xxxxx
        # xxxxx|xxxxx
        #
        # for now I will change when the colorbar is at the top, to be
        # directly usable, and leave the equivalent work to the other
        # case (bottom) as future work
        if pack_start:
            self._vertical.insert(0, size)
            self._yrefindex += 1
            locator = self.new_locator(nx=self._xrefindex, ny=0)
        else:
            # colorbar is added on top
            self._vertical.append(size)
            # if no hpad is specified, reduce to the original behavior
            if not hpad:
                locator = self.new_locator(
                    nx=self._xrefindex, ny=len(self._vertical)-1)
            else:
                # where the colorbar is located depends on whether
                # hpad is added on the left or right
                if hpad_left: nx = 1
                else: nx = 0
                locator = self.new_locator(nx=nx, ny=len(self._vertical)-1)
                # now the original axes becomes grid-based, we also need
                # to make sure the original plot spans two columns.
                # this is achieved by setting nx:nx1 --> 0:2
                locator_orig = self.new_locator(nx=0, nx1=2, ny=0)
                self._axes.set_axes_locator(locator_orig)
        ax = self._get_new_axes(**kwargs)
        ax.set_axes_locator(locator)
        return ax
