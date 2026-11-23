"""Shared typography for window spectra and component MCMC figures."""

import matplotlib
from matplotlib.ticker import MaxNLocator


FIGURE_SIZE = (8, 10)
FONT_SIZES = {'base': 16, 'title': 24, 'axis_label': 22, 'tick': 18,
              'panel_title': 16, 'annotation': 14, 'spectrum_label': 12}
PLOT_RC = {'font.family': 'serif', 'font.size': FONT_SIZES['base'],
           'mathtext.fontset': 'dejavuserif', 'axes.linewidth': 0.8,
           'xtick.direction': 'in', 'ytick.direction': 'in',
           'xtick.top': True, 'ytick.right': True, 'pdf.fonttype': 42}


def plot_context():
    """Apply the reference profile style only while drawing one figure set."""
    return matplotlib.rc_context(PLOT_RC)


def style_axes(axis, title=None, y_bins=8):
    """Use larger labels and inward ticks on all four sides."""
    axis.xaxis.set_major_locator(MaxNLocator(nbins=7))
    axis.yaxis.set_major_locator(MaxNLocator(nbins=y_bins))
    axis.xaxis.label.set_size(FONT_SIZES['axis_label'])
    axis.yaxis.label.set_size(FONT_SIZES['axis_label'])
    axis.tick_params(direction='in', top=True, right=True,
                     labelsize=FONT_SIZES['tick'])
    if title is not None:
        axis.set_title(title, fontsize=FONT_SIZES['title'], pad=10)


def label_spectrum(axis, line_label, scale, baseline):
    """Label the actual display multiplier just above the baseline at the right."""
    axis.annotate(
        f'{line_label} ×{scale:g}', xy=(0.985, baseline),
        xycoords=axis.get_yaxis_transform(), xytext=(0, 4), textcoords='offset points',
        ha='right', va='bottom', color='black', fontsize=FONT_SIZES['spectrum_label'],
        bbox=dict(facecolor='white', edgecolor='none', pad=1, alpha=0.85),
        zorder=10, gid='spectrum_scale_label')
