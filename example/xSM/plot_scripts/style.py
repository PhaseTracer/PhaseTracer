"""
Common styling for plots
========================
"""

from matplotlib import rc
import seaborn as sns


def style():  
    sns.set_theme(palette="bright", style="ticks")
    rc('text', usetex=True)
    rc('text.latex', preamble=r'\usepackage{amsmath}')
    rc('font', **{'family': 'serif', 'serif': 'Computer Modern Roman', 'size': 16})
    rc('axes', **{'grid': False, 'titlesize': 16, 'labelsize': 16})
    rc('xtick', **{'labelsize': 16})
    rc('ytick', **{'labelsize': 16})
    rc('legend', **{'fontsize': 13, 'title_fontsize': 14, 'handlelength': 1., 'framealpha': 1, 'frameon': True, 'fancybox': False})


