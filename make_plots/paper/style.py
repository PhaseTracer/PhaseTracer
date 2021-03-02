from matplotlib import rc


def style():  
    rc('text', usetex=True)
    rc('font', **{'family': 'serif', 'serif': 'Computer Modern Roman', 'size': 16})
    rc('axes', **{'grid': False, 'titlesize': 14, 'labelsize': 16})
    rc('legend', **{'fontsize': 18, 'title_fontsize': 18, 'handlelength': 1., 'frameon': True})
