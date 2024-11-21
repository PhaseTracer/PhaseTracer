#!/usr/bin/env python
"""
Plot phases
===========
"""

import re
import os
import sys
try:
    from StringIO import StringIO  # for Python 2
except ImportError:
    from io import StringIO  # for Python 3
from scipy.interpolate import interp1d
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc, rcParams
from matplotlib.legend_handler import HandlerPatch

ARROW_DELTA_T = 10.
QUIVER_ARROW = {"zorder": 10, "angles": 'xy', "scale_units": 'xy', "scale": 1, "units": 'dots', "width": 2.5, "minlength": 3.}
LINE = {"lw": 2.5, "alpha": 0.8}
GEV = ""
TITLE = ""
FONTSIZE = "x-small"
subcriticalTransitionArrowColour = '0.5'
criticalTransitionArrowColour = '0'

def make_legend_arrow(legend, orig_handle, xdescent, ydescent, width, height, fontsize):
#    #width = 17
#    p = mpatches.FancyArrowPatch((0.1 * width, 0.5 * height), (0.9 * width, 0.5 * height),
#                        arrowstyle="Simple,tail_width=0.5,head_width=4,head_length=6",
#                        mutation_scale=10,
#                        color="k")
##    trans = orig_handle.get_transform().transform
##    p.set_transform(trans + plt.gcf().dpi_scale_trans.inverted())
#    return p
    
    tail_length = 60
    p = mpatches.FancyArrow(0, 0.5*height, 0.2*tail_length, 0, width=2,
                            head_width=7, head_length=6, alpha=0.9, color="k")
    return p


# Now takes an array of 2-element arrays, where the first element is the legend label and the second axis is the colour
# of the arrow to draw. This function now allows multiple arrows to be added to the legend.
def add_arrows_to_leg(ax, legends_to_add):
    handles, labels = ax.get_legend_handles_labels()
    handlerMap = {}

    # The same patch function will be used for each arrow label.
    handlerPatch = HandlerPatch(patch_func=make_legend_arrow)

    for i in range(len(legends_to_add)):
        handles.append(plt.arrow(0., 0., 0., 0., color=legends_to_add[i][1]))
        labels.append(legends_to_add[i][0])
        handlerMap[handles[-1]] = handlerPatch

    #return ax.legend(handles, labels, handler_map=handlerMap, loc='upper left')
    return ax.legend(handles, labels, handler_map=handlerMap, fontsize=FONTSIZE)


def annotate_arrow(a, b, text=None, offset=10., color="k", **kwargs):
    q = plt.quiver(a[0], a[1], b[0] - a[0], b[1] - a[1], color=color, **QUIVER_ARROW)
    if text is not None:
        angle = np.degrees(np.arctan2(b[1] - a[1], b[0] - a[0]))
        xy = np.array([0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])])
        if (b[1] - a[1]) == 0:
            orthog = offset * np.array([0.5,0.5])
        else:
            orthog = offset * np.array([1., -(b[0] - a[0]) / (b[1] - a[1])]) / (1. + ((b[0] - a[0]) / (b[1] - a[1]))**2)
        ann = plt.annotate(
                text, xy=xy, xycoords='data',
                xytext=orthog, textcoords='offset points', fontsize=FONTSIZE)

        trans_angle = plt.gca().transData.transform_angles(np.array((angle,)), xy.reshape((1, 2)))[0]
        ann.set_rotation(trans_angle)


def arrow_plot(t, x, y, min_r=0.01, arrow_delta_t=10., arrow_size=8., dt=0.5, **kwargs):
    p = plt.plot(x, y, lw=0.8, zorder=10)

    # Arrows spaced evenly in temperature

    f_x = interp1d(t, x, fill_value='extrapolate')
    f_y = interp1d(t, y, fill_value='extrapolate')

    t_spaced = np.arange(max(t), min(t), -arrow_delta_t)
    x_spaced = f_x(t_spaced)
    y_spaced = f_y(t_spaced)
    dx = f_x(t_spaced - 0.5 * dt) - f_x(t_spaced + 0.5 * dt)
    dy = f_y(t_spaced - 0.5 * dt) - f_y(t_spaced + 0.5 * dt)

    # Normalize size of all arrows

    r = (dx**2 + dy**2)**0.5
    keep = r > min_r
    # Now only calculates this for the data points that we'll keep. This avoids dividing by zero when r = 0.
    dx[keep] *= arrow_size / r[keep]
    dy[keep] *= arrow_size / r[keep]

    plt.quiver(x_spaced[keep], y_spaced[keep], dx[keep], dy[keep],
               scale_units='dots', angles='xy', scale=1,
               color=p[0].get_color(), zorder=10, **kwargs)


def constant(x, y, atol=1e-1, **kwargs):
    mean_x = np.full(x.shape, np.mean(x))
    mean_y = np.full(y.shape, np.mean(y))
    constant_ = (np.isclose(mean_x, x, atol=atol).all()
                 and np.isclose(mean_y, y, atol=atol).all())
    if constant_:
        plt.plot(mean_x[0], mean_y[0], marker="*",
                 linestyle="None", ms=10, **kwargs)
    return constant_


def plane_phi_phi_one(pi, pj, phases, transitions, folderName="", fileName="plane"): #pdf_name="plane.pdf"):
    # pdf_name = output/xSM_MSbar_noZ2/0.pdf
    # Pass in folder name and file name: "output/xSM_MSbar_noZ2" and "0"
    # Then construct full name as: folderName + "/" + "phi_{0}_phi_{1}_".format(...) + fileName + ".pdf"
    #full_pdf_name = "phi_{0}_phi_{1}_".format(pi+1, pj+1)+pdf_name
    full_pdf_name = folderName + "/" + "phi_{0}_phi_{1}_".format(pi+1, pj+1) + fileName + ".pdf"
    print("Plotting phases on (phi_{}, phi_{}) plane for {} figure".format(pi+1, pj+1, full_pdf_name))
    fig, ax = plt.subplots()

    x_max = -np.inf
    x_min = np.inf
    y_max = -np.inf
    y_min = np.inf
    
    for i, p in enumerate(phases):
        t, x, y = p.T, p.phi[pi], p.phi[pj]
        x_max = max(max(x), x_max)
        x_min = min(min(x), x_min)
        y_max = max(max(y), y_max)
        y_min = min(min(y), y_min)
        label = r"Phase ${0}$ from $T = {1:.0f}$".format(i, t[0]) + GEV + " to ${0:.0f}$".format(t[-1]) + GEV
        if not constant(x, y, label=label):
            arrow_plot(t, x, y, label=label, arrow_delta_t=ARROW_DELTA_T)

    ax.set_xlabel(r"Field $x_{0}$ (GeV)".format(pi+1))
    ax.set_ylabel(r"Field $x_{0}$ (GeV)".format(pj+1))

    x_shift = max((x_max - x_min)*0.05,5)
    ax.set_xlim(x_min-x_shift, x_max+x_shift)
    y_shift = max((y_max - y_min)*0.05,5)
    ax.set_ylim(y_min-y_shift, y_max+y_shift)
        

    """if transitions:
        leg = add_arrows_to_leg(ax, "  FOPT(black arrow)")
    else:
        leg = ax.legend()"""

    extra_legends = []

    subcritical = False

    for t in transitions:
        if t.subcritical:
            subcritical = True
            break

    if transitions:
        extra_legends.append(("  FOPT(black arrow)", criticalTransitionArrowColour))

    if subcritical:
        extra_legends.append(("  ScPT(grey arrow)", subcriticalTransitionArrowColour))

    if len(extra_legends) > 0:
        leg = add_arrows_to_leg(ax, extra_legends)
    else:
        leg = ax.legend()


    leg.set_title(r"Colored arrows separated by $\Delta T = {0:.0f}$".format(ARROW_DELTA_T) + GEV, prop={"size": FONTSIZE})

    # Get size of legend
    plt.gcf().canvas.draw()
    #extent = leg.get_window_extent().inverse_transformed(ax.transAxes)
    extent = leg.get_window_extent().transformed(ax.transAxes.inverted())
    dy = extent.y1 - extent.y0

    # Add it to axis limit
    ylim = np.array(ax.get_ylim())
    ylim[1] += dy*(y_max-y_min)
    ax.set_ylim(ylim)

    prevColour = QUIVER_ARROW.get('color')

    if prevColour is not None:
        QUIVER_ARROW.pop('color')
    
    for t in transitions:
        if t.key > 0: continue
        false_vacuum = [t.false_vacuum[pi], t.false_vacuum[pj]]
        true_vacuum = [t.true_vacuum[pi], t.true_vacuum[pj]]
        dx = ((false_vacuum[1] - false_vacuum[0])**2 + (true_vacuum[1] - true_vacuum[0])**2)**0.5
        if dx > QUIVER_ARROW['minlength']:
            if t.subcritical:
                colour = subcriticalTransitionArrowColour
            else:
                colour = criticalTransitionArrowColour

            annotate_arrow(false_vacuum, true_vacuum, "${0:.0f}$".format(t.TC) + GEV, color=colour)

    if prevColour is not None:
        QUIVER_ARROW['color'] = prevColour

    plt.gcf().suptitle(TITLE, ha="right", x=0.9, y=0.925, fontsize="small")
    plt.savefig(full_pdf_name, bbox_inches="tight")


def plane_phi_phi(phases, transitions, folderName="", fileName="plane"): # pdf_name="plane.pdf"):
    n_field = phases[0].n_field
    print("Making plots with {} fields".format(n_field))
    if n_field <= 1:
        return
    for i in range(n_field):
        for j in range(i + 1, n_field):
            plane_phi_phi_one(i, j, phases, transitions, folderName, fileName)


def is_origin(x, atol=1e-1):
    return np.all(np.isclose(0, x, atol=atol))


def plane_T_V(phases, pdf_name="plane.pdf", forced_Tmin=-1.0, forced_Tmax=-1.0):
    print("Plotting potential against temperature for {} figure".format(pdf_name))
    fig, ax = plt.subplots()

    # Guess sensible axes limits

    p_max = -np.inf
    p_min = np.inf
    t_max = -np.inf
    t_min = np.inf

    for j, p in enumerate(phases):
        for i in range(p.n_field):
            field = p.phi[i]
            keep = np.logical_not(np.logical_or(np.isclose(field.min(), field, atol=1e-3), np.isclose(field.max(), field, atol=1e-3)))
            if any(keep):
                p_max = max(p.V[keep].max(), p_max)
                p_min = min(p.V[keep].min(), p_min)
                t_max = max(p.T[keep].max(), t_max)
                t_min = min(p.T[keep].min(), t_min)

    x = []
    y = []

    for j, p in enumerate(phases):
        ax.plot(p.T, p.V, label="Phase ${0}$ ".format(j), **LINE)
        x.extend(p.T)
        y.extend(p.V)

    #p_max = 0.1
    #p_min = -1.2
    #t_max = 2
    #t_min = 0

    # From https://stackoverflow.com/a/40289749 to automatically adjust the ylims based on the xlims.
    def correct_limit(ax, x, y):
        # ax: axes object handle
        #  x: data for entire x-axes
        #  y: data for entire y-axes
        # assumption: you have already set the x-limit as desired
        lims = ax.get_xlim()
        i = np.where((x > lims[0]) & (x < lims[1]))[0]  # Pretty obvious what this is doing. Find all data points in the
                                                        # domain given by xlim.
        ax.set_ylim(np.array(y)[i].min(), np.array(y)[i].max())  # Then find the min and max value of y for this
                                                                 # restricted set of data points.
        # I've made sure they're np arrays to allow for proper filtering by 'i'. Regular lists don't allow it.

    ax.set_xlabel(r"Temperature, $T$ (GeV)")
    ax.set_ylabel(r"Potential, $V(\phi(T), T)$ (GeV)${}^4$")
    #ax.set_ylim(p_min - 0.25 * abs(p_min), p_max + 0.25 * abs(p_max))
    # TODO: added to support user-defined temperature ranges for plotting.
    xlim_min = max(0., forced_Tmin)
    xlim_max = forced_Tmax if forced_Tmax > 0 else t_max

    # ax.set_xlim(0, T_max)
    ax.set_xlim(xlim_min, xlim_max)
    correct_limit(ax, x, y)
    #ax.set_xlim(0.75 * t_min, 1.25 * t_max)
    plt.legend()
    plt.gcf().suptitle(TITLE, ha="right", x=0.9, y=0.925, fontsize="small")
    plt.savefig(pdf_name, bbox_inches="tight")


def plane_phi_T(phases, transitions, pdf_name="plane.pdf", forced_Tmin=-1.0, forced_Tmax=-1.0):
    print("Plotting phases against temperature for {} figure".format(pdf_name))
    rcParams['figure.figsize'] = (15, 10)
    rcParams["text.usetex"] = True
#    rcParams['text.latex.preamble'] = [r"\usepackage{bm}"]
    n_field = phases[0].n_field
    fig, axs = plt.subplots(nrows=1, ncols=n_field)

    T_max = 0 # TODO: removed 250
    phi_min = np.inf
    phi_max = -np.inf

    if n_field == 1:
        axs = [axs]

    subcritical = False

    for i, ax in enumerate(axs):
        # TODO: added to fix bug where the axis limit of subsequent field directions would be affected by the plots of
        #  the previous field directions.
        phi_min = np.inf
        phi_max = -np.inf

        for j, p in enumerate(phases):

            # Plot phase
            phi = p.phi[i]
            T = p.T
            #ax.plot(phi, T, "-o", markersize=1, label="Phase ${0}$ ".format(j), **LINE)
            ax.plot(phi, T, label=f'$\\mathrm{{Phase}}$ \\boldmath$\\phi_{j+1}$ ', **LINE)
            ax.tick_params(size=8, labelsize=28)

            # Update minimum and maximum field values
            phi_min = min(phi_min, min(phi))
            phi_max = max(phi_max, max(phi))

            # Update maximum temperature
            """if is_origin(phi[np.argmax(T)]):
                if is_origin(phi[np.argmin(T)]):
                    T_max = max(T_max, min(T) + 100)
                    continue
                for T_, phi_ in zip(p.T, p.phi):
                    if is_origin(phi_) or True:
                        T_max = max(T_max, T_ + 100)
                        break"""

            #for T_, phi_ in zip(p.T, p.phi):
            #    T_max = max(T_max, T_ + 100)
            T_max = max(T_max, T.max()*1.1)

        prevColour = QUIVER_ARROW.get('color')

        for t in transitions:
            if t.key > 0:
                continue

            false_vacuum = t.false_vacuum[i]
            true_vacuum = t.true_vacuum[i]
            dx = true_vacuum - false_vacuum
            dy = 0.

            if t.subcritical:
                subcritical = True
                QUIVER_ARROW['color'] = subcriticalTransitionArrowColour
                ax.quiver(false_vacuum, t.TC, dx, dy, **QUIVER_ARROW)
            else:
                QUIVER_ARROW['color'] = criticalTransitionArrowColour
                ax.quiver(false_vacuum, t.TC, dx, dy, **QUIVER_ARROW)

        if prevColour is None and len(transitions) > 0:
            QUIVER_ARROW.pop('color')
        else:
            QUIVER_ARROW['color'] = prevColour

        # TODO: added to support user-defined temperature ranges for plotting.
        ylim_min = max(0., forced_Tmin)
        ylim_max = forced_Tmax if forced_Tmax > 0 else T_max

        #ax.set_ylim(0, T_max)
        ax.set_ylim(ylim_min, ylim_max)
        # TODO: changed to prevent issues where there was a lot of empty space if phi_max >> phi_max - phi_min.
        #phi_shift = max((phi_max - phi_min)*0.05, 0.05*max(abs(phi_max), abs(phi_min)))
        phi_shift = 0.05*(phi_max - phi_min)
        ax.set_xlim(phi_min-phi_shift, phi_max+phi_shift)

        if len(axs) == 1:
            #ax.set_xlabel(r"Field $x$ (GeV)")
            ax.set_xlabel('$\\mathrm{Field}$ $\\phi$ $\\mathrm{[GeV]}$', fontsize=32)
        else:
            #ax.set_xlabel(r"Field $x_{0}$ (GeV)".format(i+1))
            ax.set_xlabel('$\\mathrm{Field}$ ' + ('$\\phi_h$' if i == 0 else '$\\phi_s$') + ' $\\mathrm{[GeV]}$',
                fontsize=32)
        if i == 0:
            #ax.set_ylabel("Temperature, $T$ (GeV)")
            ax.set_ylabel('$\\mathrm{Temperature}$, $T$ $\\mathrm{[GeV]}$', fontsize=32)

    extra_legends = []

    #if transitions:
    #    extra_legends.append(("  $\\mathrm{FOPT}$", criticalTransitionArrowColour))

    if subcritical:
        extra_legends.append(("  $\\mathrm{ScPT}$", subcriticalTransitionArrowColour))

    if len(extra_legends) > 0:
        add_arrows_to_leg(axs[-1], extra_legends)

    plt.gcf().suptitle(TITLE, ha="right", x=0.9, y=0.925, fontsize="small")
    plt.savefig(pdf_name, bbox_inches="tight")


"""class Phase(object):
    def __init__(self, phase):
        self.T = phase[:, 0].T
        self.V = phase[:, 1].T
        self.phi = phase[:, 2:].T
        self.n_field = phase.shape[1] - 2

class Transition(object):
    def __init__(self, transition):
        self.n_field = (len(transition) - 5) // 2
        self.TC = transition[2]
        self.true_vacuum = transition[3:self.n_field+3]
        self.false_vacuum = transition[self.n_field+3:-2]
        self.key = transition[-2]
        self.subcritical = transition[-1]"""


class Phase(object):
    def __init__(self, key, phase):
        self.key = int(key)
        self.T = phase[:, 0].T
        self.V = phase[:, 1].T
        self.phi = phase[:, 2:].T
        self.n_field = phase.shape[1] - 2


class Transition(object):
    def __init__(self, transition):
        self.n_field = (len(transition) - 5) // 2
        self.TC = transition[2]
        self.true_vacuum = transition[3:self.n_field+3]
        self.false_vacuum = transition[self.n_field+3:-3]
        self.key = transition[-3]
        # The ID is used to match up the indices in a transition path (which is a list of transition ids).
        self.id = transition[-2]
        self.subcritical = transition[-1] > 0
        self.false_phase = int(transition[0])
        self.true_phase = int(transition[1])


class TransitionPath(object):
    def __init__(self, path):
        self.path = path


"""def load_data(dat_name):

    with open(dat_name) as d:
        data = d.read()

    data = re.sub(r'(?m)^#.*\n?', '', data)
    parts = data.split("\n\n")
    arrays = [np.genfromtxt(StringIO(p), delimiter=" ",
                            dtype=float) for p in parts]
    transitions = [Transition(a) for a in arrays if len(a.shape) != 2]
    phases = [Phase(a) for a in arrays if len(a.shape) == 2]
    return phases, transitions"""


def load_data(dat_name):

    with open(dat_name) as d:
        data = d.read()

    parts = data.split("\n\n")
    parts = [part.split("\n") for part in parts]

    phases = []
    transitions = []
    transitionPaths = []

    for part in parts:
        if len(part) < 2:
            continue

        id = part[0].split(" ")
        if id[1] == "phase":
            phases.append(constructPhase(part))
        elif id[1] == "transition":
            transitions.append(constructTransition(part))
        elif id[1] == "transition-path":
            transitionPaths.append(constructTransitionPath(part))
        else:
            raise(Exception("Invalid header for phase or transition: " + part[0]))

    return phases, transitions, transitionPaths


def constructPhase(text):
    return Phase(text[0].split()[2], np.array([[float(string) for string in text[i].split()] for i in range(1, len(text))]))


def constructTransition(text):
    return Transition(np.array([float(string) for string in text[1].split()]))


def constructTransitionPath(text):
    return TransitionPath(np.array([int(string) for string in text[1].split()]))

def get_tag():
    dir_path = os.path.dirname(os.path.realpath(__file__))
    owd = os.getcwd()
    os.chdir(dir_path)
    tag = os.popen('git describe --tags --abbrev=0').read().strip()
    os.chdir(owd)
    return tag

if __name__ == "__main__":
    np.seterr(all='raise')
    #fileName = "E://Software/PhaseTracer/"
    #fileName = "../bin/output/xSM_MSbar_noZ2/near-g4R-858/395.dat"

    try:
        tag = get_tag()
    except IOError:
        tag = None

    try:
        latex = os.environ["MATPLOTLIB_LATEX"]
    except KeyError:
        latex = False
    if latex:
        rc('text', usetex=True)
        rc('font', **{'family': 'serif', 'size': 14})
        GEV = r"\ensuremath{\,\textrm{GeV}}"
        TITLE = r"\texttt{PhaseTracer-" + tag + "}" if tag else r"\texttt{PhaseTracer}"
    else:
        GEV = r" GeV"
        TITLE = f"PhaseTracer-{tag}" if tag else "PhaseTracer"

    rc('axes', **{'grid': True})
    rc('grid', **{'ls': ':'})
    rc('legend', **{'fontsize': FONTSIZE})

    # Load data
    folderName = sys.argv[1]
    fileName = sys.argv[2]
    dat_name = "{}/{}.dat".format(folderName, fileName)
    phases, transitions, transitionPaths = load_data(dat_name)

    # TODO: added extra commandline parameters to control the plotted range of temperatures.
    Tmin = float(sys.argv[3]) if len(sys.argv) > 3 else -1.0
    Tmax = float(sys.argv[4]) if len(sys.argv) > 4 else -1.0

    # Plot phases on (phi_i, phi_j) plane
    #pdf_name = "{}/{}.pdf".format(folderName, fileName)
    plane_phi_phi(phases, transitions, folderName, fileName)

    # Plot phases on (phi_i, T) plane
    pdf_name = "{}/phi_T_{}.pdf".format(folderName, fileName)
    plane_phi_T(phases, transitions, pdf_name=pdf_name, forced_Tmin=Tmin, forced_Tmax=Tmax)

    # Plot phases on (T, V) plane
    pdf_name = "{}/V_T_{}.pdf".format(folderName, fileName)
    plane_T_V(phases, pdf_name, forced_Tmin=Tmin, forced_Tmax=Tmax)
