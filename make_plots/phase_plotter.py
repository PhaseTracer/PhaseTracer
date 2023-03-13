#!/usr/bin/env python
"""
Plot phases
===========
"""

import re
import os
import sys
from io import StringIO
from scipy.interpolate import interp1d
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
from matplotlib.legend_handler import HandlerPatch

ARROW_DELTA_T = 10.
QUIVER_ARROW = {"zorder": 10, "angles": 'xy', "scale_units": 'xy', "scale": 1, "units": 'dots', "width": 1.5, "minlength": 3.}
LINE = {"lw": 2, "alpha": 0.8}
GEV = ""
FONTSIZE = "x-small"

def make_legend_arrow(legend, orig_handle, xdescent, ydescent, width, height, fontsize):
    width = 17
    p = mpatches.FancyArrow(0, 0.5*height, width, 0,
                            head_width=5, head_length=4, alpha=0.9, color="k")
    return p


def add_arrow_to_leg(ax, label):
    handles, labels = ax.get_legend_handles_labels()
    arrow = plt.arrow(0., 0., 0., 0., color='k')
    return ax.legend(handles + [arrow], labels + [label],
              handler_map={arrow: HandlerPatch(
                patch_func=make_legend_arrow)},
              loc='upper right')


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
    dx *= arrow_size / r
    dy *= arrow_size / r

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


def plane_phi_phi_one(pi, pj, phases, transitions, pdf_name="plane.pdf"):
    full_pdf_name = "phi_{0}_phi_{1}_".format(pi+1, pj+1)+pdf_name
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
        

    if transitions:
        leg = add_arrow_to_leg(ax, "  FOPT(black arrow)")
    else:
        leg = ax.legend()

    leg.set_title(r"Colored arrows separated by $\Delta T = {0:.0f}$".format(ARROW_DELTA_T) + GEV, prop={"size": FONTSIZE})

    # Get size of legend
    plt.gcf().canvas.draw()
    extent = leg.get_window_extent().transformed(ax.transAxes.inverted())
    dy = extent.y1 - extent.y0

    # Add it to axis limit
    ylim = np.array(ax.get_ylim())
    ylim[1] += dy*(y_max-y_min)
    ax.set_ylim(ylim)
    
    for t in transitions:
        if t.key >0 : continue
        false_vacuum = [t.false_vacuum[pi], t.false_vacuum[pj]]
        true_vacuum = [t.true_vacuum[pi], t.true_vacuum[pj]]
        dx = ((false_vacuum[1] - false_vacuum[0])**2 + (true_vacuum[1] - true_vacuum[0])**2)**0.5
        if dx > QUIVER_ARROW['minlength']:
          annotate_arrow(false_vacuum, true_vacuum, "${0:.0f}$".format(t.TC) + GEV)

    plt.savefig(full_pdf_name, bbox_inches="tight")


def plane_phi_phi(phases, transitions, pdf_name="plane.pdf"):
    n_field = phases[0].n_field
    print("Making plots with {} fields".format(n_field))
    if n_field <= 1:
        return
    for i in range(n_field):
        for j in range(i + 1, n_field):
            plane_phi_phi_one(i, j, phases, transitions, pdf_name)


def is_origin(x, atol=1e-1):
    return np.all(np.isclose(0, x, atol=atol))


def plane_T_V(phases, pdf_name="plane.pdf"):
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

    for j, p in enumerate(phases):
        ax.plot(p.T, p.V, label="Phase ${0}$ ".format(j), **LINE)

    ax.set_xlabel(r"Temperature, $T$ (GeV)")
    ax.set_ylabel(r"Potential, $V(\phi(T), T)$ (GeV)${}^4$")
    ax.set_ylim(p_min - 0.25 * abs(p_min), p_max + 0.25 * abs(p_max))
    ax.set_xlim(0.75 * t_min, 1.25 * t_max)
    plt.legend()
    plt.savefig(pdf_name, bbox_inches="tight")


def plane_phi_T(phases, transitions, pdf_name="plane.pdf"):
    print("Plotting phases against temperature for {} figure".format(pdf_name))
    n_field = phases[0].n_field
    fig, axs = plt.subplots(nrows=1, ncols=n_field)

    T_max = 250
    phi_min = np.inf
    phi_max = -np.inf

    if n_field == 1:
        axs = [axs]

    for i, ax in enumerate(axs):
        for j, p in enumerate(phases):

            # Plot phase
            phi = p.phi[i]
            T = p.T
            ax.plot(phi, T, "-o", markersize=1, label="Phase ${0}$ ".format(j), **LINE)

            # Update minimum and maximum field values
            phi_min = min(phi_min, min(phi))
            phi_max = max(phi_max, max(phi))

            # Update maximum temperature
            if is_origin(phi[np.argmax(T)]):
                if is_origin(phi[np.argmin(T)]):
                    T_max = max(T_max, min(T) + 100)
                    continue
                for T_, phi_ in zip(p.T, p.phi):
                    if is_origin(phi_):
                        T_max = max(T_max, T_ + 100)
                        break

        for t in transitions:
            if t.key >0 : continue
            false_vacuum = t.false_vacuum[i]
            true_vacuum = t.true_vacuum[i]
            dx = true_vacuum - false_vacuum
            dy = 0.
            ax.quiver(false_vacuum, t.TC, dx, dy, **QUIVER_ARROW)

        ax.set_ylim(0, T_max)
        phi_shift = max((phi_max - phi_min)*0.05,5)
        ax.set_xlim(phi_min-phi_shift, phi_max+phi_shift)

        if len(axs) == 1:
            ax.set_xlabel(r"Field $x$ (GeV)")
        else:
            ax.set_xlabel(r"Field $x_{0}$ (GeV)".format(i+1))
        if i == 0:
            ax.set_ylabel("Temperature, $T$ (GeV)")

    if transitions:
        add_arrow_to_leg(axs[-1], "  FOPT")

    plt.savefig(pdf_name, bbox_inches="tight")


class Phase(object):
    def __init__(self, phase):
        self.T = phase[:, 0].T
        self.V = phase[:, 1].T
        self.phi = phase[:, 2:].T
        self.n_field = phase.shape[1] - 2

class Transition(object):
    def __init__(self, transition):
        self.n_field = (len(transition) - 2) // 2
        self.TC = transition[0]
        self.true_vacuum = transition[1:self.n_field+1]
        self.false_vacuum = transition[self.n_field+1:-1]
        self.key = transition[-1]

def load_data(dat_name):

    with open(dat_name) as d:
        data = d.read()

    data = re.sub(r'(?m)^#.*\n?', '', data)
    parts = data.split("\n\n")
    arrays = [np.genfromtxt(StringIO(p), delimiter=" ",
                            dtype=float) for p in parts]
    transitions = [Transition(a) for a in arrays if len(a.shape) != 2]
    phases = [Phase(a) for a in arrays if len(a.shape) == 2]
    return phases, transitions


if __name__ == "__main__":

    try:
        latex = os.environ["MATPLOTLIB_LATEX"]
    except KeyError:
        latex = False
    if latex:
        rc('text', usetex=True)
        rc('font', **{'family': 'serif', 'size': 14})
        GEV = r"\ensuremath{\,\textrm{GeV}}"
    else:
        GEV = r" GeV"


    rc('axes', **{'grid': True})
    rc('grid', **{'ls': ':'})
    rc('legend', **{'fontsize': FONTSIZE})

    # Load data
    prefix = sys.argv[1]
    dat_name = "{}.dat".format(prefix)
    phases, transitions = load_data(dat_name)

    # Plot phases on (phi_i, phi_j) plane
    pdf_name = "{}.pdf".format(prefix)
    plane_phi_phi(phases, transitions, pdf_name)

    # Plot phases on (phi_i, T) plane
    pdf_name = "phi_T_{}.pdf".format(prefix)
    plane_phi_T(phases, transitions, pdf_name)

    # Plot phases on (T, V) plane
    pdf_name = "V_T_{}.pdf".format(prefix)
    plane_T_V(phases, pdf_name)
