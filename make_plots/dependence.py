import matplotlib.pyplot as plt
import numpy as np

from h_bar_expansion import call_pt


if __name__ == "__main__":

    fig, ax = plt.subplots(1, 4, figsize=(16, 4))

    # xi-dependence

    tree_ewsb = True
    tree_level_tadpoles = True
    lambda_hs = 0.2
    mtop = 173.03
    xi = np.linspace(0., 100, 100)

    data = np.array([call_pt(lambda_hs, mtop, x, tree_level_tadpoles, tree_ewsb) for x in xi])
    
    # Plot it

    ax[0].plot(xi, data[:, 0], color="red")

    ax[1].plot(xi, data[:, 1], color="red")
    
    ax[2].plot(xi, data[:, 2], color="red", ls="--", label="Tree-level")
    ax[2].plot(xi, data[:, 3], color="red", ls="-", label="One-loop")
    ax[2].legend(fontsize="small")

    ax[3].plot(xi, data[:, 4], color="red", ls="--", label="Tree-level")
    ax[3].plot(xi, data[:, 5], color="red", ls="-", label="One-loop")
    ax[3].legend(fontsize="small")

    # Labels etc
    
    ax[0].set_ylabel("$T_C$ (GeV)")
    ax[1].set_ylabel("$\gamma$")
    ax[2].set_ylabel("$m_h$ (GeV)")
    ax[3].set_ylabel("$v$ (GeV)")

    for a in ax:
        a.set_xlabel(r"$\xi$")

    plt.tight_layout()
    name = "xi_tree.png" if tree_level_tadpoles else "xi_one_loop.png"
    plt.savefig(name)
