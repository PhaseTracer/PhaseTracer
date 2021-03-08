import matplotlib.pyplot as plt
import numpy as np

from call_pt import call_pt, transpose_list_dict
from style import style


if __name__ == "__main__":

    style()
    fig, ax_ = plt.subplots(2, 2, figsize=(8, 8))
    ax = ax_.flatten().tolist()

    # xi-dependence

    tree_ewsb = False
 
    lambda_hs = 0.2
    mtop = 173.03
    xi = np.linspace(0., 3, 100)

    for tree_level_tadpoles in [True, False]:

        data = transpose_list_dict([call_pt("h_bar_expansion", lambda_hs, mtop, x, tree_level_tadpoles, tree_ewsb) for x in xi])  
        
        # Plot it

        c = "red" if tree_level_tadpoles else "blue"
        ax[0].plot(xi, data["TC"], color=c )

        ax[1].plot(xi, data["gamma"], color=c)
        
        ax[2].plot(xi, data["mh_tree"], color=c, ls="--", label="Tree-level")
        ax[2].plot(xi, data["mh_1l"], color=c, ls="-", label="One-loop")
        ax[2].legend(fontsize="small")

        ax[3].plot(xi, data["v_tree"], color=c, ls="--", label="Tree-level")
        ax[3].plot(xi, data["v_1l"], color=c, ls="-", label="One-loop")
        ax[3].legend(fontsize="small")

    # Labels etc
    
    ax[0].set_ylabel("$T_C$ (GeV)")
    ax[1].set_ylabel("$\gamma$")
    ax[2].set_ylabel("$m_h$ (GeV)")
    ax[3].set_ylabel("$v$ (GeV)")

    for a in ax:
        a.set_xlabel(r"$\xi$")

    plt.tight_layout()
    plt.savefig("xi_detail_h_bar_expansion.pdf")
