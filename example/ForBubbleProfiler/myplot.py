import os
import random
import numpy as np
import matplotlib.pyplot as plt

E = 1.
alpha = 0.6
scale = 20

def v(phi):
  phi = phi/scale;
  return -pow(scale, 4) * E * (-alpha * phi * phi * phi * phi
                                 + phi * phi * phi
                                 + 0.5 * (4. * alpha - 3.) * phi * phi);


fig, axes = plt.subplots(1, 1, figsize=(6, 6))

x = np.linspace(-20, 50, 100)

cm=axes.plot(x, v(x))

axes.set_ylim(-2E4,1E4)

#
##axes.axhline(y=140, color='r', linestyle='--')
#
#
#axes.set_ylabel(r"$S/T$")
#axes.set_xlabel(r"$T$")
#
#cbar=plt.colorbar(cm)
#cbar.set_label(r'log$_{10}(v_f)$')
#
#
##axes.plot(data1[:,0],data1[:,1],'.')
##axes.set_ylabel(r"$v_f$")
##axes.set_xlabel(r"$v_t$")
#
plt.tight_layout()
plt.show()
