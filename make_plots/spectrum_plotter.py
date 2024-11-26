import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

sys_call = sys.argv[0]
prefix = sys.argv[1]
ii = int(sys.argv[2])
peak_freq = float(sys.argv[3])
repo = sys.argv[4] if len(sys.argv) > 4 else ""

plot_txt_file = repo + "GW_" + prefix + "_" + str(ii) + ".txt"
pdf_name = repo + prefix + "_GW_spectrum_" + str(ii) + ".pdf"
plot_title = "Gravitational Wave Spectrum " + str(ii)

data = np.loadtxt(plot_txt_file, delimiter=',')

freq = data[:, 0]
total = data[:, 1]
acoustic = data[:, 2]
turbulence = data[:, 3]
collisions = data[:, 4]

# Create a figure
fig, ax = plt.subplots(figsize=(8, 6))

SIZE_TICK = 6
SIZE_DEFAULT = 10
SIZE_LARGE = 14
plt.rc("font", size=SIZE_DEFAULT, family='serif')
plt.rc("axes", titlesize=SIZE_LARGE)
plt.rc("axes", labelsize=SIZE_LARGE)
plt.rc("xtick", labelsize=SIZE_TICK)
plt.rc("ytick", labelsize=SIZE_TICK)

# Plot the data
ax.plot(freq, acoustic, '--', color='orange', linewidth=1.5, label='Acoustic')
ax.plot(freq, turbulence, '--', color='green', linewidth=1.5, label='Turbulence')
ax.plot(freq, collisions, '--', color='red', linewidth=1.5, label='Collisions')
ax.plot(freq, total, '-', color='blue', linewidth=3.0, label='Total amplitude')

ax.legend()
# ax.set_title('Total ' + plot_title)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$\mathrm{Frequency} \ (Hz)$')
ax.set_ylabel(r'$\mathrm{Amplitude} \ (\Omega_{GW} h^2)$')

ax.grid(True, which='major', color='lightgrey', linestyle=':', linewidth=0.5)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.yaxis.set_ticks_position("left")
ax.xaxis.set_ticks_position("bottom")

locmin = LogLocator(base=10.0, subs=(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), numticks=12)
ax.xaxis.set_minor_locator(locmin)

plt.tight_layout()
plt.savefig(pdf_name)
