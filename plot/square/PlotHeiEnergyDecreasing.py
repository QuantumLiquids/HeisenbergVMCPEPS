import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Load energy data
data_nkd = np.loadtxt("../../data/energy_decrease8x8_nkd.dat", delimiter=",", dtype=float)
data_cent = np.loadtxt("../../data/energy_decrease8x8_cent.dat", delimiter=",", dtype=float)

# Calculate exact energy
exact_energy = -0.619040 * 64

# Calculate relative energy error
relative_error_nkd = np.abs((data_nkd - exact_energy) / exact_energy)
relative_error_cent = np.abs((data_cent - exact_energy) / exact_energy)

# Plot energy decrease figure
fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(data_nkd, label="Energy Data (NKD)")
ax.plot(data_cent, label="Energy Data (CENT)")
ax.axhline(y=exact_energy, color='r', linestyle='--', label="Exact Energy")
ax.set_xlabel("Iteration")
ax.set_ylabel("Energy")
ax.legend()
ax.grid(True)

# Create inset for relative energy error
ax_inset = inset_axes(ax, width="40%", height="30%", loc='upper right')
ax_inset.semilogy(relative_error_nkd, label="Relative Error (NKD)")
ax_inset.semilogy(relative_error_cent, label="Relative Error (CENT)")
ax_inset.set_xlabel("Iteration")
ax_inset.set_ylabel("Relative Energy Error")
ax_inset.legend()
ax_inset.grid(True)

# Show the plots
plt.show()
