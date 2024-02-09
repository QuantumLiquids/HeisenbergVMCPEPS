import numpy as np
import matplotlib.pyplot as plt

num_mc_chain = 100
all_data = []

for i in range(num_mc_chain):
    en_loc_data = np.loadtxt('../../data/energy_raw_data6x24J0D10/energy' + str(i))
    all_data.extend(en_loc_data)

all_data = np.array(all_data)
#all_data = all_data[all_data < 0]

fig, ax = plt.subplots()
h = ax.hist(all_data, bins=10000)
ax.set_ylim([0.1, max(h[0])])
ax.set_yscale('log')

plt.show()
