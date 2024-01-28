import numpy as np
import matplotlib.pyplot as plt

Lx = 8
Ly = 8
D = 8
Db = 15

auto_correlation_data_len = 20
bond_num = Lx*(Ly-1) + (Lx-1)*Ly
site_num = Ly*Lx 

file_id = open(f'../data/square_energy_statistics{Ly}x{Lx}D{D}-{Db}', 'rb')
energy = np.fromfile(file_id, dtype=np.float64, count=1)[0]
en_std = np.fromfile(file_id, dtype=np.float64, count=1)[0]
bond_energys = np.fromfile(file_id, dtype=np.float64, count=bond_num)
energy_auto_corr = np.fromfile(file_id, dtype=np.float64, count=auto_correlation_data_len)


print(f"Energy ± en_std: {energy} ± {en_std}")

## Print bond_energys
#print("Bond Energies:")
#for i, bond_energy in enumerate(bond_energys):
#    print(f"Bond {i+1}: {bond_energy}")


# Plot energy_auto_corr if it has values
if energy_auto_corr.size > 0:
    plt.semilogy(np.arange(auto_correlation_data_len), energy_auto_corr)
    plt.xlabel('Time')
    plt.ylabel('Energy Auto-correlation')
    plt.title('Energy Auto-correlation Plot')
    plt.show()
else:
    print("No data found for energy_auto_corr")
