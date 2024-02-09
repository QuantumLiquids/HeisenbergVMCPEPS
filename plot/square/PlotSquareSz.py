import numpy as np
import matplotlib.pyplot as plt

Lx = 8
Ly = 8
D = 8
Db = 15

site_num = Ly * Lx
auto_correlation_data_len = 20

file_id = open(f'../../data/square_one_point_functions{Ly}x{Lx}D{D}-{Db}', 'rb')
sz = np.fromfile(file_id, dtype=np.float64, count=site_num)
sz_err = np.fromfile(file_id, dtype=np.float64, count=site_num)
sz_auto_corr = np.fromfile(file_id, dtype=np.float64, count=20)

# Print sz ± sz_err
print("sz ± sz_err:")
for i in range(site_num):
    print(f"Site {i+1}: {sz[i]} ± {sz_err[i]}")


# Plot Spin Auto correlation 
plt.plot(np.arange(auto_correlation_data_len), sz_auto_corr)
print(sz_auto_corr)
plt.xlabel('Time')
plt.ylabel('Sz Auto-correlation')
plt.show()

## Create a grid of coordinates for the lattice sites
#x = np.arange(0, Lx)
#y = np.arange(0, Ly)
#X, Y = np.meshgrid(x, y)
#
#scaled_sz = 400 * sz
## Plot the arrows
#plt.quiver(X, Y, np.zeros_like(scaled_sz),  scaled_sz, cmap='coolwarm', scale=50, pivot='mid')
#plt.xlabel('X')
#plt.ylabel('Y')
#plt.title('Spin Magnitude Arrows on Square Lattice')
#plt.colorbar(label='Spin Magnitude')
#plt.show()
