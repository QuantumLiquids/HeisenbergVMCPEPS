import numpy as np
import matplotlib.pyplot as plt

Lx = 8
Ly = 8
D = 8
Db = 15

file_id = open(f'../../data/square_two_point_functions{Ly}x{Lx}D{D}-{Db}', 'rb')
data = np.fromfile(file_id, dtype=np.float64, count=Lx * 3)

# Split the data into three parts: sz-sz correlation, s_plus-s_minus correlation, and s_minus-s_plus correlation
num_points = Lx // 2
sz_sz_corr = data[:num_points]
s_plus_s_minus_corr = data[num_points:2*num_points]
s_minus_s_plus_corr = data[2*num_points:3*num_points]
sz_sz_corr_err = data[3*num_points:4*num_points]
s_plus_s_minus_corr_err = data[4*num_points:5*num_points]
s_minus_s_plus_corr_err = data[5*num_points:6*num_points]

# Plot the loaded data
x = np.arange(1, num_points+1)  # x-axis values
plt.errorbar(x, sz_sz_corr * (-1) ** x, yerr=sz_sz_corr_err, label=r'$\langle S^z(0)S^z(r)\rangle$')
plt.errorbar(x, 0.5 * s_plus_s_minus_corr * (-1) ** x, yerr=s_plus_s_minus_corr_err, label=r'$0.5 * \langle S^{+}(0)S^-(r)\rangle$')
plt.errorbar(x, 0.5 * s_minus_s_plus_corr * (-1) ** x, yerr=s_minus_s_plus_corr_err, label=r'$0.5 * \langle S^-(0)S^+(r)\rangle$')
plt.xlabel('Distance')
plt.ylabel('Spin Correlation')
plt.legend()
plt.show()
