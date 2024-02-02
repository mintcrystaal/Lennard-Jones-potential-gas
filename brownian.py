import seaborn as sns, numpy as np
import matplotlib.pyplot as plt
import math
import scipy.stats as stats
from scipy.optimize import curve_fit
from matplotlib import rc



TMAX = 7000
NPART = 100


def func_line(x, k, b):
	return k * x + b

def func_curve2(x, a, b, c):
	return a * x * x + b * x + c

def func_curve3(x, a, b, c, d):
	return a * x * x * x + b * x * x + c * x + d

# def func_real_exp(x, D, c, tau):
# 	return 6 * D * x + c * np.exp(-x) + tau


def func_real_exp(x, D, c, tau):
 	return 6 * D * x + c * np.exp(-x / tau)


coords_all_times = []
for t in range(0, TMAX):
	filename1 = "dump_files/dump.gascheck." + str(t)
	filename2 = "coords_files_n_0.25_9/coords" + str(t) + ".txt"

	# file1 = open(filename1, "r")
	file2 = open(filename2, "r")

	coords_temp = []
	# for line in file1:
	# 	a = list(line.split())
	# 	if (len(a) == 11):
	# 		id1, type1, x, y, z, fx, fy, fz, vx, vy, vz = map(float, line.split())
	# 		# coords_temp.append([x, y, z])

	for line in file2:
		x, y, z = map(float, line.split())
		coords_temp.append([x, y, z])

	coords_all_times.append(coords_temp)

t1 = 0
T_PLOT_MIN = 4000
T_PLOT_MAX = 5000
DELT = 0.001
x_for_plot = []
y_for_plot = []

for t2 in range(T_PLOT_MIN, T_PLOT_MAX):
	sum_delt_r2 = 0.0
	for i in range(NPART):
		delt_r = [0.0, 0.0, 0.0]
		for d in range(3):
			delt_r[d] = coords_all_times[t1][i][d] - coords_all_times[t2][i][d]

		delt_r_abs2 = delt_r[0] ** 2 + delt_r[1] ** 2 + delt_r[2] ** 2
		# delt_r_abs = delt_r_abs2 ** 0.5
		sum_delt_r2 += delt_r_abs2

	delt_r_avg = sum_delt_r2 / NPART
	delt_t = abs(t2 - t1) * DELT

	x_for_plot.append(delt_t)
	y_for_plot.append(delt_r_avg)


kek = np.arange(T_PLOT_MIN * DELT, (T_PLOT_MAX + 1.0) * DELT, 0.1 * DELT)
# kek = np.arange(0.0, (T_PLOT_MAX + 1.0) * DELT, 0.1 * DELT)

# LINE

rc('font', **{'family': 'Times new roman'})
rc('text', usetex=True)
# rc('text.latex',unicode=True)
rc('text.latex', preamble=r'\usepackage[utf8]{inputenc}')
rc('text.latex', preamble=r'\usepackage[russian]{babel}')

popt, pcov = curve_fit(func_line, x_for_plot, y_for_plot)

perr = np.sqrt(np.diag(pcov))
print(popt, perr, "line")

label1 = str(round(popt[0], 3)) + " * x + " + str(round(popt[1], 3))

plt.plot(x_for_plot, y_for_plot, 'ro', markersize=5)
plt.plot(kek, func_line(kek, *popt), markersize=2, color="blue", label=label1)

D_reduced = popt[0] / 6
print(D_reduced, "D_reduced", T_PLOT_MIN, "-", T_PLOT_MAX)

plt.title("$\Delta r^2(t)$")
plt.ylabel("$\Delta r^2$")
plt.xlabel("t")
plt.grid(visible=True, which='major', axis='both', alpha=1)
plt.grid(visible=True, which='minor', axis='both', alpha=0.5)
plt.legend()
plt.savefig("EinsteinSmol_line_seed11.png")
plt.show()
plt.clf()


# # POL OF DEGREE 2

# popt, pcov = curve_fit(func_curve2, x_for_plot, y_for_plot)

# perr = np.sqrt(np.diag(pcov))
# print(popt, perr, "2")

# plt.plot(x_for_plot, y_for_plot, 'ro', markersize=2)

# label2 = str(round(popt[0], 3)) + " * x^2 + " + str(round(popt[1], 3)) + " * x + " + str(round(popt[2], 3))
# plt.plot(kek, func_curve2(kek, *popt), markersize=2, color="blue", label=label2)

# plt.title(r"delt_r2_avg(t)")
# plt.ylabel("delt_r2_avg")
# plt.xlabel(r"t")
# plt.grid(visible=True, which='major', axis='both', alpha=1)
# plt.grid(visible=True, which='minor', axis='both', alpha=0.5)
# plt.legend()
# plt.savefig("delt_r2_avg(t)2.png")
# plt.show()
# plt.clf()


# # 3 DEGREE POLYNMIAL
# plt.plot(x_for_plot, y_for_plot, 'ro', markersize=2)

# popt, pcov = curve_fit(func_curve3, x_for_plot, y_for_plot)
# perr = np.sqrt(np.diag(pcov))
# # popt[0] = -popt[0]
# label3 = str(round(popt[0], 3)) + " * x^3 + " + str(round(popt[1], 3)) + " * x^2 + " + str(round(popt[2], 3)) + " * x + " + str(round(popt[3], 3))
# plt.plot(kek, func_curve3(kek, *popt), markersize=2, color="blue", label=label3)
# print(popt, perr, "3")
# tau = (popt[1] / popt[0]) / 3
# print(tau, "tau")

# D = (popt[2] + popt[3] / tau) / 6

# plt.title(r"delt_r2_avg(t)")
# plt.ylabel("delt_r2_avg")
# plt.xlabel(r"t")
# plt.grid(visible=True, which='major', axis='both', alpha=1)
# plt.grid(visible=True, which='minor', axis='both', alpha=0.5)
# plt.legend()
# plt.savefig("delt_r2_avg(t)3.png")
# plt.show()
# plt.clf()

# REAL EXP

# plt.plot(x_for_plot, y_for_plot, 'ro', markersize=2)

# popt, pcov = curve_fit(func_real_exp, x_for_plot, y_for_plot)
# perr = np.sqrt(np.diag(pcov))
# # label3 = "6 * " + str(round(popt[0], 3)) + " * t + " + str(round(popt[1], 3)) + " * exp(-t/1)" + " + " + str(round(popt[2], 3))
# label3 = "6 * " + str(round(popt[0], 3)) + " * t + " + str(round(popt[1], 3)) + " * exp(-t/" + str(round(popt[2], 3)) + ")"
# plt.plot(kek, func_real_exp(kek, *popt), markersize=2, color="blue", label=label3)
# print(popt, perr, "exponent")

# plt.title(r"delt_r2_avg(t)")
# plt.ylabel("delt_r2_avg")
# plt.xlabel(r"t")
# plt.grid(visible=True, which='major', axis='both', alpha=1)
# plt.grid(visible=True, which='minor', axis='both', alpha=0.5)
# plt.legend()
# plt.savefig("delt_r2_avg(t)_exp.png")
# plt.show()
# plt.clf()



