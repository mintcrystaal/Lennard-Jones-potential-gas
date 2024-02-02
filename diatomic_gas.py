import numpy as np
from random import randint, random, seed, choice # random() -- gives a random number from 0 to 1
import math


# const
DELT = 0.001
NSTEPS = 400
NREPEAT = 10
DELT_SHORT = DELT / NREPEAT
TMAX = DELT * NSTEPS

NMOL = 100
NPART = NMOL * 2
BOX = 10.0
DIM = 3
RBAD = 0.88

L_BOND = 1.0
K_BOND = 5000
RC = 3.0
LAMBDA__ = 0.3
RM = 1.7


# pair: 2i, 2i + 1


def force(): # compute all forces on particles
	for i in range(NPART):
		for d in range(DIM):
			f[d][i] = 0.0

	for i in range(NPART):	
		for j in range(i + 1, NPART):
			dist = [0.0, 0.0, 0.0]
			for d in range(DIM):
				dist[d] = (coords[d][i] - coords[d][j] + BOX / 2) % BOX - BOX / 2

			r2 = dist[0] ** 2 + dist[1] ** 2 + dist[2] ** 2
			r = r2 ** 0.5

			r2i = 1 / r2
			r6i = r2i ** 3
			ff = 48.0 * r2i * (r6i * r6i - 0.5 * r6i)

			for d in range(DIM):
				f[d][i] += ff * dist[d]
				f[d][j] -= ff * dist[d]


def force_bond():
	for i in range(NPART):
		for d in range(DIM):
			f[d][i] = 0.0

	for i in range(NMOL):
		pA = 2 * i
		pB = 2 * i + 1
		dist = [0.0, 0.0, 0.0]
		for d in range(DIM):
			dist[d] = (coords[d][pA] - coords[d][pB] + BOX / 2) % BOX - BOX / 2

		r2 = dist[0] ** 2 + dist[1] ** 2 + dist[2] ** 2
		r = r2 ** 0.5

		force = -K_BOND * (r - L_BOND)

		for d in range(DIM):
			f[d][pA] += ((dist[d]) / r2) * force # sin alpha * force
			f[d][pB] -= ((dist[d]) / r2) * force


def integrate(dt): # integrate Newton's equations of motion
	global temp

	sumv2 = 0.0
	for i in range(NPART):
		vi2 = 0.0
		for d in range(DIM):
			v[d][i] += f[d][i] * dt / 2
			coords[d][i] += v[d][i] * dt
			coords[d][i] %= BOX
			vi2 += v[d][i] ** 2

		sumv2 += vi2

	temp = sumv2 / (3 * NPART)


def init():
	global temp
	sumv2 = 0.0
	sumv = [0.0, 0.0, 0.0]

	for i in range(NMOL):
		flag = True
		vec = [0.0, 0.0, 0.0]
		center_mass = [0.0, 0.0, 0.0]
		while (flag): 
			for d in range(DIM):
				vec[d] = random()
				center_mass[d] = randint(0, BOX - 1) + random()

			norm = (vec[0] ** 2 + vec[1] ** 2 + vec[2] ** 2) ** 0.5

			for d in range(DIM):
				vec[d] = vec[d] * L_BOND / norm

			for d in range(DIM):
				coords[d][2 * i] = center_mass[d] + vec[d] / 2
				coords[d][2 * i + 1] = center_mass[d] - vec[d] / 2

				coords[d][2 * i] %= BOX ######################
				coords[d][2 * i + 1] %= BOX ##################

			flag = False
			for j in range(2 * i):
				check_coords = [2 * i, 2 * i + 1]

				for it_ in check_coords:
					dist = [0.0, 0.0, 0.0]

					for d in range(DIM):
						dist[d] = (coords[d][it_] - coords[d][j] + BOX / 2) % BOX - BOX / 2

					r2 = dist[0] ** 2 + dist[1] ** 2 + dist[2] ** 2
					r = r2 ** 0.5

					if (r < RBAD):
						flag = True
						break

		for d in range(DIM):
			v[d][2 * i] = (random() - 0.5)
			v[d][2 * i + 1] = v[d][2 * i] #(random() - 0.5)

		v2_0 = v[0][2 * i] ** 2 + v[1][2 * i] ** 2 + v[2][2 * i] ** 2
		v2_1 = v[0][2 * i + 1] ** 2 + v[1][2 * i + 1] ** 2 + v[2][2 * i + 1] ** 2

		for d in range(DIM):
			sumv[d] += v[d][2 * i] + v[d][2 * i + 1]
		sumv2 += v2_0 + v2_1

	for d in range(DIM):
		sumv[d] /= NPART

	temp = 2.0

	sumv2 /= NPART
	fs = (3.0 * temp / sumv2) ** 0.5

	for cur_mol in range(NMOL):
		for d in range(DIM):
			these_coords = [2 * cur_mol, 2 * cur_mol + 1]
			for i in these_coords:
				v[d][i] = (v[d][i] - sumv[d]) * fs # velocity center of mass to zero


def make_file(t):
	time_string = str(t)

	into_file = "ITEM: TIMESTEP\n" + time_string + "\n"
	into_file += "ITEM: NUMBER OF ATOMS\n" + str(NPART) + "\n"
	into_file += "ITEM: BOX BOUNDS pp pp pp\n"

	into_file += "0.0 " + str(BOX) + "\n"
	into_file += "0.0 " + str(BOX) + "\n"
	into_file += "0.0 " + str(BOX) + "\n"
	into_file += "ITEM: ATOMS id type x y z fx fy fz vx vy vz\n"

	for i in range(NPART):

		into_file += str(i + 1) + " 1 "
		into_file += str(coords[0][i]) + " " + str(coords[1][i]) + " " + str(coords[2][i]) + " "
		into_file += str(f[0][i]) + " " + str(f[1][i]) + " " + str(f[2][i]) + " "
		into_file += str(v[0][i]) + " " + str(v[1][i]) + " " + str(v[2][i]) + "\n"
		
	n_ = NPART / (BOX ** 3)	
	filename = "dump_files" + "/dump.gascheck." + time_string
	file = open(filename, "w")
	file.write(into_file)
	file.close()


def makes_file_for_graphs():
	file_graph = open("data_for_graphs_diatomic", "w")
	into_file_graph_data = ""
	for i in range(len(temp_list)):
		into_file_graph_data += str(temp_list[i]) + "\n"
	file_graph.write(into_file_graph_data)
	file_graph.close()


# main
temp = 2.0

seed(7)

a = np.empty(NPART)
a.fill(0.0)

coords = np.full((DIM, NPART), a)

v = np.full((DIM, NPART), a)
f = np.full((DIM, NPART), a)

coords_previous = np.full((DIM, NPART), a)

temp_list = np.array([], dtype = float)
sumv = [0.0, 0.0, 0.0]

init()
make_file(0)

for t in range(NSTEPS):
	force()
	integrate(DELT)

	print(temp) #####

	for small_t in range(NREPEAT):
		force_bond()
		integrate(DELT_SHORT)
		temp_list = np.append(temp_list, temp)
		make_file(t)

		print(temp) #####

	temp_list = np.append(temp_list, temp)
	make_file(t)

makes_file_for_graphs()
