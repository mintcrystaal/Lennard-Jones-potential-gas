import numpy as np
from random import randint, random, seed, choice # random() -- gives a random number from 0 to 1


# const
DELT = 0.001 # 0.0001
NPART = 100
BOX = 10.0
DELX = 1
TMAX = DELT * 500 #0.0100
RBAD = 0.88
DIM = 3


def e_cut(r):
	r6 = r ** 6
	r6i = 1 / r6
	return 4 * r6i * (r6i - 1)


def force(): # compute all forces on particles
	global e_n
	e_n = 0.0
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

			# if (r < RBAD):
			# 	print("too close :( ", CNT, i, j, r)
			# 	# exit(0)

			for d in range(DIM):
				f[d][i] += ff * dist[d]
				f[d][j] -= ff * dist[d]
	
			e_n += 4.0 * r6i * (r6i - 1)


def integrate(): # integrate Newton's equations of motion
	global e_tot, e_n, temp, sumv
	sumv = [0.0, 0.0, 0.0]

	sumvi = 0.0
	sumv2 = 0.0
	for i in range(NPART):
		new_coords = [0.0, 0.0, 0.0]

		for d in range(DIM):
			new_coords[d] = 2 * coords[d][i] - coords_previous[d][i] + (DELT ** 2) * f[d][i]
			v[d][i] = (new_coords[d] - coords_previous[d][i]) / (2 * DELT)

			if (new_coords[d] > BOX or new_coords[d] < 0):
				if (new_coords[d] % BOX > new_coords[d]):
					coords_previous[d][i] = coords[d][i] + BOX
				else:
					coords_previous[d][i] = coords[d][i] - BOX

				new_coords[d] = new_coords[d] % BOX
				coords[d][i] = new_coords[d]
			else:
				coords_previous[d][i] = coords[d][i]
				coords[d][i] = new_coords[d]

		vi2 = v[0][i] ** 2 + v[1][i] ** 2 + v[2][i] ** 2
		vi = vi2 ** 0.5

		sumvi += vi
		sumv2 += vi2

		for d in range(DIM):
			sumv[d] += v[d][i]

	temp = sumv2 / (3 * NPART)
	e_tot = (e_n + 0.5 * sumv2) / NPART


def init():
	global temp
	sumv2 = 0.0
	sumv = [0.0, 0.0, 0.0]

	for i in range(NPART):
		flag = True
		while (flag): 
			for d in range(DIM):
				coords[d][i] = randint(0, BOX - 1) + random()

			flag = False
			for j in range(i):
				dist = [0.0, 0.0, 0.0]

				for d in range(DIM):
					dist[d] = (coords[d][i] - coords[d][j] + BOX / 2) % BOX - BOX / 2

				r2 = dist[0] ** 2 + dist[1] ** 2 + dist[2] ** 2
				r = r2 ** 0.5

				if (r < RBAD):
					flag = True
					break

		vi = 0.0
		for d in range(DIM):
			v[d][i] = (random() - 0.5)

		v2 = v[0][i] ** 2 + v[1][i] ** 2 + v[2][i] ** 2
		vi = v2 ** 0.5

		for d in range(DIM):
			sumv[d] += v[d][i]
			# sumv2c[d] += v[d][i] * v[d][i]
		sumv2 += v2

	for d in range(DIM):
		sumv[d] /= NPART

	temp = 2.0 #sumv2 / (3 * NPART)

	sumv2 /= NPART
	fs = (3.0 * temp / sumv2) ** 0.5

	for i in range(NPART):
		for d in range(DIM):
			v[d][i] = (v[d][i] - sumv[d]) * fs # velocity center of mass to zero

			coords_previous[d][i] = coords[d][i] - v[d][i] * DELT # position previous time step
			
			if (coords_previous[d][i] > BOX or coords_previous[d][i] < 0):
				if (coords_previous[d][i] % BOX > coords_previous[d][i]):
					coords_previous[d][i] -= BOX
				else:
					coords_previous[d][i] += BOX


def make_file():
	time_string = str(CNT)

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
		
	filename = "dump_files/dump.gascheck." + time_string
	file = open(filename, "w")
	file.write(into_file)
	file.close()


def makes_file_for_graphs():
	file_graph = open("data_for_graphs", "w")
	into_file_graph_data = ""
	for i in range(len(e_n_list)):
		into_file_graph_data += str(e_n_list[i]) + " " + str(e_tot_list[i]) + " " + str(temp_list[i]) + " " + str(delt_e_list[i]) + " " + str(v_tot_list[i]) + "\n"
	file_graph.write(into_file_graph_data)
	file_graph.close()


# main
t = 0.0 # time
temp = 2.0
CNT = 0

# seed(4)

a = np.empty(NPART)
a.fill(0.0)

coords = np.full((DIM, NPART), a)
v = np.full((DIM, NPART), a)
f = np.full((DIM, NPART), a)

coords_previous = np.full((DIM, NPART), a)

delt_e_list = np.array([], dtype = float)

e_n = 0.0 # potential energy
e_tot = 0.0 # total energy

e_n_list = np.array([], dtype = float)
e_tot_list = np.array([], dtype = float)
temp_list = np.array([], dtype = float)
v_tot_list = np.array([], dtype = float)
sumv = [0.0, 0.0, 0.0]

init()
make_file()
CNT += 1
t += DELT
# e_n_list.append(e_n)
# e_tot_list.append(e_tot)
# temp_list.append(temp)

# print(temp, "temp")

while t <= TMAX:
	force()
	integrate()
	# print(temp, "temp")

	e_n_list = np.append(e_n_list, e_n)
	e_tot_list = np.append(e_tot_list, e_tot)
	temp_list = np.append(temp_list, temp)
	v_tot = (sumv[0] ** 2 + sumv[1] ** 2 + sumv[2] ** 2)
	# v_tot = v_tot ** 0.5
	v_tot_list = np.append(v_tot_list, v_tot)
	delt_e_list = np.append(delt_e_list, 0.0) #np.append(delt_e_list, ((e_tot - e_tot_was) / temp))
	make_file()
	t += DELT
	CNT += 1
print(CNT)
makes_file_for_graphs()
