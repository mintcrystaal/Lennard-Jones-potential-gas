import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from math import floor, ceil


DELT = 0.001


mpl.rcParams['font.size'] = 16 # Управление стилем, в данном случаем - размером шрифта 
 # Создаем фигуру
plt.figure(figsize=(7,7))

# Подписываем оси и график
plt.title(r"potential energy")
plt.ylabel("e_n(t)")
plt.xlabel(r"t")

file = open("data_for_graphs", "r")
x_coords = []
e_n_coords = []
e_tot_coords = []
temp_coords = []
delt_e_coords = []
v_tot_coords = []


t = 0.0
for line in file:
    e_n, e_tot, temp, delt_e, v_tot = map(float, line.split())
    x_coords.append(t)
    e_n_coords.append(e_n)
    e_tot_coords.append(e_tot)
    temp_coords.append(temp)
    delt_e_coords.append(delt_e)
    v_tot_coords.append(v_tot)

    t += DELT

file.close()

plt.plot(x_coords, e_n_coords, 'ro')

plt.grid(visible=True, which='major', axis='both', alpha=1)
plt.grid(visible=True, which='minor', axis='both', alpha=0.5)
plt.legend()
plt.savefig('e_n(t).png')
plt.show()
plt.clf()

plt.title(r"total energy")
plt.ylabel("e_tot(t)")
plt.xlabel(r"t")

plt.grid(visible=True, which='major', axis='both', alpha=1)
plt.grid(visible=True, which='minor', axis='both', alpha=0.5)
plt.plot(x_coords, e_tot_coords, 'ro')
plt.legend()
plt.savefig('e_tot(t).png')
plt.show()
plt.clf()

plt.title(r"temperature")
plt.ylabel("temp(t)")
plt.xlabel(r"t")

plt.grid(visible=True, which='major', axis='both', alpha=1)
plt.grid(visible=True, which='minor', axis='both', alpha=0.5)
plt.plot(x_coords, temp_coords, 'ro')
plt.legend()
plt.savefig('temp(t).png')
plt.show()

plt.title(r"dE/temp(t)")
plt.ylabel("dE/temp")
plt.xlabel(r"t")
plt.grid(visible=True, which='major', axis='both', alpha=1)
plt.grid(visible=True, which='minor', axis='both', alpha=0.5)
plt.plot(x_coords, delt_e_coords, 'ro')
plt.legend()
plt.savefig('dE(t).png')
plt.show()

plt.title(r"momentum(t)")
plt.ylabel("momentum")
plt.xlabel(r"t")
plt.grid(visible=True, which='major', axis='both', alpha=1)
plt.grid(visible=True, which='minor', axis='both', alpha=0.5)
plt.plot(x_coords, v_tot_coords, 'ro')
plt.legend()
plt.savefig('momentum(t).png')
plt.show()






# Сохраняем изображение в текущую директорию
