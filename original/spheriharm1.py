#!/usr/bin/env python3

"""
@author: daniel
"""

import sys

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.special as ss
from matplotlib import cm, colors

dmax = 6


# read your coefficients C and S from file spheriharm1.dat
C = np.zeros(shape=(dmax + 1, dmax + 1), dtype="float")
S = np.zeros(shape=(dmax + 1, dmax + 1), dtype="float")

nline = 0
with open("spheriharm1.dat", "r") as f:
    for line in f:
        data = line.split()
        if nline > 1:
            for nfield in range(1, len(data)):
                C[nline - 2, nfield - 1] = data[nfield].split("(")[1].split(",")[0]
                S[nline - 2, nfield - 1] = data[nfield].split(",")[1].split(")")[0]
        nline += 1

# set up the grid
vphi = np.linspace(0, 2 * np.pi, 90)  # longitude
vtheta = np.linspace(0, np.pi, 45)  # co-lattitude
# join the grid
theta, phi = np.meshgrid(vtheta, vphi)
sphrad = 10  # %radius of sphere

# initialize figure
norm = colors.Normalize()
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
m = cm.ScalarMappable(cmap=cm.jet)
RR2 = np.zeros(shape=(len(vphi), len(vtheta)), dtype="float")
for degree in range(int(dmax) + 1):
    # print ('Working on degree ',int(degree))
    for order in range(degree + 1):
        # compute Spherical harmonic:
        Plm = ss.lpmv(int(order), int(degree), np.cos(vtheta))
        SC = C[degree, order] * np.array(np.cos(vphi * float(order))) + S[
            degree, order
        ] * np.array(np.sin(vphi * float(order)))
        f = sp.special.factorial(int(degree) - int(order)) / sp.special.factorial(
            int(degree) + int(order)
        )
        K_norm = np.sqrt((2 * int(degree) + 1) / (4 * np.pi) * f)
        R = K_norm * np.array(np.outer(np.array(SC), np.array(Plm)))
        RR2 = np.add(RR2, R)

rho = 0.1 * np.array(sphrad + np.array(RR2) / np.amax(RR2))
x = rho * np.sin(theta) * np.cos(phi)
y = rho * np.sin(theta) * np.sin(phi)
z = rho * np.cos(theta)
ax.plot_surface(
    x, y, z, rstride=1, cstride=1, linewidth=0.5, facecolors=cm.jet(norm(RR2))
)
ax.set_title("Your Spherical Harmonic Summation", fontsize=16)
m.set_array(RR2)
cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])  # [left, bottom, width, height]
fig.colorbar(m, cax=cbar_ax)
ax.view_init(float(sys.argv[1]), float(sys.argv[2]))
# print(sys.argv[3])
plt.savefig(
    f"Imagenes/SphericalHarmonic_{sys.argv[3]}_{sys.argv[1]}_{sys.argv[2]}.png",
    dpi=300,
)

# plt.show()


# done
