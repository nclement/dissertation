from scipy import signal
import numpy as np
import math
import matplotlib.pyplot as plt

def deg2rad(x):
    return x * math.pi / 180

def bivariate_vonmises_deg(phi, psi, mu, nu, sigma):
    math.exp(sigma * math.cos(deg2rad(phi-mu)) + sigma * math.cos(deg2rad(psi - nu)))

# Construct the data matrix.
fl = open('ARG_3_3_3_r1_pdf.txt', 'r')
# Skip the first 3 lines
fl.readline()
fl.readline()
fl.readline()

mat = np.zeros(shape=(360,360))
for l in fl.readlines():
    pieces = l.split()
    i = int(pieces[0]) + 180
    j = int(pieces[1]) + 180
    mat[i][j] = float(pieces[2])

imgplot = plt.imshow(np.real(mat))
print('okay, next thing')

# Construct the vm
vm_fn = lambda phi, psi: bivariate_vonmises_deg(phi, psi, 100, 0, 5)

vm_mat = np.fromfunction(np.vectorize(vm_fn), (360, 360))
plt.imshow(vm_mat)
