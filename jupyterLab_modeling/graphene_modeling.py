import pybinding as pb
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, pi

a = 0.24595   # [nm] unit cell length
a_cc = 0.142  # [nm] carbon-carbon distance
t = -2.8      # [eV] nearest neighbour hopping



def monolayer_graphene(a, t):

    lat = pb.Lattice(a1 = [3*a/2, sqrt(3)*a/2], a2=[3*a/2, -sqrt(3)*a/2])
    lat.add_sublattices(('a', [0, 0]), ('b', [a/2, sqrt(3)*a/2]))

    lat.add_hoppings(([0,  0], 'a', 'b', t),
                     ([-1, 1], 'a', 'b', t),
                     ([-1, 0], 'a', 'b', t))                   

    return lat


lat = monolayer_graphene(a_cc, t)

#print(type(lat))
#lat.plot()
#plt.show()

model = pb.Model(lat, 
        pb.translational_symmetry()
        )
hamiltonian = model.hamiltonian.todense()



Gamma = [0, 0]
K1 = [2*pi/(3*a_cc), 2*pi/(3*sqrt(3)*a_cc)]
M = [2*pi/(3*a_cc), 0]
K2 = [2*pi/(3*a_cc), -2*pi/(3*sqrt(3)*a_cc)]

solver = pb.solver.lapack(model)

bands = solver.calc_bands(K1, Gamma, M, K2)
bands.plot(point_labels=['K', r'$\Gamma$', 'M', 'K'])
plt.show()