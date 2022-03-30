from dolfin import *
from matplotlib import interactive
from rbnics import *

import numpy as np
import matplotlib.pyplot as plt

import time

# 1. Read the mesh for this problem
# mesh = Sphere(Point(0, 0, 0), 3)
# XDMFFile("cube.xdmf").write(mesh)
parameters["ghost_mode"] = "shared_facet"
mesh = Mesh("./spheres.xml")
subdomains = MeshFunction("size_t", mesh, "./spheres_physical_region.xml")
boundaries = MeshFunction("size_t", mesh, "./spheres_facet_region.xml")

# 2. Create Finite Element space (Lagrange P1, two components)
L = VectorElement("DG", mesh.ufl_cell(), 1, dim=3)

# Make a mixed space
TH = L * L
V = FunctionSpace(mesh, TH)
L_ = FunctionSpace(mesh, L)

dx = Measure("dx")(domain=mesh, subdomain_data=subdomains)
dS = Measure("dS")(domain=mesh, subdomain_data=subdomains)
ds = Measure("ds")(domain=mesh, subdomain_data=boundaries)

n = FacetNormal(mesh)
alpha = Constant(100.)

f = Expression(("pow(x[0],2)", "pow(x[1],2)", "pow(x[2],2)"), element=L_.ufl_element())

# Define trial and test functions
(H, E) = TrialFunctions(V)
(H_, E_) = TestFunctions(V)
mu = [1e-7, 1e-6, 1e-6]
epsilon = [1, 250000, 250000]

# Bilinear form
# a_interior = (inner(mu[0] * H, H_) + inner(epsilon[0] * E, E_))*dx + (inner(E,
# curl(H_)) - inner(H, curl(E_))) * dx

#cables
# a_interior =(inner(mu[0] * H, H_) + inner(epsilon[0] * E, E_)) * dx(21) + (inner(mu[1] * H, H_) + inner(epsilon[1] * E, E_)) * dx(22) + (inner(mu[2] * H, H_) + inner(epsilon[2] * E, E_)) * dx(23) + (inner(E, curl(H_)) - inner(H, curl(E_))) * dx

#incapsulated spheres
a_interior =(inner(mu[0] * H, H_) + inner(epsilon[0] * E, E_)) * dx(3) + (inner(mu[1] * H, H_) + inner(epsilon[1] * E, E_)) * dx(4) + (inner(E, curl(H_)) - inner(H, curl(E_))) * dx

a_interfaces = inner(cross(n('+'), avg(E)), jump(H_)) * dS - inner(
    cross(n('+'), avg(H)), jump(E_)) * dS + alpha('+') * inner(
        cross(n('+'), jump(E)), jump(E_)) * dS + alpha('+') * inner(
            cross(n('+'), jump(H)), jump(H_)) * dS

a_boundary = inner(cross(n('+'), E('+')), jump(H_)) * dS - inner(
    cross(n('+'), H('+')), jump(E_)) * dS

a = a_interior + a_interfaces + a_boundary

# Linear form
source = inner(E_, f) * dx

# Assemble and apply boundary conditions
A = assemble(a)
b = assemble(source)

# Solve system
z = Function(V)
solve(A, z.vector(), b)

# save solution
H, E = split(z)
H_1 = z.sub(0).sub(0, deepcopy=True)
H_2 = z.sub(0).sub(1, deepcopy=True)
H_3 = z.sub(0).sub(2, deepcopy=True)
E_1 = z.sub(1).sub(0, deepcopy=True)
E_2 = z.sub(1).sub(1, deepcopy=True)
E_3 = z.sub(1).sub(2, deepcopy=True)

file = File("H_1.pvd")
file << H_1
file = File("H_2.pvd")
file << H_2
file = File("H_3.pvd")
file << H_3
file = File("E_1.pvd")
file << E_1
file = File("E_2.pvd")
file << E_2
file = File("E_3.pvd")
file << E_3