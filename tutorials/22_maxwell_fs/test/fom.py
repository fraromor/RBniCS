from dolfin import *
from matplotlib import interactive
from rbnics import *

import numpy as np
import matplotlib.pyplot as plt

import time

from aux import  BoundaryM, get_facet_normal

# 1. Read the mesh for this problem
mesh = Mesh("./cables.xml")
subdomains = MeshFunction("size_t", mesh,
                          "./cables_physical_region.xml")
boundaries = MeshFunction("size_t", mesh, "./cables_facet_region.xml")
# XDMFFile("cables.xdmf").write(mesh)
# XDMFFile("cables_physical_region.xdmf").write(subdomains)
# XDMFFile("cables_facet_region.xdmf").write(boundaries)

# 2. Create Finite Element space (Lagrange P1, two components)
L = VectorFunctionSpace(mesh, "DG", 1, dim=6)
X = TensorFunctionSpace(mesh, "DG", 0, shape=(6, 6))

dx = Measure("dx")(domain=mesh,
                                subdomain_data=subdomains)
dS = Measure("dS")(domain=mesh,
                        subdomain_data=boundaries)
ds = Measure("ds")(domain=mesh,
                        subdomain_data=boundaries)

n = FacetNormal(mesh)
boundaryMesh = BoundaryMesh(mesh, 'exterior')
normalCustom = get_facet_normal(boundaryMesh)
alpha = Constant(1.)

A_0_mu = Expression(
    (('1', '0', '0', '0', '0', '0'),
    ('0', '1', '0', '0', '0', '0'),
    ('0', '0', '1', '0', '0', '0'),
    ('0', '0', '0', '0', '0', '0'),
    ('0', '0', '0', '0', '0', '0'),
    ('0', '0', '0', '0', '0', '0')),
    element=X.ufl_element())
A_0_epsilon = Expression(
    (('0', '0', '0', '0', '0', '0'),
    ('0', '0', '0', '0', '0', '0'),
    ('0', '0', '0', '0', '0', '0'),
    ('0', '0', '0', '1', '0', '0'),
    ('0', '0', '0', '0', '1', '0'),
    ('0', '0', '0', '0', '0', '1')),
    element=X.ufl_element())

A_1 = Expression(
    (('0', '0', '0', '0', '0', '0'),
        ('0', '0', '0', '0', '0', '-1'),
        ('0', '0', '0', '0', '1', '0'),
        ('0', '0', '0', '0', '0', '0'),
        ('0', '0', '1', '0', '0', '0'),
        ('0', '-1', '0', '0', '0', '0')),
    element=X.ufl_element())
A_2 = Expression(
    (('0', '0', '0', '0', '0', '1'),
        ('0', '0', '0', '0', '0', '0'),
        ('0', '0', '0', '-1', '0', '0'),
        ('0', '0', '-1', '0', '0', '0'),
        ('0', '0', '0', '0', '0', '0'),
        ('1', '0', '0', '0', '0', '0')),
    element=X.ufl_element())
A_3 = Expression(
    (('0', '0', '0', '0', '-1', '0'),
        ('0', '0', '0', '1', '0', '0'),
        ('0', '0', '0', '0', '0', '0'),
        ('0', '1', '0', '0', '0', '0'),
        ('-1', '0', '0', '0', '0', '0'),
        ('0', '0', '0', '0', '0', '0')),
    element=X.ufl_element())

D = A_1 * n[0] + A_2 * n[1] + A_3 * n[2]

M = BoundaryM(normal=normalCustom, element=X.ufl_element())

f = Expression(("0", "0", "0", "0", "0", "100"),
                    element=L.ufl_element())

y = TestFunction(L)
z = TrialFunction(L)

# Bilinear form
a_a0mu_void = inner(z, dot(A_0_mu, y)) * dx(21)
a_a0epsilon_void = inner(z, dot(A_0_epsilon, y)) * dx(21)
a_a0mu_copper = inner(z, dot(A_0_mu, y)) * dx(22)
a_a0epsilon_copper = inner(z, dot(A_0_epsilon, y)) * dx(22)
a_a0mu_iron = inner(z, dot(A_0_mu, y)) * dx(23)
a_a0epsilon_iron = inner(z, dot(A_0_epsilon, y)) * dx(23)

a_a123 = - inner(z, dot(A_1, y.dx(0))) * dx - inner(z, dot(A_2, y.dx(1))) * dx - inner(z, dot(A_3, y.dx(2))) * dx

aS = alpha('+') * inner(jump(z), jump(y)) * dS
a_D = inner(dot(D('+'), avg(z)), jump(y)) * dS
a_DpM = 0.5 * inner(dot(D + M, z), y) * ds

# a_interior =  a_a123#inner(z, dot(A_0_mu, y)) * dx + inner(z, dot(A_0_epsilon, y)) * dx + a_a123
a_interior = 0.0001*a_a0mu_void + a_a0epsilon_void + 0.2*a_a0mu_copper + a_a0epsilon_copper + 0.1*a_a0mu_iron + a_a0epsilon_iron + a_a123

a_interfaces =  a_D # + aS

a_boundaries = a_DpM

a = a_interior + a_interfaces + a_boundaries

# Linear form
source = inner(y, f) * dx

# Assemble and apply boundary conditions
A = assemble(a)
b = assemble(source)

# Solve system
z = Function(L)
solve(A, z.vector(), b)

# save solution
H_1 = z.sub(0)
H_2 = z.sub(1)
H_3 = z.sub(2)
E_1 = z.sub(3)
E_2 = z.sub(4)
E_3 = z.sub(5)

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