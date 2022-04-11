from dolfin import *
from matplotlib import interactive
from rbnics import *

import numpy as np
import matplotlib.pyplot as plt


class Elastic(FriedrichsSystemProblem):
    def __init__(self, L, **kwargs):
        FriedrichsSystemProblem.__init__(self, L, **kwargs)
        # ... and also store FEniCS data structures for assembly
        assert "subdomains" in kwargs
        assert "boundaries" in kwargs
        self.subdomains, self.boundaries, self.mesh = kwargs[
            "subdomains"], kwargs["boundaries"], kwargs["mesh"]

        # Defining the function spaces
        dsu = TrialFunction(V)
        (self.s, self.u) = split(dsu)
        vq = TestFunction(V)
        (self.t, self.v) = split(vq)

        self.dx = Measure("dx")(domain=self.mesh,
                                subdomain_data=self.subdomains)
        self.dS = Measure("dS")(domain=self.mesh,
                                subdomain_data=self.subdomains)
        self.ds = Measure("ds")(domain=self.mesh,
                                subdomain_data=self.boundaries)

        self.n = FacetNormal(self.mesh)
        self.alpha = Constant(5)
        self.beta = Constant(0.1)
        self.h = CellDiameter(mesh)
        self.h_avg = (self.h('+') + self.h('-')) / 2

        self.f = Expression(
            ("-0.1*exp(-pow((x[0]-0.5), 2)-pow((x[1]-0.5), 2))",
             "-0.1*exp(-pow((x[0]-0.5), 2)-pow((x[1]-0.5), 2))"),
            element=V.sub(1).ufl_element())
        # self.f = Expression(("0", "0", "0", "0", "1", "1"),
        #                     element=L.ufl_element())

        self.I = Constant((('1', '0'), ('0', '1')))

        self.A_0_lambda = Constant(
            (('-1', '0', '0', '-1'), ('0', '0', '0', '0'),
             ('0', '0', '0', '0'), ('-1', '0', '0', '-1')))

    def compute_theta(self, term):
        mu = self.mu
        if term == "a":
            theta_a0 = 1.
            theta_a0_mu = mu[2]
            theta_a0_lambda = mu[0]
            theta_a123 = 1.
            theta_D = 1.
            theta_a123_lame2 = mu[1]
            theta_D_lame2 = mu[1]
            theta_DpM = mu[1]
            theta_S_i = 1.
            theta_S_f = 1.
            return (theta_a0, theta_a0_mu, theta_a0_lambda, theta_a123,
                    theta_D, theta_a123_lame2, theta_D_lame2, theta_DpM,
                    theta_S_i, theta_S_f)
        elif term == "f":
            theta_f0 = 1.
            return (theta_f0, )
        else:
            raise ValueError("Invalid term for compute_theta().")

    def assemble_operator(self, term):
        s = self.s
        u = self.u
        t = self.t
        v = self.v
        if term == "a":
            a_a0 = inner(s, t) * self.dx
            a_a0_mu = inner(u, v) * self.dx
            a_a0_lame1 = -tr(s) * tr(t) * self.dx

            a_a123 = 0.5 * inner(s + s.T, sym(grad(v))) * self.dx
            a_a123_lame2 = inner(u, 0.5 * div(t + t.T)) * self.dx

            # a_D = -inner(avg(u), 0.5*(dot(jump(t+t.T), self.n('+')))) * self.dS -inner(jump(v), 0.5*(dot(avg(s+s.T), self.n('+')))) * self.dS

            a_D_lame2 = -inner(avg(u), 0.5 *
                               (dot(jump(t + t.T), self.n('+')))) * self.dS
            a_D = -inner(jump(v), 0.5 *
                         (dot(avg(s + s.T), self.n('+')))) * self.dS

            a_DpM = -inner(u, 0.5 * dot(t + t.T, self.n)) * self.ds
            # a_DpM = -inner(v, 0.5 * dot(s + s.T, self.n)) * self.ds
            # a_DpM = 0.5 * inner(dot(self.DpM, z), y) * self.ds

            # stability terms (upwind): interfaces aS_i and boundaries aS_f
            aS_i = (-(self.alpha('+') / self.h_avg) *
                inner(0.5 * (dot(jump(t + t.T), self.n('+'))), 0.5 *
                      (dot(jump(s + s.T), self.n('+')))) +
                -(self.beta('+') / self.h_avg) * inner(jump(v), jump(u))) * self.dS
            aS_f = -(self.beta / self.h) * inner(
                0.5 * (dot((t + t.T), self.n)), 0.5 * (dot(
                    (s + s.T), self.n))) * self.ds

            return (a_a0, a_a0_mu, a_a0_lame1, a_a123, a_a123_lame2, a_D,
                    a_D_lame2, a_DpM, aS_i, aS_f)
        elif term == "f":
            f = self.f
            f0 = dot(v, f) * self.dx
            return (f0, )
        elif term == "inner_product":
            x0 = inner(grad(u), grad(v)) * self.dx
            return (x0, )
        else:
            raise ValueError("Invalid term for assemble_operator().")


# 1. Read the mesh for this problem
mesh = Mesh("data/elastic_block.xml")
subdomains = MeshFunction("size_t", mesh,
                          "data/elastic_block_physical_region.xml")
boundaries = MeshFunction("size_t", mesh,
                          "data/elastic_block_facet_region.xml")
# plot(mesh)
# plt.show()

# # 2. Create Finite Element space (Lagrange P1, two components)
# L_sigma = FunctionSpace(mesh, "RT", 1)
# Define function spaces and mixed (product) space
L_sigma = TensorElement("DG", mesh.ufl_cell(), 1, shape=(2, 2))
L_u = VectorElement("DG", mesh.ufl_cell(), 2, dim=2)
element = MixedElement([L_sigma, L_u])
V = FunctionSpace(mesh, element)

# 3. Allocate an object of the Friedrichs' systems class
problem = Elastic(V, mesh=mesh, subdomains=subdomains, boundaries=boundaries)
mu_range = [(1e-2, 1e-1), (1e-1, 1), (1e-4, 1e-2)]
problem.set_mu_range(mu_range)

# 4. Prepare reduction with a POD-Galerkin method
reduction_method = PODGalerkin(problem)
reduction_method.set_Nmax(20)

# 5. Perform the offline phase
reduction_method.initialize_training_set(1, sampling=LogUniformDistribution())
reduced_problem = reduction_method.offline()

# 6. Perform an online solve
online_mu = (0.1, 1, 1e-2)
# print("online", len(mu_range), len(online_mu))
reduced_problem.set_mu(online_mu)
reduction_method.truth_problem.set_mu(online_mu)

# Plot solution
s = reduction_method.truth_problem.solve()
q = abs(div(s.sub(1)))
p = plot(q)
plt.colorbar(p)
plt.title("div u mu={}".format(online_mu[0]))
plt.show()
p = plot(s.sub(1).sub(0))
plt.colorbar(p)
plt.title("x displacement mu={}".format(online_mu[0]))
plt.show()
p = plot(s.sub(1).sub(1))
plt.colorbar(p)
plt.title("y displacement")
plt.show()
p = plot(s.sub(0).sub(0))
plt.colorbar(p)
plt.title("stress 11 mu={}".format(online_mu[0]))
plt.show()
p = plot(s.sub(0).sub(3))
plt.colorbar(p)
plt.title("stress 12 mu={}")
plt.show()

# u_1 = s.sub(0)
# u_2 = s.sub(1)
# u_3 = s.sub(2)
# u_4 = s.sub(3)
# u_5 = s.sub(4)
# u_6 = s.sub(5)

# file = File("Elastic_1.pvd")
# file << u_1
# file = File("Elastic_2.pvd")
# file << u_2
# file = File("Elastic_3.pvd")
# file << u_3
# file = File("Elastic_4.pvd")
# file << u_4
# file = File("Elastic_5.pvd")
# file << u_5
# file = File("Elastic_6.pvd")
# file << u_6
print("Saved online snapshot")

# reduced_problem.solve()
# reduced_problem.export_solution(filename="online_solution")

# # 7. Perform an error analysis
# reduction_method.initialize_testing_set(100, sampling=UniformDistribution())

# reduction_method.error_analysis()

# # 8. Perform a speedup analysis
# reduction_method.speedup_analysis()