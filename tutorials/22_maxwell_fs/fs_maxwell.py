from dolfin import *
from matplotlib import interactive
from rbnics import *

import numpy as np
import matplotlib.pyplot as plt
import logging

logging.basicConfig(filename='adr.log', level=logging.DEBUG)
logging.info('Started')


class BoundaryM(UserExpression):
    def __init__(self, normal, **kwargs):
        super().__init__(kwargs)
        self.normal = normal

    def eval_cell(self, values, x, ufc_cell):
        if np.linalg.norm(np.array(x)) + DOLFIN_EPS > 5:
            t = np.array([[0, -self.normal(x)[2],
                           self.normal(x)[1]],
                          [self.normal(x)[2], 0, -self.normal(x)[0]],
                          [-self.normal(x)[1],
                           self.normal(x)[0], 0]])
            a = np.zeros((3, 3))
            values = list(np.block([[a, -t], [t.T, a]]).reshape(-1))
        else:
            values = [0] * 9

    def value_shape(self):
        return ((6, 6))


class Maxwell(FriedrichsSystemProblem):
    def __init__(self, L, **kwargs):
        FriedrichsSystemProblem.__init__(self, L, **kwargs)
        # ... and also store FEniCS data structures for assembly
        assert "subdomains" in kwargs
        assert "boundaries" in kwargs
        self.subdomains, self.boundaries, self.mesh = kwargs[
            "subdomains"], kwargs["boundaries"], kwargs["mesh"]

        # Defining the function spaces
        X = TensorFunctionSpace(self.mesh, "DG", 0, shape=(6, 6))

        self.y = TestFunction(L)
        self.z = TrialFunction(L)

        self.dx = Measure("dx")(domain=self.mesh,
                                subdomain_data=self.subdomains)
        self.dS = Measure("dS")(domain=self.mesh,
                                subdomain_data=self.subdomains)
        self.ds = Measure("ds")(domain=self.mesh,
                                subdomain_data=self.boundaries)

        self.n = FacetNormal(self.mesh)
        boundaryMesh = BoundaryMesh(self.mesh, 'exterior')
        self.normalCustom = self.get_facet_normal(boundaryMesh)
        self.alpha = Constant(10.)

        self.A_0_mu = Expression(
            (('1', '0', '0', '0', '0', '0'), ('0', '1', '0', '0', '0', '0'),
             ('0', '0', '1', '0', '0', '0'), ('0', '0', '0', '0', '0', '0'),
             ('0', '0', '0', '0', '0', '0'), ('0', '0', '0', '0', '0', '0')),
            element=X.ufl_element())
        self.A_0_epsilon = Expression(
            (('0', '0', '0', '0', '0', '0'), ('0', '0', '0', '0', '0', '0'),
             ('0', '0', '0', '0', '0', '0'), ('0', '0', '0', '1', '0', '0'),
             ('0', '0', '0', '0', '1', '0'), ('0', '0', '0', '0', '0', '1')),
            element=X.ufl_element())

        self.A_1 = Expression(
            (('0', '0', '0', '0', '0', '0'), ('0', '0', '0', '0', '0', '-1'),
             ('0', '0', '0', '0', '1', '0'), ('0', '0', '0', '0', '0', '0'),
             ('0', '0', '1', '0', '0', '0'), ('0', '-1', '0', '0', '0', '0')),
            element=X.ufl_element())
        self.A_2 = Expression(
            (('0', '0', '0', '0', '0', '1'), ('0', '0', '0', '0', '0', '0'),
             ('0', '0', '0', '-1', '0', '0'), ('0', '0', '-1', '0', '0', '0'),
             ('0', '0', '0', '0', '0', '0'), ('1', '0', '0', '0', '0', '0')),
            element=X.ufl_element())
        self.A_3 = Expression(
            (('0', '0', '0', '0', '-1', '0'), ('0', '0', '0', '1', '0', '0'),
             ('0', '0', '0', '0', '0', '0'), ('0', '1', '0', '0', '0', '0'),
             ('-1', '0', '0', '0', '0', '0'), ('0', '0', '0', '0', '0', '0')),
            element=X.ufl_element())

        self.D = self.A_1 * self.n[0] + self.A_2 * self.n[
            1] + self.A_3 * self.n[2]

        self.M = BoundaryM(normal=self.normalCustom, element=X.ufl_element())

        self.f = Expression(("0.0", "0.0", "0", "1", "1", "1"),
                            element=L.ufl_element())

    def get_facet_normal(self, bmesh):
        '''Manually calculate FacetNormal function'''

        vertices = bmesh.coordinates()
        cells = bmesh.cells()

        vec1 = vertices[cells[:, 1]] - vertices[cells[:, 0]]
        vec2 = vertices[cells[:, 2]] - vertices[cells[:, 0]]

        normals = np.cross(vec1, vec2)
        normals /= np.sqrt((normals**2).sum(axis=1))[:, np.newaxis]
        mask = np.linalg.norm(vertices[cells[:, 0]], axis=1) < 5. - DOLFIN_EPS
        normals[mask, :] = 0

        # plot
        origins = (vertices[cells[:, 1]] + vertices[cells[:, 0]] +
                   vertices[cells[:, 2]]) / 3
        origins[mask, :] = 0

        # Ensure outward pointing normal
        bmesh.init_cell_orientations(
            Expression(('x[0]', 'x[1]', 'x[2]'), degree=1))
        normals[bmesh.cell_orientations() == 1] *= -1
        # ax = plt.figure().add_subplot(projection='3d')
        # ax.quiver(origins[:, 0],
        #           origins[:, 1],
        #           origins[:, 2],
        #           normals[:, 0],
        #           normals[:, 1],
        #           normals[:, 2],
        #           length=0.5)
        # plt.show()

        V = VectorFunctionSpace(bmesh, 'DG', 0)
        norm = Function(V)
        nv = norm.vector()

        for n in (0, 1, 2):
            dofmap = V.sub(n).dofmap()
            for i in np.arange(dofmap.global_dimension()):
                dof_indices = dofmap.cell_dofs(i)
                assert len(dof_indices) == 1
                nv[dof_indices[0]] = normals[i, n]

        return norm

    def compute_theta(self, term):
        mu = self.mu
        if term == "a":
            theta_a0mu_void = mu[0]
            theta_a0epsilon_void = mu[3]
            theta_a0mu_copper = mu[1]
            theta_a0epsilon_copper = mu[4]
            theta_a0mu_iron = mu[2]
            theta_a0epsilon_iron = mu[5]
            theta_a123 = 1.
            theta_S = 1.
            theta_D = 1.
            theta_DpM = 1.
            return (theta_a0mu_void, theta_a0epsilon_void, theta_a0mu_copper,
                    theta_a0epsilon_copper, theta_a0mu_iron,
                    theta_a0epsilon_iron, theta_a123, theta_S, theta_D, theta_DpM)
        elif term == "f":
            theta_f0 = 1.
            return (theta_f0, )
        else:
            raise ValueError("Invalid term for compute_theta().")

    def assemble_operator(self, term):
        dx = self.dx
        z = self.z
        y = self.y
        if term == "a":
            a_a0mu_void = inner(z, dot(self.A_0_mu, y)) * dx(0)
            a_a0epsilon_void = inner(z, dot(self.A_0_epsilon, y)) * dx(0)
            a_a0mu_copper = inner(z, dot(self.A_0_mu, y)) * dx(2)
            a_a0epsilon_copper = inner(z, dot(self.A_0_epsilon, y)) * dx(2)
            a_a0mu_iron = inner(z, dot(self.A_0_mu, y)) * dx(1)
            a_a0epsilon_iron = inner(z, dot(self.A_0_epsilon, y)) * dx(1)
            a_a123 = -inner(z, dot(self.A_1, y.dx(0))) * dx - inner(
                z, dot(self.A_2, y.dx(1))) * dx - inner(
                    z, dot(self.A_3, y.dx(2))) * dx
            aS = self.alpha('+') * inner(jump(z), jump(y)) * dS
            a_D = inner(dot(self.D('+'), avg(z)), jump(y)) * dS
            a_DpM = 0.5 * inner(dot(self.D + self.M, z), y) * ds

            return (a_a0mu_void, a_a0epsilon_void, a_a0mu_copper,
                    a_a0epsilon_copper, a_a0mu_iron, a_a0epsilon_iron, a_a123, aS,
                    a_D, a_DpM)
        elif term == "f":
            f = self.f
            f0 = dot(y, f) * dx(2) - dot(y, f) * dx(1)
            return (f0, )
        elif term == "inner_product":
            x0 = inner(grad(z), grad(y)) * dx
            return (x0, )
        else:
            raise ValueError("Invalid term for assemble_operator().")


# 1. Read the mesh for this problem
meshSphere = Mesh("data/spheres.xml")
subdomains = MeshFunction("size_t", meshSphere,
                          "data/sphere_physical_region.xml")
boundaries = MeshFunction("size_t", meshSphere, "data/sphere_facet_region.xml")

# # 2. Create Finite Element space (Lagrange P1, two components)
L = VectorFunctionSpace(meshSphere, "DG", 1, dim=6)

# 3. Allocate an object of the Friedrichs' systems class
problem = Maxwell(L,
                  mesh=meshSphere,
                  subdomains=subdomains,
                  boundaries=boundaries)
mu_range = [(4*pi*1e-7, 4*pi*1e-7), (1.26e-6, 1.26e-6), (1e-6, 1e-5), (1, 1), (250000, 250000), (250000, 250000)]
problem.set_mu_range(mu_range)

# 4. Prepare reduction with a POD-Galerkin method
reduction_method = PODGalerkin(problem)
reduction_method.set_Nmax(20)

# 5. Perform the offline phase
reduction_method.initialize_training_set(1, sampling=UniformDistribution())
reduced_problem = reduction_method.offline()

# 6. Perform an online solve
online_mu = (4*pi*1e-7, 1.26e-6, 1e-5, 1, 250000, 250000)
print("online", len(mu_range), len(online_mu))
reduced_problem.set_mu(online_mu)
reduction_method.truth_problem.set_mu(online_mu)

# Plot solution
s = reduction_method.truth_problem.solve()
u_1 = s.sub(0)
u_2 = s.sub(1)
u_3 = s.sub(2)
u_4 = s.sub(3)
u_5 = s.sub(4)
u_6 = s.sub(5)


file = File("fs_maxwell_1.pvd")
file << u_1
file = File("fs_maxwell_2.pvd")
file << u_2
file = File("fs_maxwell_3.pvd")
file << u_3
file = File("fs_maxwell_4.pvd")
file << u_4
file = File("fs_maxwell_5.pvd")
file << u_5
file = File("fs_maxwell_6.pvd")
file << u_6
print("Saved online snapshot")

# reduced_problem.solve()
# reduced_problem.export_solution(filename="online_solution")

# # 7. Perform an error analysis
# reduction_method.initialize_testing_set(100, sampling=UniformDistribution())

# reduction_method.error_analysis()

# # 8. Perform a speedup analysis
# reduction_method.speedup_analysis()