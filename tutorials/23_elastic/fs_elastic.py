from dolfin import *
from matplotlib import interactive
from rbnics import *

import numpy as np
import matplotlib.pyplot as plt

class BoundaryM(UserExpression):
    def __init__(self, normal, **kwargs):
        super().__init__(kwargs)
        self.normal = normal

    def eval_cell(self, values, x, ufc_cell):

        if (x[0] < DOLFIN_EPS
                or 1. - x[1] < DOLFIN_EPS) or (x[1] < DOLFIN_EPS
                                               or 1. - x[1] < DOLFIN_EPS):
            N = np.array([
                [-self.normal(x)[0], 0],
                [-0.5 * self.normal(x)[1], -0.5 * self.normal(x)[0]],
                [-0.5 * self.normal(x)[1], -0.5 * self.normal(x)[0]],
                [0, -self.normal(x)[1]],
            ])
            a = np.zeros((4, 4))
            b = np.zeros((2, 2))

            # Homogeneous Neumann
            values = np.block([[a, N], [-N.T, b]]).reshape(-1)
            # # Homogeneous Dirichlet
            # values = np.block([[a, -t], [t.T, b]]).reshape(-1)
        else:
            values = [0] * 9

    def value_shape(self):
        return ((6, 6))

class BoundaryDpM(UserExpression):
    def __init__(self, normal, **kwargs):
        super().__init__(kwargs)
        self.normal = normal

    def eval_cell(self, values, x, ufc_cell):

        if (x[0] < DOLFIN_EPS
                or 1. - x[1] < DOLFIN_EPS) or (x[1] < DOLFIN_EPS
                                               or 1. - x[1] < DOLFIN_EPS):
            N = np.array([
                [-self.normal(x)[0], 0],
                [-0.5 * self.normal(x)[1], -0.5 * self.normal(x)[0]],
                [-0.5 * self.normal(x)[1], -0.5 * self.normal(x)[0]],
                [0, -self.normal(x)[1]],
            ])
            a = np.zeros((4, 4))
            b = np.zeros((2, 2))
            c = np.zeros((4, 2))

            # Homogeneous Neumann
            values = np.block([[a, c], [2*(-N.T), b]]).reshape(-1)
            # # Homogeneous Dirichlet
            # values = np.block([[a, -t], [t.T, b]]).reshape(-1)
        else:
            values = [0] * 9

    def value_shape(self):
        return ((6, 6))

# Dirichlet/Neumann homogeneous
#! TO DEBUG
class InterfaceSi(UserExpression):
    def __init__(self, normal, mesh, **kwargs):
        super().__init__(kwargs)
        self.normal = normal
        self.mesh = mesh

    def eval_cell(self, values, x, ufc_cell):
        cell = Cell(self.mesh, ufc_cell.index)
        n = cell.normal(ufc_cell.local_facet)
        N = np.array([
                [-n[0], 0],
                [-0.5 * n[1], -0.5 * n[0]],
                [-0.5 * n[1], -0.5 * n[0]],
                [0, -n[1]],
            ])
        b = np.zeros((4, 2))
        a = np.eye(2)
        values = np.block([[N.dot(N.T), b], [b.T, a]]).reshape(-1)

    def value_shape(self):
        return ((6, 6))

# Neumann only
#! TO DEBUG
class InterfaceSf(UserExpression):
    def __init__(self, normal, mesh, **kwargs):
        super().__init__(kwargs)
        self.normal = normal
        self.mesh = mesh

    def eval_cell(self, values, x, ufc_cell):
        cell = Cell(self.mesh, ufc_cell.index)
        n = cell.normal(ufc_cell.local_facet)
        N = np.array([
                [-n[0], 0],
                [-0.5 * n[1], -0.5 * n[0]],
                [-0.5 * n[1], -0.5 * n[0]],
                [0, -n[1]],
            ])
        b = np.zeros((4, 2))
        a = np.zeros((2, 2))
        values = np.block([[N.dot(N.T), b], [b.T, a]]).reshape(-1)

    def value_shape(self):
        return ((6, 6))

class Elastic(FriedrichsSystemProblem):
    def __init__(self, L, **kwargs):
        FriedrichsSystemProblem.__init__(self, L, **kwargs)
        # ... and also store FEniCS data structures for assembly
        assert "subdomains" in kwargs
        assert "boundaries" in kwargs
        self.subdomains, self.boundaries, self.mesh = kwargs[
            "subdomains"], kwargs["boundaries"], kwargs["mesh"]

        # Defining the function spaces
        X = TensorFunctionSpace(self.mesh, "CG", 1, shape=(6, 6))

        self.y = TestFunction(L)
        self.z = TrialFunction(L)

        self.dx = Measure("dx")(domain=self.mesh,
                                subdomain_data=self.subdomains)
        self.dS = Measure("dS")(domain=self.mesh,
                                subdomain_data=self.subdomains)
        self.ds = Measure("ds")(domain=self.mesh,
                                subdomain_data=self.boundaries)

        # self.n = FacetNormal(self.mesh)
        boundaryMesh = BoundaryMesh(self.mesh, 'exterior')
        self.normalCustom = self.get_facet_normal(boundaryMesh)
        self.n = self.get_internal_normal(mesh)
        self.alpha = Constant(1000.)

        self.A_0 = Constant((('1', '0', '0', '0', '0', '0'),
            ('0', '1', '0', '0', '0', '0'),
             ('0', '0', '1', '0', '0', '0'),
             ('0', '0', '0', '1', '0', '0'),
             ('0', '0', '0', '0', '0', '0'),
             ('0', '0', '0', '0', '0', '0')))
        self.A_0_mu = Constant((('0', '0', '0', '0', '0', '0'),
            ('0', '0', '0', '0', '0', '0'),
             ('0', '0', '0', '0', '0', '0'),
             ('0', '0', '0', '0', '0', '0'),
             ('0', '0', '0', '0', '1', '0'),
             ('0', '0', '0', '0', '0', '1')))
        self.A_0_lambda = Constant(
            (('-1', '0', '0', '-1', '0', '0'),
            ('0', '0', '0', '0', '0', '0'),
             ('0', '0', '0', '0', '0', '0'),
             ('-1', '0', '0', '-1', '0', '0'),
             ('0', '0', '0', '0', '0', '0'),
             ('0', '0', '0', '0', '0', '0')))

        self.A_1 = Constant(
            (('0', '0', '0', '0', '-1', '0'),
             ('0', '0', '0', '0', '0', '-0.5'),
             ('0', '0', '0', '0', '0', '-0.5'),
             ('0', '0', '0', '0', '0', '0'),
             ('-1', '0', '0', '0', '0', '0'),
             ('0', '-0.5', '-0.5', '0', '0', '0')))
        self.A_2 = Constant(
            (('0', '0', '0', '0', '0', '0'),
            ('0', '0', '0', '0', '-0.5', '0'),
             ('0', '0', '0', '0', '-0.5', '0'),
             ('0', '0', '0', '0', '0', '-1'),
             ('0', '-0.5', '-0.5', '0', '0', '0'),
             ('0', '0', '0', '-1', '0', '0')))

        self.S_i = InterfaceSi(self.n, self.mesh, element=X.ufl_element())

        # Neumann Homogeneous
        self.S_f = InterfaceSf(self.n, self.mesh, element=X.ufl_element())

        # Dirichlet Homogeneous
        # self.S_f = Expression(
        #     (('0', '0', '0', '0', '0', '0'),
        #      ('0', '0', '0', '0', '0', '0'),
        #      ('0', '0', '0', '0', '0', '0'),
        #      ('0', '0', '0', '0', '0', '0'),
        #      ('0', '0', '0', '0', '1', '0'),
        #      ('0', '0', '0', '0', '0', '1')),
        #     element=X.ufl_element())

        self.D = self.A_1 * self.n[0] + self.A_2 * self.n[1]

        self.M = BoundaryM(normal=self.normalCustom, element=X.ufl_element())
        self.DpM = BoundaryDpM(normal=self.normalCustom, element=X.ufl_element())

        self.f = Expression(
            ("0.0", "0.0", "0.0", "0.0", "0.1*exp(-pow((x[0]-0.5), 2)-pow((x[1]-0.5), 2))", "0.1*exp(-pow((x[0]-0.5), 2)-pow((x[1]-0.5), 2))"),
            element=L.ufl_element())
        # self.f = Expression(("0", "0", "0", "0", "1", "1"),
        #                     element=L.ufl_element())

    def get_facet_normal(self, bmesh):
        '''Manually calculate FacetNormal function'''

        # if not bmesh.type().dim() == 2:
        #     raise ValueError('Only works for 2-D mesh')

        vertices = bmesh.coordinates()
        cells = bmesh.cells()

        # find normals at each edge; the orientations maybe wrong
        vec = (vertices[cells[:, 1]] - vertices[cells[:, 0]])
        # print(vec.shape)
        normals = np.ones(vec.shape)
        normals[:, 0] = vec[:, 1]
        normals[:, 1] = -vec[:, 0]
        origins = (vertices[cells[:, 1]] + vertices[cells[:, 0]]) / 2
        normals /= np.sqrt((normals**2).sum(axis=1))[:, np.newaxis]
        # plt.quiver(*origins.T, normals[:, 0], normals[:, 1])
        # plt.show()
        # plt.quiver(*origins.T, normals[:, 0], normals[:, 1])

        # Ensure outward pointing normal
        bmesh.init_cell_orientations(
            Expression(('x[0]-0.5', 'x[1]-0.5'), degree=1))
        mask = np.array(bmesh.cell_orientations()) == 0
        normals[mask, :] = -1 * normals[mask, :]
        # plt.quiver(*origins.T, normals[:, 0], normals[:, 1])
        # plt.scatter(origins[:, 0], origins[:, 1], c=bmesh.cell_orientations())
        # plt.show()

        V = VectorFunctionSpace(bmesh, 'DG', 0)
        norm = Function(V)
        nv = norm.vector()

        for n in (0, 1):
            dofmap = V.sub(n).dofmap()
            for i in np.arange(dofmap.global_dimension()):
                dof_indices = dofmap.cell_dofs(i)
                assert len(dof_indices) == 1
                nv[dof_indices[0]] = normals[i, n]

        return norm

    def get_internal_normal(self, bmesh):

        # if not bmesh.type().dim() == 2:
        #     raise ValueError('Only works for 2-D mesh')

        vertices = bmesh.coordinates()
        cells = bmesh.cells()

        # find normals at each edge; the orientations maybe wrong
        vec = (vertices[cells[:, 1]] - vertices[cells[:, 0]])
        # print(vec.shape)
        normals = np.ones(vec.shape)
        normals[:, 0] = vec[:, 1]
        normals[:, 1] = -vec[:, 0]
        origins = (vertices[cells[:, 1]] + vertices[cells[:, 0]]) / 2
        normals /= np.sqrt((normals**2).sum(axis=1))[:, np.newaxis]
        plt.quiver(*origins.T, normals[:, 0], normals[:, 1])
        plt.show()
        plt.quiver(*origins.T, normals[:, 0], normals[:, 1])

        V = VectorFunctionSpace(bmesh, 'DG', 0)
        norm = Function(V)
        nv = norm.vector()

        for n in (0, 1):
            dofmap = V.sub(n).dofmap()
            for i in np.arange(dofmap.global_dimension()):
                dof_indices = dofmap.cell_dofs(i)
                assert len(dof_indices) == 1
                nv[dof_indices[0]] = normals[i, n]

        return norm

    def compute_theta(self, term):
        mu = self.mu
        if term == "a":
            theta_a0 = 1.
            theta_a0_mu = mu[2]
            theta_a0_lambda = mu[0]
            theta_a123 = mu[1]
            theta_D = mu[1]
            theta_DpM = 0.
            theta_S_i = 0.
            theta_S_f = 0.
            return (theta_a0, theta_a0_mu, theta_a0_lambda, theta_a123, theta_S_i, theta_D,
                    theta_DpM, theta_S_f)
        elif term == "f":
            theta_f0 = 1.
            return (theta_f0, )
        else:
            raise ValueError("Invalid term for compute_theta().")

    def assemble_operator(self, term):
        z = self.z
        y = self.y
        if term == "a":
            a_a0 = inner(z, dot(self.A_0, y)) * self.dx
            a_a0_mu = inner(z, dot(self.A_0_mu, y)) * self.dx
            a_a0_lambda = inner(z, dot(self.A_0_lambda, y)) * self.dx

            a_a123 = -inner(z, dot(self.A_1, y.dx(0))) * self.dx - inner(
                z, dot(self.A_2, y.dx(1))) * self.dx

            a_D = inner(dot(self.D('+'), avg(z)), jump(y)) * self.dS
            a_DpM = 0.5 * inner(dot(self.D + self.M, z), y) * self.ds
            # a_DpM = 0.5 * inner(dot(self.DpM, z), y) * self.ds

            # stability terms (upwind): interfaces aS_i and boundaries aS_f
            aS_i = self.alpha('+') * inner(dot(self.S_i('+'), jump(z)), jump(y)) * self.dS
            aS_f = self.alpha('+') * inner(dot(self.S_f, z), y) * self.ds

            return (a_a0, a_a0_mu, a_a0_lambda, a_a123, aS_i, a_D, a_DpM, aS_f)
        elif term == "f":
            f = self.f
            f0 = dot(y, f) * self.dx
            return (f0, )
        elif term == "inner_product":
            x0 = inner(grad(z), grad(y)) * self.dx
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
L = VectorFunctionSpace(mesh, "DG", 1, dim=6)

# 3. Allocate an object of the Friedrichs' systems class
problem = Elastic(L, mesh=mesh, subdomains=subdomains, boundaries=boundaries)
mu_range = [(10, 100), (50, 100), (1e-4, 1e-3)]
problem.set_mu_range(mu_range)

# 4. Prepare reduction with a POD-Galerkin method
reduction_method = PODGalerkin(problem)
reduction_method.set_Nmax(20)

# 5. Perform the offline phase
reduction_method.initialize_training_set(1, sampling=UniformDistribution())
reduced_problem = reduction_method.offline()

# 6. Perform an online solve
online_mu = (10, 50, 1e-4)
# print("online", len(mu_range), len(online_mu))
reduced_problem.set_mu(online_mu)
reduction_method.truth_problem.set_mu(online_mu)

# Plot solution
s = reduction_method.truth_problem.solve()
p = plot(s.sub(4))
plt.colorbar(p)
plt.title("x displacement mu={}".format(online_mu[0]))
plt.show()
p = plot(s.sub(5))
plt.colorbar(p)
plt.title("y displacement")
plt.show()

u_1 = s.sub(0)
u_2 = s.sub(1)
u_3 = s.sub(2)
u_4 = s.sub(3)
u_5 = s.sub(4)
u_6 = s.sub(5)

file = File("Elastic_1.pvd")
file << u_1
file = File("Elastic_2.pvd")
file << u_2
file = File("Elastic_3.pvd")
file << u_3
file = File("Elastic_4.pvd")
file << u_4
file = File("Elastic_5.pvd")
file << u_5
file = File("Elastic_6.pvd")
file << u_6
print("Saved online snapshot")

# reduced_problem.solve()
# reduced_problem.export_solution(filename="online_solution")

# # 7. Perform an error analysis
# reduction_method.initialize_testing_set(100, sampling=UniformDistribution())

# reduction_method.error_analysis()

# # 8. Perform a speedup analysis
# reduction_method.speedup_analysis()