from dolfin import *
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
        if (abs(x[0]) < DOLFIN_EPS or abs(x[0] - 1.0) < DOLFIN_EPS) \
            or (abs(x[1]) < DOLFIN_EPS or abs(x[1] - 1.0) < DOLFIN_EPS):
            values[2] = -self.normal(x)[0]
            values[5] = -self.normal(x)[1]
            values[8] = 0 # Dirichlet
            #values[8] = inner(self.beta, self.normal)(x) # Neumann
            values[6] = self.normal(x)[0]
            values[7] = self.normal(x)[1]
            values[0] = 0
            values[1] = 0
            values[3] = 0
            values[4] = 0
        else:
            values = [0] * 9


    def value_shape(self):
        return ((3, 3))

class AdvectionDiffusionReaction(FriedrichsSystemProblem):
    def __init__(self, L, **kwargs):
        FriedrichsSystemProblem.__init__(self, L, **kwargs)
        # ... and also store FEniCS data structures for assembly
        assert "subdomains" in kwargs
        assert "boundaries" in kwargs
        self.subdomains, self.boundaries, self.mesh = kwargs[
            "subdomains"], kwargs["boundaries"], kwargs["mesh"]

        # Defining the function spaces
        X = TensorFunctionSpace(self.mesh, "DG", 0, shape=(3, 3))

        self.y = TestFunction(L)
        self.z = TrialFunction(L)

        self.dx = Measure("dx")(domain=self.mesh,
                                subdomain_data=self.subdomains)
        self.dS = Measure("dS")(domain=self.mesh,
                                subdomain_data=self.boundaries)
        self.ds = Measure("ds")(domain=self.mesh,
                                subdomain_data=self.boundaries)

        self.n = FacetNormal(self.mesh)
        boundaryMesh = BoundaryMesh(mesh, 'exterior')
        self.normalCustom = self.get_facet_normal(boundaryMesh)
        self.alpha = Constant(1)

        self.A_0_kappa = Expression(
            (('1', '0', '0'), ('0', '1', '0'), ('0', '0', '0')),
            element=X.ufl_element())
        self.A_0_nu = Expression(
            (('0', '0', '0'), ('0', '0', '0'), ('0', '0', '1')),
            element=X.ufl_element())

        self.A_1_0 = Expression(
            (("0", "0", "1"), ("0", "0", "0"), ("1", "0", "0")),
            element=X.ufl_element())
        self.A_1_1 = Expression(
            (("0", "0", "0"), ("0", "0", "0"), ("0", "0", "1")),
            element=X.ufl_element())

        self.A_2_0 = Expression(
            (("0", "0", "0"), ("0", "0", "1"), ("0", "1", "0")),
            element=X.ufl_element())
        self.A_2_1 = Expression(
            (("0", "0", "0"), ("0", "0", "0"), ("0", "0", "1")),
            element=X.ufl_element())

        self.D = (self.A_1_0 + self.A_1_1) * self.n[0] + (
            self.A_2_0 + self.A_2_1) * self.n[1]

        self.M = BoundaryM(normal=self.normalCustom,
                      element=X.ufl_element())

        self.f = Expression(
            ("0.0", "0.0", "10*exp(-pow((x[0]-0.5), 2)-pow((x[1]-0.5), 2))"),
            element=L.ufl_element())

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
        origins = (vertices[cells[:, 1]] + vertices[cells[:, 0]])/2
        normals /= np.sqrt((normals**2).sum(axis=1))[:, np.newaxis]
        # plt.quiver(*origins.T, normals[:, 0], normals[:, 1])
        # plt.show()
        # plt.quiver(*origins.T, normals[:, 0], normals[:, 1])

        # Ensure outward pointing normal
        bmesh.init_cell_orientations(Expression(('x[0]-0.5', 'x[1]-0.5'), degree=1))
        mask = np.array(bmesh.cell_orientations())== 0
        normals[mask, :] = -1 * normals[mask, :]
        # plt.quiver(*origins.T, normals[:, 0], normals[:, 1])
        # plt.scatter(origins[:, 0], origins[:, 1], c=bmesh.cell_orientations())
        # plt.show()

        V = VectorFunctionSpace(bmesh, 'DG', 0)
        norm = Function(V)
        nv = norm.vector()

        for n in (0,1):
            dofmap = V.sub(n).dofmap()
            for i in np.arange(dofmap.global_dimension()):
                dof_indices = dofmap.cell_dofs(i)
                assert len(dof_indices) == 1
                nv[dof_indices[0]] = normals[i, n]

        return norm

    def compute_theta(self, term):
        mu = self.mu
        if term == "a":
            theta_a0kappa = mu[0]#1/kappa
            theta_a0nu = mu[1]#nu
            theta_a10 = 1.
            theta_a11 = mu[2]#beta1
            theta_a20 = 1.
            theta_a21 = mu[3]#beta2
            theta_D = 1.
            theta_S = 1.
            theta_DpM = 1.
            return (theta_a0kappa, theta_a0nu, theta_a10, theta_a11, theta_a20, theta_a21, theta_D, theta_S, theta_DpM)
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
            a_a0kappa = inner(z, dot(self.A_0_kappa, y)) * dx
            a_a0nu = inner(z, dot(self.A_0_nu, y)) * dx
            a_a10 = -inner(z, dot(self.A_1_0, y.dx(0))) * dx
            a_a11 = -inner(z, dot(self.A_1_1, y.dx(0))) * dx
            a_a20 = -inner(z, dot(self.A_2_0, y.dx(1))) * dx
            a_a21 = -inner(z, dot(self.A_2_1, y.dx(1))) * dx
            aD = inner(dot(self.D('+'), avg(z)), jump(y)) * dS
            aS = self.alpha('+') * inner(jump(z), jump(y)) * dS
            aDpM = 0.5 * inner(dot(self.D + self.M, z), y) * ds
            return (a_a0kappa, a_a0nu, a_a10, a_a11, a_a20, a_a21, aD, aS, aDpM)
        elif term == "f":
            f = self.f
            f0 = dot(y, f) * dx
            return (f0, )
        elif term == "inner_product":
            x0 = inner(grad(z), grad(y)) * dx
            return (x0,)
        else:
            raise ValueError("Invalid term for assemble_operator().")


# 1. Read the mesh for this problem
mesh = Mesh("data/square.xml")
subdomains = MeshFunction("size_t", mesh, "data/square_physical_region.xml")
boundaries = MeshFunction("size_t", mesh, "data/square_facet_region.xml")

# 2. Create Finite Element space (Lagrange P1, two components)
L = VectorFunctionSpace(mesh, "DG", 1, dim=3)

# 3. Allocate an object of the Friedrichs' systems class
problem = AdvectionDiffusionReaction(L, mesh= mesh, subdomains=subdomains, boundaries=boundaries)
mu_range = [(1.0, 1e2), (1e-4, 1.0), (0, 1), (0, 1)]
problem.set_mu_range(mu_range)

# 4. Prepare reduction with a POD-Galerkin method
reduction_method = PODGalerkin(problem)
reduction_method.set_Nmax(20)

# 5. Perform the offline phase
reduction_method.initialize_training_set(100, sampling=UniformDistribution())
reduced_problem = reduction_method.offline()

# 6. Perform an online solve
online_mu = (1, 1, 0.5, 0.5)
reduced_problem.set_mu(online_mu)
reduced_problem.solve()
reduced_problem.export_solution(filename="online_solution")

# 7. Perform an error analysis
reduction_method.initialize_testing_set(100, sampling=UniformDistribution())
reduction_method.error_analysis()

# 8. Perform a speedup analysis
reduction_method.speedup_analysis()