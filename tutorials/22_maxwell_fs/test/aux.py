from dolfin import *
from rbnics import *
import numpy as np

class BoundaryM(UserExpression):
    def __init__(self, normal, **kwargs):
        super().__init__(kwargs)
        self.normal = normal

    def eval_cell(self, values, x, ufc_cell):
        if np.linalg.norm(np.array(x[:2])) + DOLFIN_EPS >= 3 or x[2] + DOLFIN_EPS >= 1 or x[2] <= DOLFIN_EPS:
            t = np.array([[0, -self.normal(x)[2], self.normal(x)[1]],
                          [self.normal(x)[2], 0, -self.normal(x)[0]],
                          [-self.normal(x)[1],self.normal(x)[0], 0]])
            a = np.zeros((3, 3))
            values = np.block([[a, -t], [t.T, a]])
        else:
            values = np.zeros((6, 6))

    def value_shape(self):
        return ((6, 6))

def get_facet_normal(bmesh):
        '''Manually calculate FacetNormal function of cables domain'''

        vertices = bmesh.coordinates()
        cells = bmesh.cells()

        vec1 = vertices[cells[:, 1]] - vertices[cells[:, 0]]
        vec2 = vertices[cells[:, 2]] - vertices[cells[:, 0]]

        normals = np.zeros((cells.shape[0], 3))
        normals[:, :2] = vertices[cells[:, 0]][:, :2]
        normals /= np.sqrt((normals**2).sum(axis=1))[:, np.newaxis]
        origins = (vertices[cells[:, 1]] + vertices[cells[:, 0]] +
                   vertices[cells[:, 2]]) / 3

        # mask = np.linalg.norm(vertices[cells[:, 0]], axis=1) < 3. - DOLFIN_EPS #or origins[:, 2] < DOLFIN_EPS or origins[:, 2] > 1 - DOLFIN_EPS
        # normals[mask, :] = 0
        # origins[mask, :] = 0
        normals[origins[:, 2] < DOLFIN_EPS] = np.array([0, 0, -1])
        normals[origins[:, 2] > 1 - DOLFIN_EPS] = np.array([0, 0, 1])

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
        #           length=0.2, normalize=True)
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