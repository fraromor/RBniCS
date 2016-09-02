# Copyright (C) 2015-2016 by the RBniCS authors
#
# This file is part of RBniCS.
#
# RBniCS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# RBniCS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with RBniCS. If not, see <http://www.gnu.org/licenses/>.
#
## @file scm.py
#  @brief Implementation of the successive constraints method for the approximation of the coercivity constant
#
#  @author Francesco Ballarin <francesco.ballarin@sissa.it>
#  @author Gianluigi Rozza    <gianluigi.rozza@sissa.it>
#  @author Alberto   Sartori  <alberto.sartori@sissa.it>

from math import sqrt
from dolfin import adjoint, Function, DirichletBC
from RBniCS.problems.base import ParametrizedProblem
from RBniCS.backends import AffineExpansionStorage, EigenSolver, sum, product
from RBniCS.utils.decorators import sync_setters, Extends, override

@Extends(ParametrizedProblem)
class ParametrizedHermitianEigenProblem(ParametrizedProblem):
    ###########################     CONSTRUCTORS     ########################### 
    ## @defgroup Constructors Methods related to the construction of the EIM object
    #  @{

    ## Default initialization of members
    @override
    @sync_setters("truth_problem", "set_mu", "mu")
    @sync_setters("truth_problem", "set_mu_range", "mu_range")
    def __init__(self, truth_problem, term, multiply_by_theta, constrain_eigenvalue, spectrum, eigensolver_parameters):
        # Call the parent initialization
        ParametrizedProblem.__init__(self, folder_prefix="") # this class does not export anything
        self.truth_problem = truth_problem
        
        # We need to discard dofs related to bcs in eigenvalue computations. To avoid having to create a PETSc submatrix
        # we simply zero rows and columns and replace the diagonal element with an eigenvalue that for sure
        # will not be the one we are interested in
        self.constrain_eigenvalue = constrain_eigenvalue
        # Matrices/vectors resulting from the truth discretization: condensed version discard
        # Dirichlet DOFs
        self.term = term
        assert isinstance(self.term, (tuple, str))
        if isinstance(self.term, tuple):
            assert len(self.term) == 2
            isinstance(self.term[0], str)
            isinstance(self.term[1], int)
        self.multiply_by_theta = multiply_by_theta
        assert isinstance(self.multiply_by_theta, bool)
        self.operator__condensed = AffineExpansionStorage()
        self.inner_product__condensed = AffineExpansionStorage() # even though it will contain only one matrix
        self.spectrum = spectrum
        self.eigensolver_parameters = eigensolver_parameters
        
        # Avoid useless computations
        self._solve__previous_mu = None
        self._solve__previous_eigenvalue = None
        self._solve__previous_eigenvector = None
        
    #  @}
    ########################### end - CONSTRUCTORS - end ###########################
    
    def init(self):
        # Condense the symmetric part of the required term
        if isinstance(self.term, tuple):
            forms = (self.truth_problem.assemble_operator(self.term[0])[ self.term[1] ], )
        else:
            assert isinstance(self.term, str)
            forms = self.truth_problem.assemble_operator(self.term)
        symmetric_forms = [ 0.5*(form + adjoint(form)) for form in forms]
        symmetric_forms = tuple(symmetric_forms)
        self.operator__condensed = AffineExpansionStorage(symmetric_forms)
        self.clear_constrained_dofs(self.operator__condensed, self.constrain_eigenvalue)
        
        # Condense the inner product matrix
        self.inner_product__condensed = AffineExpansionStorage(self.truth_problem.assemble_operator("inner_product"))
        self.clear_constrained_dofs(self.inner_product__condensed, 1.)
        
    # Clear constrained dofs
    def clear_constrained_dofs(self, operator, diag_value):
        dirichlet_bc = self.truth_problem.dirichlet_bc
        V = self.truth_problem.V
        for op in operator:
            if len(dirichlet_bc) > 0:
                dummy = Function(V)
                for bc_list in dirichlet_bc:
                    for bc in bc_list:
                        bc.zero(op)
                        bc.zero_columns(op, dummy.vector(), diag_value)
    
    def solve(self):
        if self._solve__previous_mu == self.mu:
            return (self._solve__previous_eigenvalue, self._solve__previous_eigenvector)
        else:
            if self.multiply_by_theta:
                assert isinstance(self.term, str) # method untested otherwise
                O = sum(product(self.truth_problem.compute_theta(self.term), self.operator__condensed))
            else:
                assert isinstance(self.term, tuple) # method untested otherwise
                theta = (1.,)
                assert len(theta) == len(self.operator__condensed)
                O = sum(product(theta, self.operator__condensed))
            assert len(self.inner_product__condensed) == 1
            X = self.inner_product__condensed[0]
            
            eigensolver = EigenSolver(O, X, self.truth_problem.V)
            eigensolver_parameters = dict()
            eigensolver_parameters["problem_type"] = "gen_hermitian"
            assert self.spectrum is "largest" or self.spectrum is "smallest"
            eigensolver_parameters["spectrum"] = self.spectrum + " real"
            if self.eigensolver_parameters is not None:
                eigensolver_parameters.update(self.eigensolver_parameters)
            eigensolver.set_parameters(eigensolver_parameters)
            eigensolver.solve(1)
            
            r, c = eigensolver.get_eigenvalue(0) # real and complex part of the eigenvalue
            r_vector, c_vector = eigensolver.get_eigenvector(0) # real and complex part of the eigenvectors
            
            from numpy import isclose
            assert isclose(c, 0), "The required eigenvalue is not real"
            assert not isclose(r, self.constrain_eigenvalue), "The required eigenvalue is too close to the one used to constrain Dirichlet boundary conditions"
            #assert r >= 0 or isclose(r, 0), "The required eigenvalue is not positive"
            
            r_vector = Function(self.truth_problem.V, r_vector)
            
            self._solve__previous_mu = self.mu
            self._solve__previous_eigenvalue = r
            self._solve__previous_eigenvector = r_vector
            
            return (r, r_vector)
        
