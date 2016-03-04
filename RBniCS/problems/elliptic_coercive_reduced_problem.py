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
## @file elliptic_coercive_reduced_problem.py
#  @brief Implementation of projection based reduced order models for elliptic coervice problems: base class
#
#  @author Francesco Ballarin <francesco.ballarin@sissa.it>
#  @author Gianluigi Rozza    <gianluigi.rozza@sissa.it>
#  @author Alberto   Sartori  <alberto.sartori@sissa.it>

from __future__ import print_function
import types
from RBniCS.elliptic_coercive_problem import EllipticCoerciveProblem

#~~~~~~~~~~~~~~~~~~~~~~~~~     ELLIPTIC COERCIVE REDUCED ORDER MODEL BASE CLASS     ~~~~~~~~~~~~~~~~~~~~~~~~~# 
## @class EllipticCoerciveReducedOrderModelBase
#
# Base class containing the interface of a projection based ROM
# for elliptic coercive problems.
class EllipticCoerciveReducedProblem(EllipticCoerciveProblem):
    
    ###########################     CONSTRUCTORS     ########################### 
    ## @defgroup Constructors Methods related to the construction of the reduced order model object
    #  @{
    
    ## Default initialization of members.
    def __init__(self, truth_problem):
        # Call to parent
        EllipticCoerciveProblem.__init__(self)
        
        self.name = truth_problem.name
        self.current_stage=None
        
        # $$ ONLINE DATA STRUCTURES $$ #
        # 1. Online reduced space dimension
        self.N = 0
        self.N_bc = 0
        # 3a. Number of terms in the affine expansion
        self.Qa = 0
        self.Qf = 0
        # 3b. Theta multiplicative factors of the affine expansion
        self.theta_a = tuple()
        self.theta_f = tuple()
        # 3c. Reduced order operators
        self.operator_a = AffineExpansionOnlineStorage()
        self.operator_f = AffineExpansionOnlineStorage()
        # Solution
        self._solution = OnlineVector()
        self._output = 0
        self._compute_error.__func__.previous_mu = None
        
        # $$ OFFLINE DATA STRUCTURES $$ #
        # 3. High fidelity problem
        self.truth_problem = truth_problem
        # 6. Basis functions matrix
        self.Z = BasisFunctionsMatrix()
        # 9. I/O
        self.basis_folder = "basis"
        self.reduced_operators_folder = "reduced_operators"
        
    #  @}
    ########################### end - CONSTRUCTORS - end ########################### 
    
    ###########################     SETTERS     ########################### 
    ## @defgroup Setters Set properties of the reduced order approximation
    #  @{
    
    ## OFFLINE/ONLINE: set the current value of the parameter. Overridden to propagate to truth problem.
    def setmu(self, mu):
        self.mu = mu
        self.truth_problem.setmu(mu)
    
    #  @}
    ########################### end - SETTERS - end ########################### 
    
    ###########################     ONLINE STAGE     ########################### 
    ## @defgroup OnlineStage Methods related to the online stage
    #  @{
    
    ## Initialize data structures required for the online phase
    def init(self, current_stage="online"):
        self.current_stage = current_stage
        if current_stage == "online":
            self.assemble_operator("a")
            self.assemble_operator("f")
            self.Qa = len(self.operator_a)
            self.Qf = len(self.operator_f)
            # Also load basis functions
            self.Z.load(self.basis_folder, "basis")
            # To properly initialize N and N_bc, detect how many theta terms
            # are related to boundary conditions
            try:
                theta_bc = self.compute_theta("dirichlet_bc")
            except RuntimeError: # there were no Dirichlet BCs
                self.N = len(self.Z)
            else: # there were Dirichlet BCs
                if not theta_bc or theta_bc.count(0.) == len(theta_bc):
                    self.N = len(self.Z)
                else:
                    self.N = len(self.Z) - len(theta_bc)
                    self.N_bc = len(theta_bc)
        elif current_stage == "offline":
            self.Qa = self.truth_problem.Qa
            self.Qf = self.truth_problem.Qf
            self.operator_a = AffineExpansionOnlineStorage(self.Qa)
            self.operator_f = AffineExpansionOnlineStorage(self.Qf)
            # Store the lifting functions in self.Z
            self.assemble_operator("dirichlet_bc")
        else:
            raise RuntimeError("Invalid stage in init().")
            
    # Perform an online solve. self.N will be used as matrix dimension if the default value is provided for N.
    def solve(self, N=None, with_plot=True):
        self.init()
        if N is None:
            N = self.N
        N += self.N_bc
        uN = self._solve(N)
        reduced_solution = self.Z*uN
        if with_plot == True:
            self._plot(reduced_solution, title = "Reduced solution. mu = " + str(self.mu), interactive = True)
        return reduced_solution
    
    # Perform an online solve (internal)
    def _solve(self, N):
        self.theta_a = self.compute_theta("a")
        self.theta_f = self.compute_theta("f")
        try:
            theta_bc = self.compute_theta("dirichlet_bc")
        except RuntimeError: # there were no Dirichlet BCs
            theta_bc = ()
        assembled_operator_a = sum(product(self.theta_a, self.operator_a[:N, :N]))
        assembled_operator_f = sum(product(self.theta_f, self.operator_f[:N]))
        solve(assembled_operator_a == assembled_operator_f, self._solution, self.theta_bc)
        return self._solution
        
    # Perform an online evaluation of the (compliant) output
    def output(self):
        N = self._solution.size
        self.theta_f = self.compute_theta("f")
        assembled_operator_f = sum(product(self.theta_f, self.operator_f[:N]))
        self._output = transpose(assembled_operator_f)*self._solution
        return self._output
        
    #  @}
    ########################### end - ONLINE STAGE - end ########################### 

    ###########################     OFFLINE STAGE     ########################### 
    ## @defgroup OfflineStage Methods related to the offline stage
    #  @{
        
    ## Assemble the reduced order affine expansion.
    def build_reduced_operators(self):
        self.assemble_operator("a")
        self.assemble_operator("f")
        
    ## Postprocess a snapshot before adding it to the basis/snapshot matrix, for instance removing
    # non-homogeneous Dirichlet boundary conditions
    def postprocess_snapshot(self, snapshot):
        try:
            theta_bc = self.compute_theta("dirichlet_bc")
        except RuntimeError: # there were no Dirichlet BCs
            pass # nothing to be done
        else: # there were Dirichlet BCs
            assert N_bc == len(theta_bc)
            snapshot -= self.Z[:N_bc]*theta_bc
        
    #  @}
    ########################### end - OFFLINE STAGE - end ########################### 
    
    ###########################     ERROR ANALYSIS     ########################### 
    ## @defgroup ErrorAnalysis Error analysis
    #  @{
    
    # Compute the error of the reduced order approximation with respect to the full order one
    # for the current value of mu
    def compute_error(self, N=None):
        if self._compute_error.__func__.previous_mu != self.mu:
            self.truth_problem.solve()
            self.truth_problem.output()
            # Do not carry out truth solves anymore for the same parameter
            self._compute_error.__func__.previous_mu = self.mu
        # Compute the error on the solution
        reduced_solution = self.solve(N, False)
        reduced_solution -= self.truth_problem._solution # store the error as a function in the reduced solution
        error = reduced_solution
        error_norm_squared = self.compute_scalar_product(error, self._error_inner_product_matrix(), error) # norm of the error
        # Compute the error on the output
        error_output = abs(self.truth_problem._output - self.online_output())
        return (sqrt(error_norm_squared), error_output)
        
    # Internal method for error computation: returns the inner product matrix to be used.
    def _error_inner_product_matrix(self):
        self.theta_a = self.compute_theta("a") # not really necessary, for symmetry with the parabolic case
        assembled_operator_a = sum(product(self.theta_a, self.operator_a)) # use the energy norm (skew part will discarded by the scalar product)
        return assembled_operator_a
        
    #  @}
    ########################### end - ERROR ANALYSIS - end ########################### 
    
    ###########################     PROBLEM SPECIFIC     ########################### 
    ## @defgroup ProblemSpecific Problem specific methods
    #  @{

    ## Return theta multiplicative terms of the affine expansion of the problem.
    def compute_theta(self, term):
        return self.truth_problem.compute_theta(theta)
        
    ## Assemble the reduced order affine expansion
    def assemble_operator(self, term):
        if self.current_stage == "online": # load from file
            # Note that we return the loaded operator (even though we do not need to
            # catch this return value in init()) because we want this interface
            # to be compatible with the one in EllipticCoerciveProblem, i.e.
            # we would like to be able to use a reduced problem also as a 
            # truth problem
            if term == "a":
                self.operator_a.load(self.reduced_operators_folder, "operator_a")
                return self.operator_a
            elif term == "f":
                self.operator_f.load(self.reduced_operators_folder, "operator_f")
                return self.operator_f
            elif term == "dirichlet_bc":
                raise RuntimeError("There should be no need to assemble Dirichlet BCs when querying online reduced problems.")
            else:
                raise RuntimeError("Invalid term for assemble_operator().")
        elif self.current_stage == "offline":
            # There is no need to return anything because the previous remark cannot hold here
            # (we are still training the reduced order model, we cannot possibly use it 
            #  anywhere else)
            if term == "a":
                for qa in range(self.Qa):
                    self.operator_a[qa] = transpose(self.Z)*self.truth_problem.operator_a[qa]*self.Z
                self.operator_a.save(self.reduced_operators_folder, "operator_a")
            elif term == "f":
                for qf in range(self.Qf):
                    self.operator_f[qf] = transpose(self.Z)*self.truth_problem.operator_f[qf]
                self.operator_f.save(self.reduced_operators_folder, "operator_f")
            elif term == "dirichlet_bc":
                try:
                    theta_bc = self.compute_theta("dirichlet_bc")
                except RuntimeError: # there were no Dirichlet BCs
                    return
                Q_dirichlet_bcs = len(theta_bc)
                # By convention, an homogeneous Dirichlet BC has all theta terms equal to 0.
                # In this case, no additional basis functions will need to be added.
                if theta_bc.count(0.) == Q_dirichlet_bcs:
                    return
                # Temporarily override compute_theta method to return only one nonzero 
                # theta term related to boundary conditions
                standard_compute_theta = self.truth_problem.compute_theta
                for i in range(Q_dirichlet_bcs):
                    def modified_compute_theta(self, term):
                        if term == "dirichlet_bc":
                            modified_theta_bc = ()
                            for j in range(Q_dirichlet_bcs):
                                if j != i:
                                    modified_theta_bc += (0.,)
                                else:
                                    modified_theta_bc += (theta_bc[i],)
                            return modified_theta
                        else:
                            return standard_compute_theta()
                    self.truth_problem.compute_theta = types.MethodType(modified_compute_theta, self.truth_problem)
                    # ... and store the solution of the truth problem corresponding to that boundary condition
                    # as lifting function
                    print("Computing and storing lifting function n.", i, " in the basis matrix")
                    lifting = self.truth_problem.solve()
                    lifting.vector()[:] /= theta_bc[i]
                    self.Z.enrich(lifting)
                # Restore the standard compute_theta method
                self.truth_problem.compute_theta = standard_compute_theta
                # Save basis functions matrix, that contains up to now only lifting functions
                self.Z.save(self.basis_folder, "basis")
                self.N_bc = Q_dirichlet_bcs
                # Note that, however, self.N is not increased, so it will actually contain the number
                # of basis functions without the lifting ones
            else:
                raise RuntimeError("Invalid term for assemble_operator().")
        else:
            raise RuntimeError("Invalid stage in assemble_operator().")
    
    ## Return a lower bound for the coercivity constant
    def get_alpha_lb(self):
        return self.truth_problem.get_alpha_lb()
                    
    #  @}
    ########################### end - PROBLEM SPECIFIC - end ########################### 

