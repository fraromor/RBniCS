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
## @file elliptic_coercive_pod_galerkin_reduction.py
#  @brief Implementation of a POD-Galerkin ROM for elliptic coervice problems
#
#  @author Francesco Ballarin <francesco.ballarin@sissa.it>
#  @author Gianluigi Rozza    <gianluigi.rozza@sissa.it>
#  @author Alberto   Sartori  <alberto.sartori@sissa.it>

from __future__ import print_function
from numpy import log, exp, mean # for error analysis
import os # for path and makedir
from RBniCS.proper_orthogonal_decomposition import ProperOrthogonalDecomposition
from RBniCS.elliptic_coercive_base import EllipticCoerciveBase

#~~~~~~~~~~~~~~~~~~~~~~~~~     ELLIPTIC COERCIVE POD BASE CLASS     ~~~~~~~~~~~~~~~~~~~~~~~~~# 
## @class EllipticCoercivePODGalerkinReduction
#
# Base class containing the interface of a POD-Galerkin ROM
# for elliptic coercive problems
class EllipticCoercivePODGalerkinReduction(EllipticCoerciveReductionMethodBase):
    """This class implements a reduced order method based on a POD (Proper
    Orthogonal Decomposition) Galerkin approach. In particular, it
    implements the offline phase and the error analysis proper for the
    POD approach.
    
    This class provides the following methods:
    
    ##  Methods related to the offline stage
    - offline()
    - update_snapshot_matrix()
    - apply_POD()

    ## Error analysis
    - error_analysis()

    A typical usage of this class is reported in tutorial 2.

    """
    
    ###########################     CONSTRUCTORS     ########################### 
    ## @defgroup Constructors Methods related to the construction of the POD-Galerkin ROM object
    #  @{
    
    ## Default initialization of members
    def __init__(self, truth_problem):
        # Call the parent initialization
        EllipticCoerciveReductionMethodBase.__init__(self, truth_problem)
                
        # $$ OFFLINE DATA STRUCTURES $$ #
        # 6bis. Declare a POD object
        self.POD = ProperOrthogonalDecomposition(self.compute_scalar_product, self.S)
        # 9. I/O
        self.xi_train_folder = "xi_train__pod"
        self.xi_test_folder = "xi_test__pod"
        self.snapshots_folder = "snapshots__pod"
        self.post_processing_folder = "post_processing__pod"
        
    #  @}
    ########################### end - CONSTRUCTORS - end ########################### 
    
    ###########################     OFFLINE STAGE     ########################### 
    ## @defgroup OfflineStage Methods related to the offline stage
    #  @{
    
    ## Initialize data structures required for the offline phase
    def _init_offline(self):
        # Call the parent initialization
        need_to_do_offline_stage = EllipticCoerciveReductionMethodBase._init_offline(self)
        
        # Also create folders for snapshots and postprocessing
        folders = (self.snapshots_folder, self.post_processing_folder)
        for f in folders:
            if not os.path.exists(f):
                os.makedirs(f)
        
        return need_to_do_offline_stage
    
    ## Perform the offline phase of the reduced order model
    def offline(self):
        need_to_do_offline_stage = self._init_offline()
        if not need_to_do_offline_stage:
            return self.reduced_problem
        
        print("==============================================================")
        print("=             Offline phase begins                           =")
        print("==============================================================")
        print("")
        
        for run in range(len(self.xi_train)):
            print("############################## run = ", run, " ######################################")
            
            self.truth_problem.setmu(self.xi_train[run])
            
            print("truth solve for mu = ", self.mu)
            snapshot = self.truth_problem.solve()
            self.truth_problem.export_solution(snapshot, self.snapshots_folder, "truth_" + str(run))
            self.reduced_problem.postprocess_snapshot(snapshot)
            
            print("update snapshot matrix")
            self.update_snapshot_matrix(snapshot)

            print("")
            run += 1
            
        print("############################## perform POD ######################################")
        (Z, N) = self.POD.apply(self.Nmax)
        self.reduced_problem.Z.enrich(Z)
        self.reduced_problem.N += N
        self.reduced_problem.Z.save(self.basis_folder, "basis")
        self.POD.print_eigenvalues()
        self.POD.save_eigenvalues_file(self.post_processing_folder, "eigs")
        self.POD.save_retained_energy_file(self.post_processing_folder, "retained_energy")
        
        print("")
        print("build reduced operators")
        self.reduced_problem.build_reduced_operators()
        
        print("")
        print("==============================================================")
        print("=             Offline phase ends                             =")
        print("==============================================================")
        print("")
        
        return self.reduced_problem
        
    ## Update the snapshot matrix
    def update_snapshot_matrix(self, snapshot):
        self.POD.store_snapshot(snapshot)
        
    #  @}
    ########################### end - OFFLINE STAGE - end ########################### 
    
    ###########################     ERROR ANALYSIS     ########################### 
    ## @defgroup ErrorAnalysis Error analysis
    #  @{
    
    # Compute the error of the reduced order approximation with respect to the full order one
    # over the test set
    def error_analysis(self, N=None):
        if N is None:
            N = self.reduced_problem.N
            
        print("==============================================================")
        print("=             Error analysis begins                          =")
        print("==============================================================")
        print("")
        
        error_u = MultiIndexArray((N, len(self.xi_test)))
        error_s = MultiIndexArray((N, len(self.xi_test)))
        
        for run in range(len(self.xi_test)):
            print("############################## run = ", run, " ######################################")
            
            self.reduced_problem.setmu(self.xi_test[run])
                        
            for n in range(N): # n = 0, 1, ... N - 1
                (error_u[n, run], error_s[n, run]) = self.reduced_problem.compute_error(n + 1, True)
        
        # Print some statistics
        print("")
        print("N \t gmean(err_u)")
        for n in range(N): # n = 0, 1, ... N - 1
            mean_error_u = exp(mean(log((error_u[n, :]))))
            print(str(n+1) + " \t " + str(mean_error_u))
        
        print("")
        print("N \t gmean(err_s)")
        for n in range(N): # n = 0, 1, ... N - 1
            mean_error_s = exp(mean(log((error_s[n, :]))))
            print(str(n+1) + " \t " + str(mean_error_s))
        
        print("")
        print("==============================================================")
        print("=             Error analysis ends                            =")
        print("==============================================================")
        print("")
        
    #  @}
    ########################### end - ERROR ANALYSIS - end ########################### 
