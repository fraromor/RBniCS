# Copyright (C) 2015-2017 by the RBniCS authors
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

from rbnics.problems.base import NonlinearReducedProblem
from rbnics.problems.base import ParametrizedReducedDifferentialProblem
from navier_stokes_tensor3_problem import NavierStokesProblem
from rbnics.backends import product, sum
from rbnics.backends.online import OnlineFunction
from rbnics.utils.decorators import Extends, override

def NavierStokesReducedProblem(StokesReducedProblem_DerivedClass):
    
    NavierStokesReducedProblem_Base = NonlinearReducedProblem(StokesReducedProblem_DerivedClass)
    
    @Extends(NavierStokesReducedProblem_Base)
    class NavierStokesReducedProblem_Class(NavierStokesReducedProblem_Base):
        
        class ProblemSolver(NavierStokesReducedProblem_Base.ProblemSolver):
            def residual_eval(self, solution):
                self.store_solution(solution)
                problem = self.problem
                N = self.N
                assembled_operator = dict()
                for term in ("a", "b", "bt", "c", "f", "g"):
                    assert self.terms_order[term] in (1, 2)
                    if self.terms_order[term] == 2:
                        assembled_operator[term] = sum(product(problem.compute_theta(term), problem.operator[term][:N, :N]))
                    elif self.terms_order[term] == 1:
                        assembled_operator[term] = sum(product(problem.compute_theta(term), problem.operator[term][:N]))
                    else:
                        raise AssertionError("Invalid value for order of term " + term)
                return (
                      assembled_operator["a"] + assembled_operator["c"]
                    + assembled_operator["b"] + assembled_operator["bt"]
                    - assembled_operator["f"] - assembled_operator["g"]
                )
                
            def jacobian_eval(self, solution):
                self.store_solution(solution)
                problem = self.problem
                N = self.N
                assembled_operator = dict()
                for term in ("da", "db", "dbt", "dc"):
                    assert self.terms_order[term] is 2
                    assembled_operator[term] = sum(product(problem.compute_theta(term), problem.operator[term][:N, :N]))
                return (
                      assembled_operator["da"] + assembled_operator["dc"]
                    + assembled_operator["db"] + assembled_operator["dbt"]
                )
        
    # return value (a class) for the decorator
    return NavierStokesReducedProblem_Class

