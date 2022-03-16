# Copyright (C) 2015-2021 by the RBniCS authors
#
# This file is part of RBniCS.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

from rbnics.problems.base import LinearProblem, ParametrizedDifferentialProblem
from rbnics.backends import product, sum, transpose

FriedrichsSystemProblem_Base = LinearProblem(ParametrizedDifferentialProblem)


# Base class containing the definition of Friedrichs' systems problems
class FriedrichsSystemProblem(FriedrichsSystemProblem_Base):

    # Default initialization of members
    def __init__(self, V, **kwargs):
        # Call to parent
        FriedrichsSystemProblem_Base.__init__(self, V, **kwargs)

        # Form names for
        self.terms = ["a", "f", "s"]
        self.terms_order = {"a": 2, "f": 1, "s": 1}
        self.components = ["z"]

    class ProblemSolver(FriedrichsSystemProblem_Base.ProblemSolver):
        def matrix_eval(self):
            problem = self.problem
            return sum(product(problem.compute_theta("a"), problem.operator["a"]))

        def vector_eval(self):
            problem = self.problem
            return sum(product(problem.compute_theta("f"), problem.operator["f"]))
