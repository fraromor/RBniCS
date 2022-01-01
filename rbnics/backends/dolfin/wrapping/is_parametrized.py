# Copyright (C) 2015-2022 by the RBniCS authors
#
# This file is part of RBniCS.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

from dolfin.function.expression import BaseExpression
from rbnics.backends.dolfin.wrapping.is_problem_solution import is_problem_solution
from rbnics.backends.dolfin.wrapping.is_problem_solution_dot import is_problem_solution_dot
from rbnics.backends.dolfin.wrapping.is_problem_solution_type import is_problem_solution_type
from rbnics.backends.dolfin.wrapping.pull_back_to_reference_domain import (
    is_pull_back_expression, is_pull_back_expression_parametrized)
from rbnics.utils.decorators import ModuleWrapper


def basic_is_parametrized(backend, wrapping):

    def _basic_is_parametrized(expression_or_form, iterator):
        for node in iterator(expression_or_form):
            # ... parametrized expressions
            if isinstance(node, BaseExpression):
                if is_pull_back_expression(node) and is_pull_back_expression_parametrized(node):
                    return True
                else:
                    parameters = node._parameters
                    if "mu_0" in parameters:
                        return True
            # ... problem solutions related to nonlinear terms
            elif wrapping.is_problem_solution_type(node):
                if wrapping.is_problem_solution(node) or wrapping.is_problem_solution_dot(node):
                    return True
        return False

    return _basic_is_parametrized


backend = ModuleWrapper()
wrapping = ModuleWrapper(is_problem_solution, is_problem_solution_dot, is_problem_solution_type)
is_parametrized = basic_is_parametrized(backend, wrapping)
