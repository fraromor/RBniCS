# Copyright (C) 2015-2022 by the RBniCS authors
#
# This file is part of RBniCS.
#
# SPDX-License-Identifier: LGPL-3.0-or-later


def functions_list_mul_online_matrix(functions_list, online_matrix, FunctionsListType):
    raise RuntimeError("TODO")  # TODO


def functions_list_mul_online_vector(functions_list, online_vector):
    from rbnics.backends.online.numpy.function import Function
    assert len(functions_list) > 0
    output = Function(functions_list[0].vector().N)
    for (i, fun_i) in enumerate(functions_list):
        output.vector()[:] += fun_i.vector() * online_vector[i]
    return output
