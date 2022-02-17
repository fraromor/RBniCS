# Copyright (C) 2015-2022 by the RBniCS authors
#
# This file is part of RBniCS.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

from mpi4py.MPI import COMM_WORLD


def get_mpi_comm(arg):
    if isinstance(arg, int):
        return COMM_WORLD
    else:
        basis_functions = arg
        return basis_functions.mpi_comm
