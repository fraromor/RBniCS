# Copyright (C) 2015-2022 by the RBniCS authors
#
# This file is part of RBniCS.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

def function_extend_or_restrict(function, function_components, V, V_components, weight, copy_,
                                extended_or_restricted_function=None):

    from rbnics.backends.online.numpy.copy import copy
    if copy_:
        return copy(function)
    else:
        return function
