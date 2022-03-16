# Copyright (C) 2015-2021 by the RBniCS authors
#
# This file is part of RBniCS.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

from rbnics.problems.friedrichs_systems.fs_problem import FriedrichsSystemProblem
from rbnics.reduction_methods.base import DifferentialProblemReductionMethod, LinearPODGalerkinReduction
from rbnics.reduction_methods.friedrichs_systems.friedrichs_systems_reduction_method import FriedrichsSystemsReductionMethod
from rbnics.utils.decorators import ReductionMethodFor

FriedrichsSystemsPODGalerkinReduction_Base = LinearPODGalerkinReduction(
    FriedrichsSystemsReductionMethod(DifferentialProblemReductionMethod))


# Base class containing the interface of a POD-Galerkin ROM
# for Friedrichs' systems problems
@ReductionMethodFor(FriedrichsSystemProblem, "PODGalerkin")
class FriedrichsSystemsPODGalerkinReduction(FriedrichsSystemsPODGalerkinReduction_Base):
    pass
