#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/problems/Solver.py
##
## @brief Python PyLith abstract base class for solver.
##
## Factory: solver

from pylith.utils.PetscComponent import PetscComponent

# VALIDATORS ///////////////////////////////////////////////////////////

# Solver class
class Solver(PetscComponent):
  """
  Python abstract base class for solver.

  Factory: solver.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
    """
    Python object for managing Solver facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Solver facilities and properties.
    ##
    ## \b Properties
    ## @li \b use_cuda Use CUDA in solve if supported by solver.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="solver"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="solver")
    self._createModuleObj()
    return

  # PRIVATE METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    PetscComponent._configure(self)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def solver():
  """
  Factory associated with Solver.
  """
  return Solver()


# End of file 
