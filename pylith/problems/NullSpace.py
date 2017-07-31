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

## @file pylith/problems/NullSpace.py
##
## @brief Python object for creating the null space for the solver.
##
## This object does not construct the null space. Omitting the null
## space is useful if the problem is a small, overconstrained problem.
##
## Factory: null_space

from pylith.utils.PetscComponent import PetscComponent
from problems import NullSpace as ModuleNullSpace

# NullSpace class
class NullSpace(PetscComponent, ModuleNullSpace):
  """
  Python object for creating the null space for the solver.

  Factory: null_space
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="nullspace"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="null_space")
    self._createModuleObject()
    return


  def verifyConfiguration(self, materials):
    """
    Verify configuration
    """
    return

  def initialize(self, materials):
    """
    Initialize object.
    """
    return

  # PRIVATE METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    PetscComponent._configure(self)
    return


  def _createModuleObj(self):
    """
    Create handle to C++ NullSpace.
    """
    ModuleNullSpace.__init__(self)
    return
    
  
# FACTORIES ////////////////////////////////////////////////////////////

def null_space():
  """
  Factory associated with NullSpace..
  """
  return NullSpace()


# End of file 
