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

## @file pylith/problems/NullSpaceSingleBody.py
##
## @brief Python object for creating the null space for a domain with
## a single body.
##
## Factory: null_space

from NullSpace import NullSpace
from problems import NullSpaceBodies as ModuleNullSpaceBodies

# NullSpaceSingleBody class
class NullSpaceSingleBody(NullSpace, ModuleNullSpaceBodies):
  """
  Python object for creating the null space for a domain with a single body.

  Factory: null_space
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="nullspacesinglebody"):
    """
    Constructor.
    """
    NullSpace.__init__(self, name)
    return


  def initialize(self, materials):
    """
    Initialize object. Setup a single body containing all of the materials.
    """
    numMaterials = len(materials.components())
    numMaterialsInBodies = [numMaterials]
    numBodies = 1
    bodyIds = []
    for material in materials.components():
        bodyIds.append[material.id()]
    ModuleNullSpaceBodies.bodies(self, bodyIds, numMaterials, numMaterialsInBodies, numBodies)
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
    Create handle to C++ NullSpaceSingleBody.
    """
    ModuleNullSpaceSingleBody.__init__(self)
    return
    
  
# FACTORIES ////////////////////////////////////////////////////////////

def null_space():
  """
  Factory associated with NullSpaceSingleBody.
  """
  return NullSpaceNone()


# End of file 
