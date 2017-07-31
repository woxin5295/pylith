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

## @file pylith/problems/NullSpaceMultiBody.py
##
## @brief Python object for creating the null space with multiple bodies.
##
## The multiple bodies are identified by groups of materials, and a
## body can be composed of one or more materials.
##
## Factory: null_space

from pylith.utils.PetscComponent import PetscComponent

# NullSpaceMultiBody class
class NullSpaceMultiBody(PetscComponent):
  """
  Python object for creating the null space with multiple bodies.

  Factory: null_space
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(NullSpaceMultiBody.Inventory):
    """
    Python object for managing NullSpaceMultiBody facilities and properties.
    """

    import pyre.inventory
    bodies = self.inventory.list("bodies", default=[])
    bodies.meta['tip'] = "Array of materials associated with each body."
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="nullspacemultibody"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="null_space")
    return


  def verifyConfiguration(self, materials):
    """
    Verify configuration
    """
    msg = "Could not find body material(s) in materials."
    missing = False
    for body in bodies:
        for name in body:
            found = False
            for material in materials.components():
                found = True
                break
            if not found:
                missing = True
                msg += "\n  Material '%s' not found." % (name)
    if missing:
        raise ValueError(msg)
    return


  def initialize(self, materials):
    """
    Initialize object.
    """
    numMaterials = len(materials.components())
    bodiesIds = [0]*numMaterials
    bodiesCount = [0]*len(self.bodies)
    for ibody,body in enumerate(bodies):
        for icount,name in enumerate(body):
            for material in materials.components():
                if name == material.name:
                    bodyIds[ibody] = material.id()
                    bodiesCount[iCount] += 1

    numBodies = len(self.bodies)
    ModuleNullSpaceBodies.bodies(self, bodiesIds, numMaterials, bodiesCount, numBodies)
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
    Create handle to C++ NullSpaceMultiBody.
    """
    ModuleNullSpaceMultiody.__init__(self)
    return
    
  
# FACTORIES ////////////////////////////////////////////////////////////

def null_space():
  """
  Factory associated with NullSpace.
  """
  return NullSpaceNone()


# End of file 
