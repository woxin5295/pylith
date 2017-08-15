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

from NullSpace import NullSpace
from pyre.components.Component import Component

# ITEM FACTORIES ///////////////////////////////////////////////////////

# RigidBody
class RigidBody(Component):

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing RigidBody facilities and properties.
    """

    import pyre.inventory
    materials = pyre.inventory.list("materials", default=[])
    materials.meta['tip'] = "Array of materials in body."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="rigidbody"):
    """
    Constructor.
    """
    Component.__init__(self, name)
    return

  def _configure(self):
    """
    Set members based using inventory.    
    """
    Component._configure(self)
    self.materials = self.inventory.materials
    return
  
def bodyFactory(name):
  """
  Factory for rigid bodies.
  """
  from pyre.inventory import facility
  return facility(name, family="rigid_body", factory=RigidBody)


# NullSpaceMultiBody class
class NullSpaceMultiBody(NullSpace):
  """
  Python object for creating the null space with multiple bodies.

  Factory: null_space
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(NullSpace.Inventory):
    """
    Python object for managing NullSpaceMultiBody facilities and properties.
    """

    import pyre.inventory
    bodies = pyre.inventory.facilityArray("bodies", itemFactory=bodyFactory)
    bodies.meta['tip'] = "Array of bodies separated by faults."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="nullspacemultibody"):
    """
    Constructor.
    """
    NullSpace.__init__(self, name)
    return


  def verifyConfiguration(self):
    """
    Verify configuration
    """
    msg = "Could not find body material(s) in materials."
    missing = False
    for body in self.bodies.components():
        for name in body.materials:
            found = False
            for material in self.materials.facilities():
                if name == material.name:
                    found = True
                    break
            if not found:
                missing = True
                msg += "\n  Material '%s' not found." % (name)
    if missing:
        raise ValueError(msg)
    return


  def initialize(self, fields, formulation):
    """
    Initialize object.
    """
    import numpy
    
    numBodies = len(self.bodies.components())
    numMaterialsInBodies = numpy.zeros((numBodies,), dtype=numpy.int32)
    numMaterials = len(self.materials.components())
    bodiesIds = numpy.zeros((numMaterials,), dtype=numpy.int32)
    index = 0
    for ibody,body in enumerate(self.bodies.components()):
        for icount,name in enumerate(body.materials):
            for material in self.materials.components():
                materialFacility = material.aliases[-1]
                if name == materialFacility:
                    bodiesIds[index] = material.id()
                    numMaterialsInBodies[ibody] += 1
                    index += 1

    formulation.createNullSpaceBodies(fields, numMaterialsInBodies, bodiesIds)
    return

  # PRIVATE METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    NullSpace._configure(self)
    self.bodies = self.inventory.bodies
    return


# FACTORIES ////////////////////////////////////////////////////////////

def null_space():
  """
  Factory associated with NullSpace.
  """
  return NullSpaceMultiBody()


# End of file 
