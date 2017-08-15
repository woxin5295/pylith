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

# NullSpaceSingleBody class
class NullSpaceSingleBody(NullSpace):
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


  def initialize(self, fields, formulation):
    """
    Initialize object.
    """
    import numpy
    numBodies = 1
    numMaterials = len(self.materials.components())
    bodiesIds = numpy.zeros((numMaterials,), dtype=numpy.int32)
    numMaterialsInBodies = numpy.array([numMaterials], dtype=numpy.int32)
    for imaterial,material in enumerate(self.materials.components()):
        bodiesIds[imaterial] = material.id()

    formulation.createNullSpaceBodies(fields, numMaterialsInBodies, bodiesIds)
    return

# FACTORIES ////////////////////////////////////////////////////////////

def null_space():
  """
  Factory associated with NullSpaceSingleBody.
  """
  return NullSpaceSingleBody()


# End of file 
