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

from pyre.components.Component import Component

# NullSpace class
class NullSpace(Component):
  """
  Python object for creating the null space for the solver.

  Factory: null_space
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="nullspace"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="null_space")
    return


  def preinitialize(self, materials):
    """
    Do minimum initialization.
    """
    self.materials = materials
    return

  
  def verifyConfiguration(self):
    """
    Verify configuration
    """
    return

  def initialize(self, fields, formulation):
    """
    Initialize object.
    """
    return

  # PRIVATE METHODS /////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def null_space():
  """
  Factory associated with NullSpace..
  """
  return NullSpace()


# End of file 
