#!/usr/bin/env python
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2014 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

## @file unittests/pytests/faults/TestTractionPerturbation.py

## @brief Unit testing of TractionPerturbation object.

import unittest

from pylith.faults.TractionPerturbation import TractionPerturbation

# ----------------------------------------------------------------------
class TestTractionPerturbation(unittest.TestCase):
  """
  Unit testing of TractionPerturbation object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    tract = TractionPerturbation()
    return


  def test_configure(self):
    """
    Test initialize().
    """
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    from pyre.units.time import second

    ioInitial = SimpleIOAscii()
    ioInitial.inventory.filename = "tri3_initialtractions.spatialdb"
    ioInitial._configure()
    dbInitial = SimpleDB()
    dbInitial.inventory.iohandler = ioInitial
    dbInitial.inventory.label = "initial tractions"
    dbInitial._configure()
    
    ioChange = SimpleIOAscii()
    ioChange.inventory.filename = "tri3_changetractions.spatialdb"
    ioChange._configure()
    dbChange = SimpleDB()
    dbChange.inventory.iohandler = ioChange
    dbChange.inventory.label = "traction change"
    dbChange._configure()
    
    tract = TractionPerturbation()
    tract.inventory.dbInitial = dbInitial
    tract.inventory.dbChange = dbChange
    tract._configure()
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.faults.TractionPerturbation import traction_perturbation
    fn = traction_perturbation()
    return


# End of file 
