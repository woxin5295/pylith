// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestPointForceTri3.hh" // Implementation of class methods

#include "data/PointForceDataTri3.hh" // USES DirichletDataTri3

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::bc::TestPointForceTri3 );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::bc::TestPointForceTri3::setUp(void)
{ // setUp
  _data = new PointForceDataTri3();
} // setUp


// End of file 