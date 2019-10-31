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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestElasticity.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

#include "journal/debug.h" // USES journal::debug_t

// ---------------------------------------------------------------------------------------------------------------------
// Test setters/getters.
void
pylith::materials::TestElasticity::testAccessors(void) {
    PYLITH_METHOD_BEGIN;

    pylith::materials::Elasticity material;

    CPPUNIT_ASSERT_EQUAL_MESSAGE("Default inertia flag", false, material.useInertia());
    material.useInertia(true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Inertia flag", true, material.useInertia());

    CPPUNIT_ASSERT_EQUAL_MESSAGE("Default body force flag", false, material.useBodyForce());
    material.useBodyForce(true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Body force flag", true, material.useBodyForce());

    CPPUNIT_ASSERT_MESSAGE(":TODO: @brad Finish implementing test.")

    PYLITH_METHOD_END;
} // testAccessors


// End of file
