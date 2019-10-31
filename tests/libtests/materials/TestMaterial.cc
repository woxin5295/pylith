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

#include "TestMaterial.hh" // Implementation of class methods

#include "pylith/testing/MaterialStub.hh" // USES MaterialStub

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/spatialdb/GravityField.hh" // USES GravityField

#include "journal/debug.h" // USES journal::debug_t

CPPUNIT_TEST_SUITE_REGISTRATION(pylith::materials::TestMaterial);

// ---------------------------------------------------------------------------------------------------------------------
// Test set/getMaterialId(), set/getDescriptiveLabel(), setGravityField().
void
pylith::materials::TestMaterial::testAccessors(void) {
    PYLITH_METHOD_BEGIN;

    pylith::materials::MaterialStub material;

    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in default material id.", 0, material.getMaterialId());
    const int matId = 20;
    material.setMaterialId(matId);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in material id.", matId, material.getMaterialId());

    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in default material label.", std::string(""), std::string(material.getDescriptiveLabel()));
    const std::string& label = "my material";
    material.setDescriptiveLabel(label.c_str());
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in material label.", label, std::string(material.getDescriptiveLabel()));

    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in default gravity field.", (spatialdata::spatialdb::GravityField*)NULL, material._gravityField);
    spatialdata::spatialdb::GravityField gravity;
    material.setGravityField(&gravity);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Mismatch in gravity field.", &gravity, material._gravityField);

    PYLITH_METHOD_END;
} // testAccessors


// ---------------------------------------------------------------------------------------------------------------------
// Test createConstraint().
void
pylith::materials::TestMaterial::testCreateConstraint(void) {
    PYLITH_METHOD_BEGIN;

    pylith::topology::Mesh mesh;
    pylith::topology::Field solution(mesh);

    pylith::materials::MaterialStub material;
    CPPUNIT_ASSERT_EQUAL((pylith::feassemble::Constraint*)NULL, material.createConstraint(solution));

    PYLITH_METHOD_END;
} // testCreateConstraint


// End of file
