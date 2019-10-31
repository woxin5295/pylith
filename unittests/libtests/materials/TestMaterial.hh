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

/**
 * @file unittests/libtests/materials/TestMaterial.hh
 *
 * @brief C++ abstract base class for testing material objects.
 */

#if !defined(pylith_materials_testmaterial_hh)
#define pylith_materials_testmaterial_hh

#include <cppunit/extensions/HelperMacros.h>
#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/materials/materialsfwd.hh" // forward declarations
#include "pylith/topology/topologyfwd.hh" // forward declarations
#include "pylith/topology/Field.hh" // HASA FieldBase::Discretization

#include "spatialdata/spatialdb/spatialdbfwd.hh" // HOLDSA UserFunctionDB
#include "spatialdata/geocoords/geocoordsfwd.hh" // HOLDSA CoordSys
#include "spatialdata/units/unitsfwd.hh" // HOLDSA Nondimensional

/// Namespace for pylith package
namespace pylith {
    namespace materials {
        class TestMaterial;
    } // materials
} // pylith

/// C++ abstract base class for testing material objects.
class pylith::materials::TestMaterial : public CppUnit::TestFixture, public pylith::utils::GenericComponent {
    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestMaterial);

    CPPUNIT_TEST(testAccessors);
    CPPUNIT_TEST(testCreateConstraint);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Test set/getMaterialId(), set/getDescriptiveLabel(), setGravityField().
    void testAccessors(void);

    /// Test createConstraint().
    void testCreateConstraint(void);

}; // class TestMaterial

#endif // pylith_materials_testmaterial_hh

// End of file
