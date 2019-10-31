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
 * @file unittests/libtests/materials/TestElasticity.hh
 *
 * @brief C++ class for testing material objects.
 */

#if !defined(pylith_materials_testmaterial_hh)
#define pylith_materials_testmaterial_hh

#include <cppunit/extensions/HelperMacros.h>
#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/materials/materialsfwd.hh" // forward declarations
#include "pylith/topology/topologyfwd.hh" // forward declarations

/// Namespace for pylith package
namespace pylith {
    namespace materials {
        class TestElasticity;
    } // materials
} // pylith

class pylith::materials::TestElasticity : public CppUnit::TestFixture, public pylith::utils::GenericComponent {
    // CPPUNIT TEST SUITE /////////////////////////////////////////////////
    CPPUNIT_TEST_SUITE(TestElasticity);

    CPPUNIT_TEST(testAccessors);
    CPPUNIT_TEST(testVerifyConfiguration);
    CPPUNIT_TEST(testCreateIntegrator);
    CPPUNIT_TEST(testCreateAuxiliaryField);
    CPPUNIT_TEST(testCreateDerivedField);

    CPPUNIT_TEST_SUITE_END();

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Test getters/setters.
    void testAccessors(void);

    /// Test verifyConfiguration().
    void testVerifyConfiguration(void);

    /// Test createIntegrator().
    void testCreateIntegrator(void);

    /// Test createAuxiliaryField().
    void testCreateAuxiliaryField(void);

    /// Test createDerivedField().
    void testCreateDerivedField(void);

}; // class TestElasticity

#endif // pylith_materials_testmaterial_hh

// End of file
