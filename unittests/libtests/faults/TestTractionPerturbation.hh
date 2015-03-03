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
 * @file unittests/libtests/faults/TestTractionPerturbation.hh
 *
 * @brief C++ TestTractionPerturbation object
 *
 * C++ unit testing for TractionPerturbation.
 */

#if !defined(pylith_faults_testtractionperturbation_hh)
#define pylith_faults_testtractionperturbation_hh

#include "pylith/faults/faultsfwd.hh" // USES TractionPerturbation
#include "pylith/topology/topologyfwd.hh" // USES Mesh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace faults {
    class TestTractionPerturbation;
  } // faults
} // pylith

/// C++ unit testing for TractionPerturbation
class pylith::faults::TestTractionPerturbation : public CppUnit::TestFixture
{ // class TestTractionPerturbation

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestTractionPerturbation );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testLabel );
  CPPUNIT_TEST( testHasParameter );
  CPPUNIT_TEST( testInitialize );
  CPPUNIT_TEST( testCalculate );
  CPPUNIT_TEST( testParameterFields );
  CPPUNIT_TEST( testVertexField );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test label().
  void testLabel(void);

  /// Test hasParameter().
  void testHasParameter(void);

  /// Test initialize() with 2-D mesh.
  void testInitialize(void);

  /// Test calculate() with 2-D mesh.
  void testCalculate(void);

  /// Test parameterFields() with 2-D mesh.
  void testParameterFields(void);

  /// Test VertexField() with 2-D mesh.
  void testVertexField(void);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Initialize TractionPerturbation.
   *
   * @param mesh Finite-element mesh of domain.
   * @param faultMesh Finite-element mesh of fault.
   * @param tract Traction perturbation.
   */
  static
  void _initialize(topology::Mesh* mesh,
		   topology::Mesh* faultMesh,
		   TractionPerturbation* tract);

}; // class TestTractionPerturbation

#endif // pylith_faults_testtractionperturbation_hh


// End of file 
