// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2019 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/testing/MaterialStub.hh
 *
 * @brief Minimal C++ implementation of Material to allow testing of
 * basic Material functionality and use of Material objects in other
 * tests.
 */

#if !defined(pylith_materials_materialstub_hh)
#define pylith_materials_materialstub_hh

#include "pylith/testing/testingfwd.hh" // forward declarations

#include "pylith/testing/PhysicsStub.hh" // ISA PhysicsStub

class pylith::material::MaterialStub : public pylith::problems::PhysicsStub {

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    MaterialStub(void);

    /// Destructor
    ~MaterialStub(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    MaterialStub(const MaterialStub&); ///< Not implemented.
    const MaterialStub& operator=(const MaterialStub&); ///< Not implemented

}; // MaterialStub

#endif // pylith_materials_materialstub_hh


// End of file
