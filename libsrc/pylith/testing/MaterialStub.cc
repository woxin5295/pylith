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
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "MaterialStub.hh" // Implementation of class methods

// ---------------------------------------------------------------------------------------------------------------------
// Constructor.
pylith::materials::MaterialStub::MaterialStub(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor
pylith::materials::MaterialStub::~MaterialStub(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::materials::MaterialStub::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    pylith::problems::PhysicsStub::deallocate();
    pylith::materials::Material::deallocate();

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Create constraint and set kernels.
pylith::feassemble::Constraint*
pylith::materials::MaterialStub::createConstraint(const pylith::topology::Field& solution) {
    return pylith::materials::Material::createConstraint(solution);
} // createConstraint


// End of file
