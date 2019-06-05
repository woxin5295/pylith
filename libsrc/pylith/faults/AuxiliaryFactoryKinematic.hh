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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/faults/AuxiliaryFactoryKinematic.hh
 *
 * @brief C++ helper class for setting up auxiliary subfields for faults.
 */

#if !defined(pylith_faults_auxiliaryfactorykinematic_hh)
#define pylith_faults_auxiliaryfactorykinematic_hh

#include "faultsfwd.hh" // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

class pylith::faults::AuxiliaryFactoryKinematic : public pylith::feassemble::AuxiliaryFactory {
    friend class TestAuxiliaryFactoryKinematic; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactoryKinematic(void);

    /// Destructor.
    ~AuxiliaryFactoryKinematic(void);

    /// Add fault strike direction subfield to auxiliary field.
    void addStrikeDir(void);

    /// Add fault up-dip direction modulus subfield to auxiliary field.
    void addUpDipDir(void);

    /// Add fault normal direction subfield to auxiliary field.
    void addNormalDir(void);

    /// Add slip subfield to auxiliary field.
    void addSlip(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactoryKinematic(const AuxiliaryFactoryKinematic &); ///< Not implemented.
    const AuxiliaryFactoryKinematic& operator=(const AuxiliaryFactoryKinematic&); ///< Not implemented

}; // class AuxiliaryFactoryKinematic

#endif // pylith_faults_auxiliaryfactorykinematic_hh

// End of file