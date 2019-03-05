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
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/bc/FaultCohesiveStub.hh
 *
 * @brief Minimal C++ implementation of FaultCohesive to allow unit tests of
 * other objects needing meshes with cohesive cells.
 */

#if !defined(pylith_faults_faultcohesivestub_hh)
#define pylith_faults_faultcohesivestub_hh

#include "pylith/faults/FaultCohesive.hh" // ISA FaultCohesive

class pylith::faults::FaultCohesiveStub : public pylith::faults::FaultCohesive {
    friend class TestFaultCohesiveStub; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    FaultCohesiveStub(void);

    /// Destructor.
    ~FaultCohesiveStub(void);

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    void verifyConfiguration(const pylith::topology::Field& solution) const;

    /** Create integrator and set kernels.
     *
     * @param[in] solution Solution field.
     * @returns Integrator if applicable, otherwise NULL.
     */
    pylith::feassemble::Integrator* createIntegrator(const pylith::topology::Field& solution);

    /** Create constraint and set kernels.
     *
     * @param[in] solution Solution field.
     * @returns Constraint if applicable, otherwise NULL.
     */
    pylith::feassemble::Constraint* createConstraint(const pylith::topology::Field& solution);

    /** Create auxiliary field.
     *
     * @param[in] solution Solution field.
     * @param[in\ physicsMesh Finite-element mesh associated with physics.
     *
     * @returns Auxiliary field if applicable, otherwise NULL.
     */
    pylith::topology::Field* createAuxiliaryField(const pylith::topology::Field& solution,
                                                  const pylith::topology::Mesh& physicsMesh);

    /** Create derived field.
     *
     * @param[in] solution Solution field.
     * @param[in\ physicsMesh Finite-element mesh associated with physics.
     *
     * @returns Derived field if applicable, otherwise NULL.
     */
    pylith::topology::Field* createDerivedField(const pylith::topology::Field& solution,
                                                const pylith::topology::Mesh& physicsMesh);

    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    pylith::feassemble::AuxiliaryFactory* _getAuxiliaryFactory(void);

    /** Update kernel constants.
     *
     * @param[in] dt Current time step.
     */
    void _updateKernelConstants(const PylithReal dt);

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    FaultCohesiveStub(const FaultCohesiveStub&); ///< Not implemented.
    const FaultCohesiveStub& operator=(const FaultCohesiveStub&); ///< Not implemented.

}; // class FaultCohesiveStub

class pylith::faults::FaultCohesiveStubException {
    // PUBLIC ENUMS ////////////////////////////////////////////////////////////////////////////////////////////////////
public:

    enum MethodEnum {
        VERIFY_CONFIGURATION=0,
        CREATE_INTEGRATOR=1,
        CREATE_CONSTRAINT=2,
        CREATE_AUXILIARY_FIELD=3,
        CREATE_DERIVED_FIELD=4,
    }; // MethodEnum

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** Constructor.
     *
     * @param[in] value Method called.
     */
    FaultCohesiveStubException(const MethodEnum value);

    /// Destructor
    ~FaultCohesiveStubException(void);

    /** Get method called.
     *
     * @returns Method called.
     */
    MethodEnum getMethodCalled(void) const;

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
private:

    MethodEnum _methodCalled;
}; // FaultCohesiveStubException

#endif // pylith_faults_faultcohesivestub_hh

// End of file
