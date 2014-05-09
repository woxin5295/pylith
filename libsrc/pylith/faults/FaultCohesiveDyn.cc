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
// Copyright (c) 2010-2014 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "FaultCohesiveDyn.hh" // implementation of object methods

#include "CohesiveTopology.hh" // USES CohesiveTopology
#include "TractionPerturbation.hh" // HOLDSA TractionPerturbation

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Jacobian.hh" // USES Jacobian
#include "pylith/topology/SolutionFields.hh" // USES SolutionFields
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh
#include "pylith/topology/Stratum.hh" // USES Stratum, StratumIS

#include "pylith/friction/FrictionModel.hh" // USES FrictionModel

#include "pylith/feassemble/Quadrature.hh" // USES Quadrature
#include "pylith/feassemble/CellGeometry.hh" // USES CellGeometry

#include "pylith/utils/EventLogger.hh" // USES EventLogger
#include "pylith/utils/constdefs.h" // USES PYLITH_MAXFLOAT
#include "pylith/utils/macrodefs.h" // USES CALL_MEMBER_FN

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/spatialdb/SpatialDB.hh" // USES SpatialDB
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <cmath> // USES pow(), sqrt()
#include <strings.h> // USES strcasecmp()
#include <cstring> // USES strlen()
#include <cstdlib> // USES atoi()
#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

#include <iostream> // TEMPORARY

//#define DETAILED_EVENT_LOGGING

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::FaultCohesiveDyn::FaultCohesiveDyn(void) :
  _zeroTolerance(1.0e-10),
  _openFreeSurf(true),
  _tractionPerturbation(0),
  _friction(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::FaultCohesiveDyn::~FaultCohesiveDyn(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void pylith::faults::FaultCohesiveDyn::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  FaultCohesiveLagrange::deallocate();

  _tractionPerturbation = 0; // :TODO: Use shared pointer
  _friction = 0; // :TODO: Use shared pointer

  PYLITH_METHOD_END;
} // deallocate

// ----------------------------------------------------------------------
// Sets the spatial database for the inital tractions
void
pylith::faults::FaultCohesiveDyn::tractionPerturbation(TractionPerturbation* tract)
{ // tractionPerturbation
  _tractionPerturbation = tract;
} // tractionPerturbation

// ----------------------------------------------------------------------
// Get the friction (constitutive) model.  
void
pylith::faults::FaultCohesiveDyn::frictionModel(friction::FrictionModel* const model)
{ // frictionModel
  _friction = model;
} // frictionModel

// ----------------------------------------------------------------------
// Nondimensional tolerance for detecting near zero values.
void
pylith::faults::FaultCohesiveDyn::zeroTolerance(const PylithScalar value)
{ // zeroTolerance
  if (value < 0.0) {
    std::ostringstream msg;
    msg << "Tolerance (" << value << ") for detecting values near zero for "
      "fault " << label() << " must be nonnegative.";
    throw std::runtime_error(msg.str());
  } // if
  
  _zeroTolerance = value;
} // zeroTolerance

// ----------------------------------------------------------------------
// Set flag used to determine when fault is traction free when it
// opens or it still imposes any initial tractions.
void
pylith::faults::FaultCohesiveDyn::openFreeSurf(const bool value)
{ // openFreeSurf
  _openFreeSurf = value;
} // openFreeSurf

// ----------------------------------------------------------------------
// Initialize fault. Determine orientation and setup boundary
void
pylith::faults::FaultCohesiveDyn::initialize(const topology::Mesh& mesh,
					     const PylithScalar upDir[3])
{ // initialize
  PYLITH_METHOD_BEGIN;

  assert(upDir);
  assert(_quadrature);
  assert(_normalizer);

  FaultCohesiveLagrange::initialize(mesh, upDir);

  // Get initial tractions using a spatial database.
  if (_tractionPerturbation) {
    const topology::Field& orientation = _fields->get("orientation");
    _tractionPerturbation->initialize(*_faultMesh, orientation, *_normalizer);
  } // if

  // Setup fault constitutive model.
  assert(_friction);
  assert(_faultMesh);
  assert(_fields);
  _friction->normalizer(*_normalizer);
  _friction->initialize(*_faultMesh, _quadrature);

  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  assert(cs);

  topology::Field& dispRel = _fields->get("relative disp");

  // Create field for relative velocity associated with fault vertex
  _fields->add("relative velocity", "relative_velocity");
  topology::Field& velRel = _fields->get("relative velocity");
  velRel.cloneSection(dispRel);
  velRel.vectorFieldType(topology::FieldBase::VECTOR);
  velRel.scale(_normalizer->lengthScale() / _normalizer->timeScale());

  // Create field for contact traction associated with fault vertex
  _fields->add("contact traction", "traction");
  topology::Field& tractionContact = _fields->get("contact traction");
  tractionContact.cloneSection(dispRel);
  tractionContact.vectorFieldType(topology::FieldBase::VECTOR);
  tractionContact.scale(_normalizer->pressureScale());

  PYLITH_METHOD_END;
} // initialize

// ----------------------------------------------------------------------
// Integrate contributions to residual term (r) for operator.
void
pylith::faults::FaultCohesiveDyn::integrateResidual(const topology::Field& residual,
						    const PylithScalar t,
						    topology::SolutionFields* const fields)
{ // integrateResidual
  PYLITH_METHOD_BEGIN;

  assert(fields);
  assert(_fields);
  assert(_logger);

  // Cohesive cells with conventional vertices N and P, and constraint
  // vertex L make contributions to the assembled residual:
  //
  // \vec{T_c}_p = \vec{T_f}_p - \vec{l}_p 
  //
  // Locked (l_p > 0)
  // Sliding (l_p = 0)
  //
  // DOF P: \int_{S_f^+} \tensor{N}_m^T \cdot \tensor{N}_p \cdot \vec{T_c}_p dS
  // DOF N: -\int_{S_f^+} \tensor{N}_m^T \cdot \tensor{N}_p \cdot \vec{T_c}_p dS
  // DOF L: \int_S_f \tensor{R} \cdot \tensor{N}_p^T \cdot \vec{l}_p^{fault} \cdot 
  //                 \tensor{R} \cdot (-\tensor{N}_{n^+} \cdot \vec{u}_{n^+}
  //                                   +\tensor{N}_{n^-} \cdot \vec{u}_{n^-}) dS


  const int setupEvent = _logger->eventId("FaIR setup");
  const int geometryEvent = _logger->eventId("FaIR geometry");
  const int computeEvent = _logger->eventId("FaIR compute");
  const int restrictEvent = _logger->eventId("FaIR restrict");
  const int updateEvent = _logger->eventId("FaIR update");

  _logger->eventBegin(restrictEvent);
  _updateRelMotion(*fields);
  _logger->eventEnd(restrictEvent);

  _logger->eventBegin(setupEvent);

  // Get cell geometry information that doesn't depend on cell
  const int spaceDim = _quadrature->spaceDim();
  const int indexN = spaceDim - 1;

  // Get sections associated with cohesive cells
  PetscDM residualDM = residual.dmMesh();assert(residualDM);
  PetscSection residualGlobalSection = NULL;
  PetscErrorCode err = DMGetDefaultGlobalSection(residualDM, &residualGlobalSection);PYLITH_CHECK_ERROR(err);assert(residualGlobalSection);

  topology::VecVisitorMesh residualVisitor(residual);
  PetscScalar* residualArray = residualVisitor.localArray();

  topology::Field& dispT = fields->get("disp(t)");
  topology::VecVisitorMesh dispTVisitor(dispT);
  const PetscScalar* dispTArray = dispTVisitor.localArray();

  topology::Field& dispTIncr = fields->get("dispIncr(t->t+dt)");
  topology::VecVisitorMesh dispTIncrVisitor(dispTIncr);
  const PetscScalar* dispTIncrArray = dispTIncrVisitor.localArray();

  scalar_array tractionPerturbVertex(spaceDim);
  topology::VecVisitorMesh* tractionPerturbVisitor = 0;
  PetscScalar* tractionPerturbArray = NULL;
  if (_tractionPerturbation) {
    _tractionPerturbation->calculate(t);
    
    const topology::Fields* params = _tractionPerturbation->parameterFields();assert(params);
    const topology::Field& tractions = params->get("value");

    tractionPerturbVisitor = new topology::VecVisitorMesh(tractions);
    tractionPerturbArray = tractionPerturbVisitor->localArray();
  } // if
  tractionPerturbVertex = 0.0;

  scalar_array tractionInternalVertex(spaceDim);
  scalar_array tractionRheologyVertex(spaceDim);
  scalar_array tractionInternalGlobalVertex(spaceDim);
  scalar_array tractionRheologyGlobalVertex(spaceDim);
  scalar_array slipVertex(spaceDim);
  scalar_array slipRateVertex(spaceDim);

  topology::Field& orientation = _fields->get("orientation");
  topology::VecVisitorMesh orientationVisitor(orientation);
  const PetscScalar* orientationArray = orientationVisitor.localArray();

  topology::Field& area = _fields->get("area");
  topology::VecVisitorMesh areaVisitor(area);
  const PetscScalar* areaArray = areaVisitor.localArray();

  topology::Field& tractionContact = _fields->get("contact traction");
  topology::VecVisitorMesh tractionContactVisitor(tractionContact);
  PetscScalar* tractionContactArray = tractionContactVisitor.localArray();

  _logger->eventEnd(setupEvent);
#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  // Loop over fault vertices
  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    // Compute contribution only if Lagrange constraint is local.
    PetscInt goff = 0;
    err = PetscSectionGetOffset(residualGlobalSection, v_lagrange, &goff);PYLITH_CHECK_ERROR(err);
    if (goff < 0) {
      continue;
    } // if

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(restrictEvent);
#endif

    // Get prescribed traction perturbation at fault vertex.
    if (_tractionPerturbation) {
      const PetscInt toff = tractionPerturbVisitor->sectionOffset(v_fault);
      assert(spaceDim == tractionPerturbVisitor->sectionDof(v_fault));
      for (PetscInt iDim = 0; iDim < spaceDim; ++iDim) {
        tractionPerturbVertex[iDim] = tractionPerturbArray[toff+iDim];
      } // for
    } // if/else

    // Get orientation associated with fault vertex.
    const PetscInt ooff = orientationVisitor.sectionOffset(v_fault);
    assert(spaceDim*spaceDim == orientationVisitor.sectionDof(v_fault));

    // Get area associated with fault vertex.
    const PetscInt aoff = areaVisitor.sectionOffset(v_fault);
    assert(1 == areaVisitor.sectionDof(v_fault));

    // Get area associated with fault vertex.
    const PetscInt coff = tractionContactVisitor.sectionOffset(v_fault);
    assert(spaceDim == tractionContactVisitor.sectionDof(v_fault));

    // Get disp(t) at conventional vertices and Lagrange vertex.
    const PetscInt dtnoff = dispTVisitor.sectionOffset(v_negative);
    assert(spaceDim == dispTVisitor.sectionDof(v_negative));

    const PetscInt dtpoff = dispTVisitor.sectionOffset(v_positive);
    assert(spaceDim == dispTVisitor.sectionDof(v_positive));

    const PetscInt dtloff = dispTVisitor.sectionOffset(v_lagrange);
    assert(spaceDim == dispTVisitor.sectionDof(v_lagrange));

    // Get dispIncr(t->t+dt) at conventional vertices and Lagrange vertex.
    const PetscInt dinoff = dispTIncrVisitor.sectionOffset(v_negative);
    assert(spaceDim == dispTIncrVisitor.sectionDof(v_negative));

    const PetscInt dipoff = dispTIncrVisitor.sectionOffset(v_positive);
    assert(spaceDim == dispTIncrVisitor.sectionDof(v_positive));

    const PetscInt diloff = dispTIncrVisitor.sectionOffset(v_lagrange);
    assert(spaceDim == dispTIncrVisitor.sectionDof(v_lagrange));

    const PetscInt rnoff = residualVisitor.sectionOffset(v_negative);
    assert(spaceDim == residualVisitor.sectionDof(v_negative));

    const PetscInt rpoff = residualVisitor.sectionOffset(v_positive);
    assert(spaceDim == residualVisitor.sectionDof(v_positive));

    const PetscInt rloff = residualVisitor.sectionOffset(v_lagrange);
    assert(spaceDim == residualVisitor.sectionDof(v_lagrange));

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif
    
    // Compute slip, slip rate, and Lagrange multiplier at time t+dt
    // in fault coordinate system.
    slipVertex = 0.0;
    slipRateVertex = 0.0;
    tractionInternalVertex = 0.0;
    for (PetscInt iDim = 0; iDim < spaceDim; ++iDim) {
      tractionInternalVertex[iDim] = tractionPerturbVertex[iDim];
      
      for (PetscInt jDim = 0; jDim < spaceDim; ++jDim) {
        slipVertex[iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * (dispTArray[dtpoff+jDim] + dispTIncrArray[dipoff+jDim] - dispTArray[dtnoff+jDim] - dispTIncrArray[dinoff+jDim]);
        slipRateVertex[iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * (dispTIncrArray[dipoff+jDim] - dispTIncrArray[dinoff+jDim]) / _dt;
        tractionInternalVertex[iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * (dispTArray[dtloff+jDim] + dispTIncrArray[diloff+jDim]);
      } // for
      if (fabs(slipRateVertex[iDim]) < _zeroTolerance) {
        slipRateVertex[iDim] = 0.0;
      } // if
    } // for
    if (fabs(slipVertex[indexN]) < _zeroTolerance) {
      slipVertex[indexN] = 0.0;
    } // if
    
#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif
    if (slipVertex[indexN] < _zeroTolerance || !_openFreeSurf) { 
      // If no opening or flag indicates to still impose tractions
      // when fault is open, then assemble contributions into field
      
      const PylithScalar jacobianShearVertex = 0.0;
      bool needNewJacobianVertex = false;
      switch (spaceDim) { // switch
      case 2 :
	needNewJacobianVertex = _calcRheologyTraction2D(&tractionRheologyVertex, t, slipVertex, slipRateVertex, tractionInternalVertex, jacobianShearVertex);
	break;
      case 3 :
	needNewJacobianVertex = _calcRheologyTraction3D(&tractionRheologyVertex, t, slipVertex, slipRateVertex, tractionInternalVertex, jacobianShearVertex);
	break;
      } // switch
      if (needNewJacobianVertex) {
	_needNewJacobian = true;
      } // if

      tractionRheologyGlobalVertex = 0.0;
      tractionInternalGlobalVertex = 0.0;
      for (PetscInt iDim = 0; iDim < spaceDim; ++iDim) {
	tractionContactArray[coff+iDim] = tractionRheologyVertex[iDim] - tractionInternalVertex[iDim];

	for (PetscInt jDim = 0; jDim < spaceDim; ++jDim) {
	  tractionRheologyGlobalVertex += orientationArray[ooff+jDim*spaceDim+iDim] * tractionRheologyVertex[jDim];
	  tractionInternalGlobalVertex += orientationArray[ooff+jDim*spaceDim+iDim] * tractionInternalVertex[jDim];
	} // for
      } // for

      for (PetscInt iDim = 0; iDim < spaceDim; ++iDim) {
	const PylithScalar tractionTerm = areaArray[aoff] * (tractionRheologyGlobalVertex[iDim] - tractionInternalGlobalVertex[iDim]);
	residualArray[rnoff+iDim] += tractionTerm;
	residualArray[rpoff+iDim] -= tractionTerm;

	const PylithScalar dispRelVertex = dispTArray[dtpoff+iDim] + dispTIncrArray[dipoff+iDim] - dispTArray[dtnoff+iDim] - dispTIncrArray[dinoff+iDim];
	residualArray[rloff+iDim] += areaArray[aoff] * tractionInternalGlobalVertex[iDim] * dispRelVertex;
      } // for
    } else { // opening, normal traction should be zero
      for (PetscInt iDim = 0; iDim < spaceDim; ++iDim) {
	tractionContactArray[coff+iDim] = 0.0;
      } // for

      std::ostringstream msg;
      if (fabs(tractionInternalVertex[indexN]) > _zeroTolerance) {
        msg << "ERROR! Fault opening with nonzero traction."
                  << ", v_fault: " << v_fault
                  << ", opening: " << slipVertex[indexN]
                  << ", normal traction: " << tractionInternalVertex[indexN]
                  << std::endl;
        throw std::runtime_error(msg.str());
      } // if
    }  // if/else

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  PetscLogFlops(numVertices*spaceDim*8);
  delete tractionPerturbVisitor; tractionPerturbVisitor = 0;

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventEnd(computeEvent);
#endif

  PYLITH_METHOD_END;
} // integrateResidual

// ----------------------------------------------------------------------
// Update state variables as needed.
void
pylith::faults::FaultCohesiveDyn::updateStateVars(const PylithScalar t,
						  topology::SolutionFields* const fields)
{ // updateStateVars
  PYLITH_METHOD_BEGIN;

  assert(fields);
  assert(_fields);

  _updateRelMotion(*fields);

  const int spaceDim = _quadrature->spaceDim();

  // Allocate arrays for vertex values
  scalar_array tractionTpdtVertex(spaceDim); // Fault coordinate system

  // Get fields.
  topology::Field& dispT = fields->get("disp(t)");
  topology::VecVisitorMesh dispTVisitor(dispT);
  const PetscScalar* dispTArray = dispTVisitor.localArray();

  topology::Field& dispTIncr = fields->get("dispIncr(t->t+dt)");
  topology::VecVisitorMesh dispTIncrVisitor(dispTIncr);
  const PetscScalar* dispTIncrArray = dispTIncrVisitor.localArray();

  scalar_array slipVertex(spaceDim);
  topology::Field& dispRel = _fields->get("relative disp");
  topology::VecVisitorMesh dispRelVisitor(dispRel);
  const PetscScalar* dispRelArray = dispRelVisitor.localArray();

  scalar_array slipRateVertex(spaceDim);
  topology::Field& velRel = _fields->get("relative velocity");
  topology::VecVisitorMesh velRelVisitor(velRel);
  const PetscScalar* velRelArray = velRelVisitor.localArray();

  topology::Field& orientation = _fields->get("orientation");
  topology::VecVisitorMesh orientationVisitor(orientation);
  const PetscScalar* orientationArray = orientationVisitor.localArray();

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    // Get relative displacement
    const PetscInt droff = dispRelVisitor.sectionOffset(v_fault);
    assert(spaceDim == dispRelVisitor.sectionDof(v_fault));

    // Get relative velocity
    const PetscInt vroff = velRelVisitor.sectionOffset(v_fault);
    assert(spaceDim == velRelVisitor.sectionDof(v_fault));

    // Get orientation
    const PetscInt ooff = orientationVisitor.sectionOffset(v_fault);
    assert(spaceDim*spaceDim == orientationVisitor.sectionDof(v_fault));

    // Get Lagrange multiplier values from disp(t), and dispIncr(t->t+dt)
    const PetscInt dtloff = dispTVisitor.sectionOffset(v_lagrange);
    assert(spaceDim == dispTVisitor.sectionDof(v_lagrange));

    const PetscInt diloff = dispTIncrVisitor.sectionOffset(v_lagrange);
    assert(spaceDim == dispTIncrVisitor.sectionDof(v_lagrange));

    // Compute slip, slip rate, and fault traction (Lagrange
    // multiplier) at time t+dt in fault coordinate system.
    slipVertex = 0.0;
    slipRateVertex = 0.0;
    tractionTpdtVertex = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      for (int jDim=0; jDim < spaceDim; ++jDim) {
        slipVertex[iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * dispRelArray[droff+jDim];
        slipRateVertex[iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * velRelArray[vroff+jDim];
        tractionTpdtVertex[iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * (dispTArray[dtloff+jDim]+dispTIncrArray[diloff+jDim]);
      } // for
    } // for

    // Get friction properties and state variables.
    _friction->retrievePropsStateVars(v_fault);

    // Use fault constitutive model to compute traction associated with
    // friction.
    switch (spaceDim) { // switch
    case 2: { // case 2
      const PylithScalar slipMag = fabs(slipVertex[0]);
      const PylithScalar slipRateMag = fabs(slipRateVertex[0]);
      const PylithScalar tractionNormal = tractionTpdtVertex[1];
      _friction->updateStateVars(t, slipMag, slipRateMag, tractionNormal, v_fault);
      break;
    } // case 2
    case 3: { // case 3
      const PylithScalar slipMag = 
	sqrt(slipVertex[0]*slipVertex[0] + slipVertex[1]*slipVertex[1]);
      const PylithScalar slipRateMag = 
	sqrt(slipRateVertex[0]*slipRateVertex[0] + 
	     slipRateVertex[1]*slipRateVertex[1]);
      const PylithScalar tractionNormal = tractionTpdtVertex[2];
      _friction->updateStateVars(t, slipMag, slipRateMag, tractionNormal, v_fault);
      break;
    } // case 3
    default:
      assert(0);
      throw std::logic_error("Unknown spatial dimension in FaultCohesiveDyn::updateStateVars().");
    } // switch
  } // for

  PYLITH_METHOD_END;
} // updateStateVars

// ----------------------------------------------------------------------
// Adjust solution from solver with lumped Jacobian to match Lagrange
// multiplier constraints.
void
pylith::faults::FaultCohesiveDyn::adjustSolnLumped(topology::SolutionFields* const fields,
						   const PylithScalar t,
						   const topology::Field& jacobian)
{ // adjustSolnLumped
  PYLITH_METHOD_BEGIN;

  /// Member prototype for _constrainSolnSpaceXD()
  typedef void (pylith::faults::FaultCohesiveDyn::*constrainSolnSpace_fn_type)
    (scalar_array*,
     const PylithScalar,
     const scalar_array&,
     const scalar_array&,
     const scalar_array&,
     const PylithScalar,
     const bool);

  assert(fields);
  assert(_quadrature);

  // Cohesive cells with conventional vertices i and j, and constraint
  // vertex k require three adjustments to the solution:
  //
  //   * DOF k: Compute increment in Lagrange multipliers
  //            dl_k = S^{-1} (-C_ki (A_i^{-1} r_i - C_kj A_j^{-1} r_j + u_i - u_j) - d_k)
  //            S = C_ki (A_i^{-1} + A_j^{-1}) C_ki^T
  //
  //   * Adjust Lagrange multipliers to match friction criterion
  //
  //   * DOF k: Adjust displacement increment (solution) to create slip
  //     consistent with Lagrange multiplier constraints
  //            du_i = +A_i^-1 C_ki^T dlk
  //            du_j = -A_j^-1 C_kj^T dlk

  const int setupEvent = _logger->eventId("FaAS setup");
  const int geometryEvent = _logger->eventId("FaAS geometry");
  const int computeEvent = _logger->eventId("FaAS compute");
  const int restrictEvent = _logger->eventId("FaAS restrict");
  const int updateEvent = _logger->eventId("FaAS update");

  _logger->eventBegin(setupEvent);

  // Get cell information and setup storage for cell data
  const int spaceDim = _quadrature->spaceDim();

  // Allocate arrays for vertex values
  scalar_array tractionTpdtVertex(spaceDim);
  scalar_array lagrangeTpdtVertex(spaceDim);
  scalar_array dTractionTpdtVertex(spaceDim);
  scalar_array dLagrangeTpdtVertex(spaceDim);

  // Update time step in friction (can vary).
  _friction->timeStep(_dt);

  // Get section information
  scalar_array slipVertex(spaceDim);
  scalar_array dispRelVertex(spaceDim);
  topology::VecVisitorMesh dispRelVisitor(_fields->get("relative disp"));
  PetscScalar* dispRelArray = dispRelVisitor.localArray();

  scalar_array slipRateVertex(spaceDim);

  topology::VecVisitorMesh areaVisitor(_fields->get("area"));
  const PetscScalar* areaArray = areaVisitor.localArray();

  topology::VecVisitorMesh orientationVisitor(_fields->get("orientation"));
  const PetscScalar* orientationArray = orientationVisitor.localArray();

  topology::VecVisitorMesh dispTVisitor(fields->get("disp(t)"));
  const PetscScalar* dispTArray = dispTVisitor.localArray();

  scalar_array dispIncrVertexN(spaceDim);
  scalar_array dispIncrVertexP(spaceDim);
  scalar_array lagrangeTIncrVertex(spaceDim);
  topology::VecVisitorMesh dispTIncrVisitor(fields->get("dispIncr(t->t+dt)"));
  PetscScalar* dispTIncrArray = dispTIncrVisitor.localArray();

  topology::VecVisitorMesh dispTIncrAdjVisitor(fields->get("dispIncr adjust"));
  PetscScalar* dispTIncrAdjArray = dispTIncrAdjVisitor.localArray();

  topology::VecVisitorMesh jacobianVisitor(jacobian);
  const PetscScalar* jacobianArray = jacobianVisitor.localArray();

  topology::VecVisitorMesh residualVisitor(fields->get("residual"));
  const PetscScalar* residualArray = residualVisitor.localArray();

  PetscDM solnDM = fields->get("dispIncr(t->t+dt)").dmMesh();assert(solnDM);
  PetscSection solnGlobalSection = NULL;
  PetscErrorCode err = DMGetDefaultGlobalSection(solnDM, &solnGlobalSection);PYLITH_CHECK_ERROR(err);

  constrainSolnSpace_fn_type constrainSolnSpaceFn;
  switch (spaceDim) { // switch
  case 2: 
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDyn::_constrainSolnSpace2D;
    break;
  case 3:
    constrainSolnSpaceFn = 
      &pylith::faults::FaultCohesiveDyn::_constrainSolnSpace3D;
    break;
  default :
    assert(0);
    throw std::logic_error("Unknown spatial dimension in FaultCohesiveDyn::adjustSolnLumped.");
  } // switch

  _logger->eventEnd(setupEvent);

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventBegin(computeEvent);
#endif

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_lagrange = _cohesiveVertices[iVertex].lagrange;
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventBegin(restrictEvent);
#endif

    // Get residual at cohesive cell's vertices.
    const PetscInt rloff = residualVisitor.sectionOffset(v_lagrange);
    assert(spaceDim == residualVisitor.sectionDof(v_lagrange));

    // Get jacobian at cohesive cell's vertices.
    const PetscInt jnoff = jacobianVisitor.sectionOffset(v_negative);
    assert(spaceDim == jacobianVisitor.sectionDof(v_negative));

    const PetscInt jpoff = jacobianVisitor.sectionOffset(v_positive);
    assert(spaceDim == jacobianVisitor.sectionDof(v_positive));

    // Get disp(t) at Lagrange vertex.
    const PetscInt dtloff = dispTVisitor.sectionOffset(v_lagrange);
    assert(spaceDim == dispTVisitor.sectionDof(v_lagrange));

    // Get dispIncr(t) at cohesive cell's vertices.
    const PetscInt dinoff = dispTIncrVisitor.sectionOffset(v_negative);
    assert(spaceDim == dispTIncrVisitor.sectionDof(v_negative));

    const PetscInt dipoff = dispTIncrVisitor.sectionOffset(v_positive);
    assert(spaceDim == dispTIncrVisitor.sectionDof(v_positive));

    const PetscInt diloff = dispTIncrVisitor.sectionOffset(v_lagrange);
    assert(spaceDim == dispTIncrVisitor.sectionDof(v_lagrange));

    // Get relative displacement at fault vertex.
    const PetscInt droff = dispRelVisitor.sectionOffset(v_fault);
    assert(spaceDim == dispRelVisitor.sectionDof(v_fault));

    // Get area at fault vertex.
    const PetscInt aoff = areaVisitor.sectionOffset(v_fault);
    assert(1 == areaVisitor.sectionDof(v_fault));
    const PetscScalar areaVertex = areaArray[aoff];
    assert(areaVertex > 0.0);

    // Get fault orientation at fault vertex.
    const PetscInt ooff = orientationVisitor.sectionOffset(v_fault);
    assert(spaceDim*spaceDim == orientationVisitor.sectionDof(v_fault));

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(restrictEvent);
    _logger->eventBegin(computeEvent);
#endif

    // Adjust solution as in prescribed rupture, updating the Lagrange
    // multipliers and the corresponding displacment increments.
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      assert(jacobianArray[jpoff+iDim] > 0.0);
      assert(jacobianArray[jnoff+iDim] > 0.0);
      const PylithScalar S = (1.0/jacobianArray[jpoff+iDim] + 1.0/jacobianArray[jnoff+iDim]) * areaVertex*areaVertex;
      assert(S > 0.0);
      lagrangeTIncrVertex[iDim] = 1.0/S * (-residualArray[rloff+iDim] + areaVertex * (dispTIncrArray[dipoff+iDim] - dispTIncrArray[dinoff+iDim]));
      dispIncrVertexN[iDim] =  areaVertex / jacobianArray[jnoff+iDim]*lagrangeTIncrVertex[iDim];
      dispIncrVertexP[iDim] = -areaVertex / jacobianArray[jpoff+iDim]*lagrangeTIncrVertex[iDim];
    } // for

    // Compute slip, slip rate, and Lagrange multiplier at time t+dt
    // in fault coordinate system.
    slipVertex = 0.0;
    slipRateVertex = 0.0;
    tractionTpdtVertex = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      for (int jDim=0; jDim < spaceDim; ++jDim) {
        slipVertex[iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * dispRelArray[droff+jDim];
        tractionTpdtVertex[iDim] += orientationArray[ooff+iDim*spaceDim+jDim] * (dispTArray[dtloff+jDim] + lagrangeTIncrVertex[jDim]);
      } // for
    } // for
    // Jacobian is diagonal and isotropic, so it is invariant with
    // respect to rotation and contains one unique term.
    const PylithScalar jacobianShearVertex = -1.0 / (areaVertex * (1.0 / jacobianArray[jnoff+0] + 1.0 / jacobianArray[jpoff+0]));
    
    // Get friction properties and state variables.
    _friction->retrievePropsStateVars(v_fault);

    // Use fault constitutive model to compute traction associated with
    // friction.
    dTractionTpdtVertex = 0.0;

    const bool iterating = false; // No iteration for friction in lumped soln
    CALL_MEMBER_FN(*this, constrainSolnSpaceFn)(&dTractionTpdtVertex, t, slipVertex, slipRateVertex, tractionTpdtVertex, jacobianShearVertex, iterating);

    // Rotate traction back to global coordinate system.
    dLagrangeTpdtVertex = 0.0;
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      for (int jDim=0; jDim < spaceDim; ++jDim) {
        dLagrangeTpdtVertex[iDim] += orientationArray[ooff+jDim*spaceDim+iDim] * dTractionTpdtVertex[jDim];
      } // for
    } // for

#if 0 // debugging
    std::cout << "dispIncrP: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << dispIncrVertexP[iDim];
    std::cout << ", dispIncrN: "; 
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << dispIncrVertexN[iDim];
    std::cout << ", slipVertex: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << slipVertex[iDim];
    std::cout << ", slipRateVertex: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << slipRateVertex[iDim];
    std::cout << ", orientationVertex: ";
    for (int iDim=0; iDim < spaceDim*spaceDim; ++iDim)
      std::cout << "  " << orientationArray[ooff+iDim];
    std::cout << ", tractionVertex: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << tractionTpdtVertex[iDim];
    std::cout << ", lagrangeTVertex: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << lagrangeTVertex[iDim];
    std::cout << ", lagrangeTIncrVertex: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << lagrangeTIncrVertex[iDim];
    std::cout << ", dTractionTpdtVertex: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << dTractionTpdtVertex[iDim];
    std::cout << ", dLagrangeTpdtVertex: ";
    for (int iDim=0; iDim < spaceDim; ++iDim)
      std::cout << "  " << dLagrangeTpdtVertex[iDim];
    std::cout << std::endl;
#endif

    // Compute change in displacement.
    for (int iDim=0; iDim < spaceDim; ++iDim) {
      assert(jacobianArray[jpoff+iDim] > 0.0);
      assert(jacobianArray[jnoff+iDim] > 0.0);

      dispIncrVertexN[iDim] += areaVertex * dLagrangeTpdtVertex[iDim] / jacobianArray[jnoff+iDim];
      dispIncrVertexP[iDim] -= areaVertex * dLagrangeTpdtVertex[iDim] / jacobianArray[jpoff+iDim];

      // Update increment in Lagrange multiplier.
      lagrangeTIncrVertex[iDim] += dLagrangeTpdtVertex[iDim];
    } // for

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(computeEvent);
    _logger->eventBegin(updateEvent);
#endif

    // Compute contribution to adjusting solution only if Lagrange
    // constraint is local (the adjustment is assembled across processors).
    PetscInt goff;
    err = PetscSectionGetOffset(solnGlobalSection, v_lagrange, &goff);PYLITH_CHECK_ERROR(err);
    if (goff >= 0) {
      const PetscInt dianoff = dispTIncrAdjVisitor.sectionOffset(v_negative);
      assert(spaceDim == dispTIncrAdjVisitor.sectionDof(v_negative));

      const PetscInt diapoff = dispTIncrAdjVisitor.sectionOffset(v_positive);
      assert(spaceDim == dispTIncrAdjVisitor.sectionDof(v_positive));

      // Adjust displacements to account for Lagrange multiplier values
      // (assumed to be zero in preliminary solve).
      // Update displacement field
      for(PetscInt d = 0; d < spaceDim; ++d) {
        dispTIncrAdjArray[dianoff+d] += dispIncrVertexN[d];
        dispTIncrAdjArray[diapoff+d] += dispIncrVertexP[d];
      } // for
    } // if

    // The Lagrange multiplier and relative displacement are NOT
    // assembled across processors, so update even if Lagrange vertex
    // is not local.

    // Set Lagrange multiplier value. Value from preliminary solve is
    // bogus due to artificial diagonal entry in Jacobian of 1.0.
    for(PetscInt d = 0; d < spaceDim; ++d) {
      dispTIncrArray[diloff+d] = lagrangeTIncrVertex[d];
    } // for

#if defined(DETAILED_EVENT_LOGGING)
    _logger->eventEnd(updateEvent);
#endif
  } // for
  PetscLogFlops(numVertices*spaceDim*(17 + // adjust solve
                                      9 + // updates
                                      spaceDim*9));

#if !defined(DETAILED_EVENT_LOGGING)
  _logger->eventEnd(computeEvent);
#endif

#if 0 // DEBUGGING
  //dLagrangeTpdtSection->view("AFTER dLagrange");
  //dispIncrSection->view("AFTER DISP INCR (t->t+dt)");
#endif

  PYLITH_METHOD_END;
} // adjustSolnLumped

// ----------------------------------------------------------------------
// Get vertex field associated with integrator.
const pylith::topology::Field&
pylith::faults::FaultCohesiveDyn::vertexField(const char* name,
					      const topology::SolutionFields* fields)
{ // vertexField
  PYLITH_METHOD_BEGIN;

  assert(_faultMesh);
  assert(_quadrature);
  assert(_normalizer);
  assert(_fields);
  assert(_friction);

  const int cohesiveDim = _faultMesh->dimension();
  const int spaceDim = _quadrature->spaceDim();

  const topology::Field& orientation = _fields->get("orientation");

  PylithScalar scale = 0.0;
  int fiberDim = 0;
  if (0 == strcasecmp("slip", name)) {
    const topology::Field& dispRel = _fields->get("relative disp");
    _allocateBufferVectorField();
    topology::Field& buffer =  _fields->get("buffer (vector)");
    buffer.copy(dispRel);
    buffer.label("slip");
    FaultCohesiveLagrange::globalToFault(&buffer, orientation);
    PYLITH_METHOD_RETURN(buffer);

  } else if (0 == strcasecmp("slip_rate", name)) {
    const topology::Field& velRel = _fields->get("relative velocity");
    _allocateBufferVectorField();
    topology::Field& buffer = _fields->get("buffer (vector)");
    buffer.copy(velRel);
    buffer.label("slip_rate");
    FaultCohesiveLagrange::globalToFault(&buffer, orientation);
    PYLITH_METHOD_RETURN(buffer);

  } else if (cohesiveDim > 0 && 0 == strcasecmp("strike_dir", name)) {
    _allocateBufferVectorField();
    topology::Field& buffer = _fields->get("buffer (vector)");
    buffer.copySubfield(orientation, "strike_dir");
    PYLITH_METHOD_RETURN(buffer);

  } else if (2 == cohesiveDim && 0 == strcasecmp("dip_dir", name)) {
    _allocateBufferVectorField();
    topology::Field& buffer = _fields->get("buffer (vector)");
    buffer.copySubfield(orientation, "dip_dir");
    PYLITH_METHOD_RETURN(buffer);

  } else if (0 == strcasecmp("normal_dir", name)) {
    _allocateBufferVectorField();
    topology::Field& buffer = _fields->get("buffer (vector)");
    buffer.copySubfield(orientation, "normal_dir");
    PYLITH_METHOD_RETURN(buffer);

  } else if (0 == strcasecmp("traction", name)) {
    const topology::Field& tractionContact = _fields->get("contact traction");
    _allocateBufferVectorField();
    topology::Field& buffer = _fields->get("buffer (vector)");
    buffer.copy(tractionContact);
    buffer.label("traction");
    PYLITH_METHOD_RETURN(buffer);

  } else if (_friction->hasPropStateVar(name)) {
    PYLITH_METHOD_RETURN(_friction->getField(name));

  } else if (_tractionPerturbation && _tractionPerturbation->hasParameter(name)) {
    const topology::Field& param = _tractionPerturbation->vertexField(name, fields);
    if (param.vectorFieldType() == topology::FieldBase::VECTOR) {
      _allocateBufferVectorField();
      topology::Field& buffer = _fields->get("buffer (vector)");
      buffer.copy(param);
      FaultCohesiveLagrange::globalToFault(&buffer, orientation);
      PYLITH_METHOD_RETURN(buffer);
    } else {
      PYLITH_METHOD_RETURN(param);
    } // if/else

  } else {
    std::ostringstream msg;
    msg << "Request for unknown vertex field '" << name << "' for fault '"
        << label() << "'.";
    throw std::runtime_error(msg.str());
  } // else

  // Should never get here.
  throw std::logic_error("Unknown field in FaultCohesiveDyn::vertexField().");

  // Satisfy return values
  assert(_fields);
  const topology::Field& buffer = _fields->get("buffer (vector)");

  PYLITH_METHOD_RETURN(buffer);
} // vertexField

// ----------------------------------------------------------------------
// Update relative displacement and velocity (slip and slip rate)
// associated with Lagrange vertex k corresponding to diffential
// velocity between conventional vertices i and j.
void
pylith::faults::FaultCohesiveDyn::_updateRelMotion(const topology::SolutionFields& fields)
{ // _updateRelMotion
  PYLITH_METHOD_BEGIN;

  assert(_fields);

  const int spaceDim = _quadrature->spaceDim();

  // Get section information
  topology::VecVisitorMesh dispTVisitor(fields.get("disp(t)"));
  const PetscScalar* dispTArray = dispTVisitor.localArray();

  topology::VecVisitorMesh dispTIncrVisitor(fields.get("dispIncr(t->t+dt)"));
  const PetscScalar* dispTIncrArray = dispTIncrVisitor.localArray();

  topology::VecVisitorMesh velocityVisitor(fields.get("velocity(t)"));
  const PetscScalar* velocityArray = velocityVisitor.localArray();

  topology::VecVisitorMesh dispRelVisitor(_fields->get("relative disp"));
  PetscScalar* dispRelArray = dispRelVisitor.localArray();

  topology::VecVisitorMesh velRelVisitor(_fields->get("relative velocity"));
  PetscScalar* velRelArray = velRelVisitor.localArray();

  const int numVertices = _cohesiveVertices.size();
  for (int iVertex=0; iVertex < numVertices; ++iVertex) {
    const int v_fault = _cohesiveVertices[iVertex].fault;
    const int v_negative = _cohesiveVertices[iVertex].negative;
    const int v_positive = _cohesiveVertices[iVertex].positive;

    // Get displacement offsets.
    const PetscInt dtnoff = dispTVisitor.sectionOffset(v_negative);
    assert(spaceDim == dispTVisitor.sectionDof(v_negative));
    
    const PetscInt dtpoff = dispTVisitor.sectionOffset(v_positive);
    assert(spaceDim == dispTVisitor.sectionDof(v_positive));

    const PetscInt dinoff = dispTIncrVisitor.sectionOffset(v_negative);
    assert(spaceDim == dispTIncrVisitor.sectionDof(v_negative));

    const PetscInt dipoff = dispTIncrVisitor.sectionOffset(v_positive);
    assert(spaceDim == dispTIncrVisitor.sectionDof(v_positive));

    // Get velocity offsets.
    const PetscInt vnoff = velocityVisitor.sectionOffset(v_negative);
    assert(spaceDim == velocityVisitor.sectionDof(v_negative));

    const PetscInt vpoff = velocityVisitor.sectionOffset(v_positive);
    assert(spaceDim == velocityVisitor.sectionDof(v_positive));

    // Relative displacement/velocity offsets.
    const PetscInt droff = dispRelVisitor.sectionOffset(v_fault);
    assert(spaceDim == dispRelVisitor.sectionDof(v_fault));

    const PetscInt vroff = velRelVisitor.sectionOffset(v_fault);
    assert(spaceDim == velRelVisitor.sectionDof(v_fault));

    for(PetscInt d = 0; d < spaceDim; ++d) {
      const PylithScalar dispValue = dispTArray[dtpoff+d] + dispTIncrArray[dipoff+d] - dispTArray[dtnoff+d] - dispTIncrArray[dinoff+d];
      dispRelArray[droff+d] = fabs(dispValue) > _zeroTolerance ? dispValue : 0.0;

      const PylithScalar velValue = velocityArray[vpoff+d] - velocityArray[vnoff+d];
      velRelArray[vroff+d] = fabs(velValue) > _zeroTolerance ? velValue : 0.0;
    } // for

  } // for
  PetscLogFlops(numVertices*spaceDim*spaceDim*4);

  PYLITH_METHOD_END;
} // _updateRelMotion


// ----------------------------------------------------------------------
// Compute fault rheology traction for 2-D.
bool
pylith::faults::FaultCohesiveDyn::_calcRheologyTraction2D(scalar_array* tractionRheology,
							  const PylithScalar t,
							  const scalar_array& slip,
							  const scalar_array& slipRate,
							  const scalar_array& tractionInternal,
							  const PylithScalar jacobianShear)
{ // _calcRheologyTraction2D
  assert(tractionRheology);

  PylithScalar slipMag = fabs(slip[0]);
  const PylithScalar slipRateMag = fabs(slipRate[0]);

  const PylithScalar tractionNormal = tractionInternal[1];
  const PylithScalar tractionShear = tractionInternal[0];
  const PylithScalar tractionShearMag = fabs(tractionShear);

  bool needNewJacobian = false;

  if (fabs(slip[1]) < _zeroTolerance && tractionNormal < -_zeroTolerance) {
    // if in compression and no opening
    PylithScalar frictionStress = _friction->calcFriction(t, slipMag, slipRateMag, tractionNormal);

#if 1 // New Newton stuff
    if (tractionShearMag > 0.0 && 0.0 != jacobianShear) {
      assert(jacobianShear < 0.0);
      // Use Newton to get better update
      const int maxiter = 32;
      PylithScalar slipMagCur = slipMag;
      PylithScalar slipRateMagCur = slipRateMag;
      PylithScalar tractionShearMagCur = tractionShearMag;
      const PylithScalar slipMag0 = fabs(slip[0] - slipRate[0] * _dt);
      for (int iter=0; iter < maxiter; ++iter) {
	const PylithScalar frictionDeriv = _friction->calcFrictionDeriv(t, slipMagCur, slipRateMagCur, tractionNormal);
	slipMag = slipMagCur;
	if (slipMag > 0.0) {
	  // Use Newton (in log slip space) to get better update in slip & traction.
	  // D_{i+1} = exp(ln(D_i) - (T-T_f)/(D_i * (jacobian - frictionDeriv))
	  slipMagCur = exp(log(slipMag) - (tractionShearMagCur - frictionStress) / (slipMag * (jacobianShear - frictionDeriv)));
	} else {
	  // Use Newton (in linear slip space) to get better update in slip & traction.
	  // D_{i+1} = D_i - (T-T_f)/(jacobian - frictionDeriv)
	  slipMagCur = slipMag - (tractionShearMagCur - frictionStress) / (jacobianShear - frictionDeriv);
	} // if
	tractionShearMagCur += (slipMagCur - slipMag) * jacobianShear;
	slipRateMagCur = (slipMagCur - slipMag0) / _dt;
	frictionStress = _friction->calcFriction(t, slipMagCur, slipRateMagCur, tractionNormal);
	if (fabs(tractionShearMagCur - frictionStress) < _zeroTolerance) {
	  break;
	} // if
      } // for
    } // if
#endif

    // Set direction of rheology traction
    if (tractionShear >= 0.0) {
      (*tractionRheology)[0] = +frictionStress;
    } else {
      (*tractionRheology)[0] = -frictionStress;
    } // if/else
    
    // Determine if flipping between sliding and locked (means need new Jacobian)
    if (slipMag < _zeroTolerance && frictionStress > 0.0 && tractionShearMag < _zeroTolerance) {
      needNewJacobian = true;
    } else if (slipMag > _zeroTolerance && tractionShearMag > _zeroTolerance) {
      needNewJacobian = true;
    } // if/else

  } else {
    (*tractionRheology)[0] = 0.0;
    (*tractionRheology)[1] = 0.0;
  } // if/else

  PetscLogFlops(8);

  return needNewJacobian;
} // _calcRheologyTraction2D


// ----------------------------------------------------------------------
// Compute fault rheology traction for 3-D.
bool
pylith::faults::FaultCohesiveDyn::_calcRheologyTraction3D(scalar_array* tractionRheology,
							  const PylithScalar t,
							  const scalar_array& slip,
							  const scalar_array& slipRate,
							  const scalar_array& tractionInternal,
							  const PylithScalar jacobianShear)
{ // _calcRheologyTraction3D
  assert(tractionRheology);

  PylithScalar slipMag = sqrt(slip[0] * slip[0] + slip[1] * slip[1]);
  const PylithScalar slipRateMag = sqrt(slipRate[0]*slipRate[0] + slipRate[1]*slipRate[1]);
  
  const PylithScalar tractionNormal = tractionInternal[2];
  const PylithScalar tractionShearMag = sqrt(tractionInternal[0] * tractionInternal[0] + tractionInternal[1] * tractionInternal[1]);

  bool needNewJacobian = false;

  if (fabs(slip[2]) < _zeroTolerance && tractionNormal < -_zeroTolerance) {
    // if in compression and no opening
    PylithScalar frictionStress = _friction->calcFriction(t, slipMag, slipRateMag, tractionNormal);

#if 1 // New Newton stuff
    if (tractionShearMag > 0.0 && 0.0 != jacobianShear) {
      assert(jacobianShear < 0.0);
      // Use Newton to get better update
      const int maxiter = 32;
      PylithScalar slipMagCur = slipMag;
      PylithScalar slipRateMagCur = slipRateMag;
      PylithScalar tractionShearMagCur = tractionShearMag;
      const PylithScalar slipMag0 = sqrt(pow(slip[0]-slipRate[0]*_dt, 2) + pow(slip[1]-slipRate[1]*_dt, 2));
      for (int iter=0; iter < maxiter; ++iter) {
	const PylithScalar frictionDeriv = _friction->calcFrictionDeriv(t, slipMagCur, slipRateMagCur, tractionNormal);
	slipMag = slipMagCur;
	if (slipMag > 0.0) {
	  // Use Newton (in log slip space) to get better update in slip & traction.
	  // D_{i+1} = exp(ln(D_i) - (T-T_f)/(D_i * (jacobian - frictionDeriv))
	  slipMagCur = exp(log(slipMag) - (tractionShearMagCur - frictionStress) / (slipMag * (jacobianShear - frictionDeriv)));
	} else {
	  // Use Newton (in linear slip space) to get better update in slip & traction.
	  // D_{i+1} = D_i - (T-T_f)/(jacobian - frictionDeriv)
	  slipMagCur = slipMag - (tractionShearMagCur - frictionStress) / (jacobianShear - frictionDeriv);
	} // if
	tractionShearMagCur += (slipMagCur - slipMag) * jacobianShear;
	slipRateMagCur = (slipMagCur - slipMag0) / _dt;
	frictionStress = _friction->calcFriction(t, slipMagCur, slipRateMagCur, tractionNormal);
	if (fabs(tractionShearMagCur - frictionStress) < _zeroTolerance) {
	  break;
	} // if
      } // for
    } // if
#endif

    // Set direction of rheology traction
    if (tractionShearMag >= 0.0) {
      (*tractionRheology)[0] = +frictionStress * tractionInternal[0] / tractionShearMag;
      (*tractionRheology)[1] = +frictionStress * tractionInternal[1] / tractionShearMag;
    } else {
      (*tractionRheology)[0] = -frictionStress * tractionInternal[0] / tractionShearMag;
      (*tractionRheology)[1] = -frictionStress * tractionInternal[1] / tractionShearMag;
    } // if/else
    
    // Determine if flipping between sliding and locked (means need new Jacobian)
    if (slipMag < _zeroTolerance && frictionStress > 0.0 && tractionShearMag < _zeroTolerance) {
      needNewJacobian = true;
    } else if (slipMag > _zeroTolerance && tractionShearMag > _zeroTolerance) {
      needNewJacobian = true;
    } // if/else

  } else {
    (*tractionRheology)[0] = 0.0;
    (*tractionRheology)[1] = 0.0;
    (*tractionRheology)[2] = 0.0;
  } // if/else

  PetscLogFlops(8);

  return needNewJacobian;
} // _calcRheologyTraction3D


// ----------------------------------------------------------------------
// Constrain solution space in 2-D.
void
pylith::faults::FaultCohesiveDyn::_constrainSolnSpace2D(scalar_array* dTractionTpdt,
							const PylithScalar t,
							const scalar_array& slip,
							const scalar_array& slipRate,
							const scalar_array& tractionTpdt,
							const PylithScalar jacobianShear,
							const bool iterating)
{ // _constrainSolnSpace2D
  assert(dTractionTpdt);

  PylithScalar slipMag = fabs(slip[0]);
  const PylithScalar slipRateMag = fabs(slipRate[0]);

  const PylithScalar tractionNormal = tractionTpdt[1];
  const PylithScalar tractionShearMag = fabs(tractionTpdt[0]);
  
  if (fabs(slip[1]) < _zeroTolerance && tractionNormal < -_zeroTolerance) {
    // if in compression and no opening
    PylithScalar frictionStress = _friction->calcFriction(t, slipMag, slipRateMag, tractionNormal);

    if (tractionShearMag > frictionStress || (iterating && slipRateMag > 0.0)) {
      // traction is limited by friction, so have sliding OR
      // friction exceeds traction due to overshoot in slip

      if (tractionShearMag > 0.0) {
#if 1 // New Newton stuff
	if (0.0 != jacobianShear) {
	  assert(jacobianShear < 0.0);
	  // Use Newton to get better update
	  const int maxiter = 32;
	  PylithScalar slipMagCur = slipMag;
	  PylithScalar slipRateMagCur = slipRateMag;
	  PylithScalar tractionShearMagCur = tractionShearMag;
	  const PylithScalar slipMag0 = fabs(slip[0] - slipRate[0] * _dt);
	  for (int iter=0; iter < maxiter; ++iter) {
	    const PylithScalar frictionDeriv = _friction->calcFrictionDeriv(t, slipMagCur, slipRateMagCur, tractionNormal);
	    slipMag = slipMagCur;
	    if (slipMag > 0.0) {
	      // Use Newton (in log slip space) to get better update in slip & traction.
	      // D_{i+1} = exp(ln(D_i) - (T-T_f)/(D_i * (jacobian - frictionDeriv))
	      slipMagCur = exp(log(slipMag) - (tractionShearMagCur - frictionStress) / (slipMag * (jacobianShear - frictionDeriv)));
	    } else {
	      // Use Newton (in linear slip space) to get better update in slip & traction.
	      // D_{i+1} = D_i - (T-T_f)/(jacobian - frictionDeriv)
	      slipMagCur = slipMag - (tractionShearMagCur - frictionStress) / (jacobianShear - frictionDeriv);
	    } // if
	    tractionShearMagCur += (slipMagCur - slipMag) * jacobianShear;
	    slipRateMagCur = (slipMagCur - slipMag0) / _dt;
	    frictionStress = _friction->calcFriction(t, slipMagCur, slipRateMagCur, tractionNormal);
	    if (fabs(tractionShearMagCur - frictionStress) < _zeroTolerance) {
	      break;
	    } // if
	  } // for
	} // if
#endif

	// Update traction increment based on value required to stick
	// versus friction
	const PylithScalar dlp = -(tractionShearMag - frictionStress) * tractionTpdt[0] / tractionShearMag;
	(*dTractionTpdt)[0] = dlp;
      } else {
	// No shear stress and no friction.
      } // if/else
    } else {
      // friction exceeds value necessary to stick
      // no changes to solution
      if (iterating) {
	assert(0.0 == slipRateMag);
      } // if
    } // if/else
  } else {
    // if in tension, then traction is zero.
    (*dTractionTpdt)[0] = -tractionTpdt[0];
    (*dTractionTpdt)[1] = -tractionTpdt[1];
  } // else

  PetscLogFlops(8);
} // _constrainSolnSpace2D

// ----------------------------------------------------------------------
// Constrain solution space in 3-D.
void
pylith::faults::FaultCohesiveDyn::_constrainSolnSpace3D(scalar_array* dTractionTpdt,
							const PylithScalar t,
							const scalar_array& slip,
							const scalar_array& slipRate,
							const scalar_array& tractionTpdt,
							const PylithScalar jacobianShear,
							const bool iterating)
{ // _constrainSolnSpace3D
  assert(dTractionTpdt);

  PylithScalar slipMag = sqrt(slip[0] * slip[0] + slip[1] * slip[1]);
  const PylithScalar slipRateMag = sqrt(slipRate[0]*slipRate[0] + slipRate[1]*slipRate[1]);
  
  const PylithScalar tractionNormal = tractionTpdt[2];
  const PylithScalar tractionShearMag = sqrt(tractionTpdt[0] * tractionTpdt[0] + tractionTpdt[1] * tractionTpdt[1]);
  
  if (fabs(slip[2]) < _zeroTolerance && tractionNormal < -_zeroTolerance) {
    // if in compression and no opening
    PylithScalar frictionStress = _friction->calcFriction(t, slipMag, slipRateMag, tractionNormal);

    if (tractionShearMag > frictionStress || (iterating && slipRateMag > 0.0)) {
      // traction is limited by friction, so have sliding OR
      // friction exceeds traction due to overshoot in slip
      
      if (tractionShearMag > 0.0) {
#if 1 // New Newton stuff
	if (0.0 != jacobianShear) {
	  assert(jacobianShear < 0.0);
	  // Use Newton to get better update
	  const int maxiter = 32;
	  PylithScalar slipMagCur = slipMag;
	  PylithScalar slipRateMagCur = slipRateMag;
	  PylithScalar tractionShearMagCur = tractionShearMag;
	  const PylithScalar slipMag0 = sqrt(pow(slip[0]-slipRate[0]*_dt, 2) + pow(slip[1]-slipRate[1]*_dt, 2));
	  for (int iter=0; iter < maxiter; ++iter) {
	    const PylithScalar frictionDeriv = _friction->calcFrictionDeriv(t, slipMagCur, slipRateMagCur, tractionNormal);
	    slipMag = slipMagCur;
	    if (slipMag > 0.0) {
	      // Use Newton (in log slip space) to get better update in slip & traction.
	      // D_{i+1} = exp(ln(D_i) - (T-T_f)/(D_i * (jacobian - frictionDeriv))
	      slipMagCur = exp(log(slipMag) - (tractionShearMagCur - frictionStress) / (slipMag * (jacobianShear - frictionDeriv)));
	    } else {
	      // Use Newton (in linear slip space) to get better update in slip & traction.
	      // D_{i+1} = D_i - (T-T_f)/(jacobian - frictionDeriv)
	      slipMagCur = slipMag - (tractionShearMagCur - frictionStress) / (jacobianShear - frictionDeriv);
	    } // if
	    tractionShearMagCur += (slipMagCur - slipMag) * jacobianShear;
	    slipRateMagCur = (slipMagCur - slipMag0) / _dt;
	    frictionStress = _friction->calcFriction(t, slipMagCur, slipRateMagCur, tractionNormal);
	    if (fabs(tractionShearMagCur - frictionStress) < _zeroTolerance) {
	      break;
	    } // if
	  } // for
	} // if
#endif

	// Update traction increment based on value required to stick
	// versus friction
	const PylithScalar dlp = -(tractionShearMag - frictionStress) * tractionTpdt[0] / tractionShearMag;
	const PylithScalar dlq = -(tractionShearMag - frictionStress) * tractionTpdt[1] / tractionShearMag;
	
	(*dTractionTpdt)[0] = dlp;
	(*dTractionTpdt)[1] = dlq;
      } else {
	// No shear stress and no friction.
      } // if/else	
      
    } else {
      // else friction exceeds value necessary, so stick
      // no changes to solution
      if (iterating) {
	assert(0.0 == slipRateMag);
      } // if
    } // if/else
  } else {
    // if in tension, then traction is zero.
    (*dTractionTpdt)[0] = -tractionTpdt[0];
    (*dTractionTpdt)[1] = -tractionTpdt[1];
    (*dTractionTpdt)[2] = -tractionTpdt[2];
  } // else

  PetscLogFlops(22);
} // _constrainSolnSpace3D


// End of file 
