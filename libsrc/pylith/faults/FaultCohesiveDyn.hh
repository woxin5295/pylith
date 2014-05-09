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

/** @file libsrc/faults/FaultCohesiveDyn.hh
 *
 * @brief C++ implementation for a fault surface with spontaneous
 * (dynamic) slip implemented with cohesive cells.
 */

#if !defined(pylith_faults_faultcohesivedyn_hh)
#define pylith_faults_faultcohesivedyn_hh

// Include directives ---------------------------------------------------
#include "FaultCohesiveLagrange.hh" // ISA FaultCohesiveLagrange

#include "pylith/friction/frictionfwd.hh" // HOLDSA Friction model
#include "pylith/utils/petscfwd.h" // HASA PetscKSP

// FaultCohesiveDyn -----------------------------------------------------
/**
 * @brief C++ implementation for a fault surface with spontaneous
 * (dynamic) slip implemented with cohesive cells.
 *
 */
class pylith::faults::FaultCohesiveDyn : public FaultCohesiveLagrange
{ // class FaultCohesiveDyn
  friend class TestFaultCohesiveDyn; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  FaultCohesiveDyn(void);

  /// Destructor.
  virtual
  ~FaultCohesiveDyn(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Sets the traction perturbation for prescribed tractions.
   *
   * @param tract Spatial and temporal variation of tractions.
   */
  void tractionPerturbation(TractionPerturbation* tract);
  
  /** Set the friction (constitutive) model.
   *
   * @param model Fault constutive model.
   */
  void frictionModel(friction::FrictionModel* const model);

  /** Nondimensional tolerance for detecting near zero values.
   *
   * @param value Nondimensional tolerance
   */
  void zeroTolerance(const PylithScalar value);

  /** Set flag used to determine when fault is traction free when it
   * opens or it still imposes any initial tractions.
   *
   * If true, acts as a frictional contact. If false, one can simulate
   * a dike opening.
   *
   * @param value Nondimensional tolerance
   */
  void openFreeSurf(const bool value);

  /** Initialize fault. Determine orientation and setup boundary
   * condition parameters.
   *
   * @param mesh Finite-element mesh.
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction; applies to fault surfaces in 2-D and 3-D).
   */
  void initialize(const topology::Mesh& mesh,
		  const PylithScalar upDir[3]);

  /** Integrate contributions to residual term (r) for operator that
   * do not require assembly across processors.
   *
   * Initial tractions (if specified) are already assembled and
   * contribute to the residual like Neumann boundary conditions.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   */
  virtual
  void integrateResidual(const topology::Field& residual,
			 const PylithScalar t,
			 topology::SolutionFields* const fields);

  /** Update state variables as needed.
   *
   * @param t Current time
   * @param fields Solution fields
   * @param mesh Finite-element mesh
   */
  void updateStateVars(const PylithScalar t,
		       topology::SolutionFields* const fields);

  /** Adjust solution from solver with lumped Jacobian to match Lagrange
   *  multiplier constraints.
   *
   * @param fields Solution fields.
   * @param t Current time.
   * @param jacobian Jacobian of the system.
   */
  void adjustSolnLumped(topology::SolutionFields* fields,
			const PylithScalar t,
			const topology::Field& jacobian);

  /** Get vertex field associated with integrator.
   *
   * @param name Name of cell field.
   * @param fields Solution fields.
   * @returns Vertex field.
   */
  const topology::Field& vertexField(const char* name,
				     const topology::SolutionFields* fields =0);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Update relative displacement and velocity associated with
   * Lagrange vertex k corresponding to diffential velocity between
   * conventional vertices i and j.
   *
   * @param fields Solution fields.
   */
  void _updateRelMotion(const topology::SolutionFields& fields);

  /** Compute fault rheology traction for 2-D.
   *
   * @param tractionRheology Traction vector associated with fault rheology.
   * @param t Current time.
   * @param slip Fault slip vector.
   * @param slipRate Fault slip rate vector.
   * @param tractionInternal Traction vector associated with internal forces.
   * @param jacobianShear Derivative of shear traction with respect to slip (elasticity).
   *
   * @returns True if Jacobian needs recomputing, false otherwise.
   */
  bool _calcRheologyTraction2D(scalar_array* tractionRheology,
			       const PylithScalar t,
			       const scalar_array& slip,
			       const scalar_array& slipRate,
			       const scalar_array& tractionInternal,
			       const PylithScalar jacobianShear);

  /** Compute fault rheology traction for 3-D.
   *
   * @param tractionRheology Traction vector associated with fault rheology.
   * @param t Current time.
   * @param slip Fault slip vector.
   * @param slipRate Fault slip rate vector.
   * @param tractionInternal Traction vector associated with internal forces.
   * @param jacobianShear Derivative of shear traction with respect to slip (elasticity).
   *
   * @returns True if Jacobian needs recomputing, false otherwise.
   */
  bool _calcRheologyTraction3D(scalar_array* tractionRheology,
			       const PylithScalar t,
			       const scalar_array& slip,
			       const scalar_array& slipRate,
			       const scalar_array& tractionInternal,
			       const PylithScalar jacobianShear);

  /** Constrain solution space in 2-D.
   *
   * @param dLagrangeTpdt Adjustment to Lagrange multiplier.
   * @param t Current time.
   * @param slip Slip assoc. w/Lagrange multiplier vertex.
   * @param slipRate Slip rate assoc. w/Lagrange multiplier vertex.
   * @param tractionTpdt Fault traction assoc. w/Lagrange multiplier vertex.
   * @param jacobianShear Derivative of shear traction with respect to slip (elasticity).
   * @param iterating True if iterating on solution.
   */
  void _constrainSolnSpace2D(scalar_array* dLagrangeTpdt,
			     const PylithScalar t,
			     const scalar_array& slip,
			     const scalar_array& slipRate,
			     const scalar_array& tractionTpdt,
			     const PylithScalar jacobianShear,
			     const bool iterating =true);

  /** Constrain solution space in 3-D.
   *
   * @param dLagrangeTpdt Adjustment to Lagrange multiplier.
   * @param t Current time.
   * @param slip Slip assoc. w/Lagrange multiplier vertex.
   * @param slipRate Slip rate assoc. w/Lagrange multiplier vertex.
   * @param tractionTpdt Fault traction assoc. w/Lagrange multiplier vertex.
   * @param jacobianShear Derivative of shear traction with respect to slip (elasticity).
   * @param iterating True if iterating on solution.
   */
  void _constrainSolnSpace3D(scalar_array* dLagrangeTpdt,
			     const PylithScalar t,
			     const scalar_array& slip,
			     const scalar_array& slipRate,
			     const scalar_array& tractionTpdt,
			     const PylithScalar jacobianShear,
			     const bool iterating =true);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /// Flag to control whether to continue to impose initial tractions
  /// on the fault surface when it opens. If it is a frictional
  /// contact, then it should be a free surface.
  bool _openFreeSurf;

  /// Minimum resolvable value accounting for roundoff errors
  PylithScalar _zeroTolerance;

  /// Prescribed traction variation.
  TractionPerturbation* _tractionPerturbation;

  /// To identify constitutive model
  friction::FrictionModel* _friction;

// NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  FaultCohesiveDyn(const FaultCohesiveDyn&);

  /// Not implemented
  const FaultCohesiveDyn& operator=(const FaultCohesiveDyn&);

}; // class FaultCohesiveDyn

#endif // pylith_faults_faultcohesivedyn_hh


// End of file 
