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

/** @file libsrc/materials/Material.hh
 *
 * @brief C++ abstract base class for materials.
 */

#if !defined(pylith_materials_material_hh)
#define pylith_materials_material_hh

#include "materialsfwd.hh" // forward declarations

#include "pylith/problems/Physics.hh" // ISA Physics

#include <string> // HASA std::string

// Material -------------------------------------------------------------
/** @brief C++ abstract base class for materials.
 *
 * Interface definition for a material. The material corresponds to the governing equation for a domain. The rheology is
 * a component of a material.
 *
 * An individual material must abide by specific rules for the interface, especially the order of the subfields in the
 * solution.
 */

class pylith::materials::Material : public virtual pylith::problems::Physics {
    friend class TestMaterial; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    Material(void);

    /// Destructor.
    virtual ~Material(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set value of label material-id used to identify material cells.
     *
     * @param value Material identifier
     */
    void setMaterialId(const int value);

    /** Get value of label material-id used to identify material cells.
     *
     * @returns Material identifier
     */
    int getMaterialId(void) const;

    /** Set descriptive label for material.
     *
     * @param value Label of material.
     */
    void setDescriptiveLabel(const char* value);

    /** Get descruptive label of material.
     *
     * @returns Label of material
     */
    const char* getDescriptiveLabel(void) const;

    /** Set gravity field.
     *
     * @param g Gravity field.
     */
    void setGravityField(spatialdata::spatialdb::GravityField* const g);

    /** Create constraint and set kernels.
     *
     * @param[in] solution Solution field.
     * @returns Constraint if applicable, otherwise NULL.
     */
    virtual
    pylith::feassemble::Constraint* createConstraint(const pylith::topology::Field& solution);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    spatialdata::spatialdb::GravityField* _gravityField; ///< Gravity field for gravitational body forces.

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    int _materialId; ///< Value of material-id label in mesh.
    std::string _descriptiveLabel; ///< Descriptive label for material.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    Material(const Material&); ///< Not implemented.
    const Material& operator=(const Material&); ///< Not implemented

}; // Material

#endif // pylith_materials_material_hh

// End of file
