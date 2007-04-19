// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

/* Original mesh
 *
 *              1
 *             /|\
 *            / | \
 *           /  |  \
 *          /   |   \
 *         0    |    3
 *          \   |   /
 *           \  |  /
 *            \ | /
 *             \|/
 *              2
 *
 *
 * After adding cohesive elements
 *
 *              1 -- 4
 *             /|    |\
 *            / |    | \
 *           /  |    |  \
 *          /   |    |   \
 *         0    |    |    3
 *          \   |    |   /
 *           \  |    |  /
 *            \ |    | /
 *             \|    |/
 *              2 -- 5
 */

#include "CohesiveDataTri3.hh"

const int pylith::faults::CohesiveDataTri3::_numVertices = 6;

const int pylith::faults::CohesiveDataTri3::_spaceDim = 2;

const int pylith::faults::CohesiveDataTri3::_numCells = 3;

const int pylith::faults::CohesiveDataTri3::_cellDim = 2;

const double pylith::faults::CohesiveDataTri3::_vertices[] = {
 -1.0,  0.0,
  0.0,  1.0,
  0.0, -1.0,
  1.0,  0.0,
  0.0,  1.0,
  0.0, -1.0
};

const int pylith::faults::CohesiveDataTri3::_numCorners[] = {
  3,
  3,
  4
};

const int pylith::faults::CohesiveDataTri3::_cells[] = {
  0,  2,  1,
  4,  5,  3,
  1,  2,  4, 5
};

const int pylith::faults::CohesiveDataTri3::_materialIds[] = {
  0,  0,
  1
};

const int pylith::faults::CohesiveDataTri3::_numGroups = 2;

const int pylith::faults::CohesiveDataTri3::_groupSizes[] = 
  { 2, 3 };

const int pylith::faults::CohesiveDataTri3::_groups[] = {
  1, 2,
  1, 2, 3
};

const char* pylith::faults::CohesiveDataTri3::_groupNames[] = {
  "fault", "output"
};

const char* pylith::faults::CohesiveDataTri3::_groupTypes[] = {
  "vertex", "vertex"
};

const char* pylith::faults::CohesiveDataTri3::_filename = "data/meshTri3A.txt";

pylith::faults::CohesiveDataTri3::CohesiveDataTri3(void)
{ // constructor
  numVertices = _numVertices;
  spaceDim = _spaceDim;
  numCells = _numCells;
  cellDim = _cellDim;
  vertices = const_cast<double*>(_vertices);
  numCorners = const_cast<int*>(_numCorners);
  cells = const_cast<int*>(_cells);
  materialIds = const_cast<int*>(_materialIds);
  groups = const_cast<int*>(_groups);
  groupSizes = const_cast<int*>(_groupSizes);
  groupNames = const_cast<char**>(_groupNames);
  groupTypes = const_cast<char**>(_groupTypes);
  numGroups = _numGroups;
  filename = const_cast<char*>(_filename);
} // constructor

pylith::faults::CohesiveDataTri3::~CohesiveDataTri3(void)
{}


// End of file
