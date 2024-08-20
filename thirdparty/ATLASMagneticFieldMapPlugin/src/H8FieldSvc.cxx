/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

//
// H8FieldSvc.cxx
//
// Masahiro Morii, Harvard University
//

#include <fstream>
#include <iostream>

// ISF_Services include
#include "H8FieldSvc.h"

// PathResolver
// #include "PathResolver/PathResolver.h"

// StoreGate
// #include "StoreGate/StoreGateSvc.h"

// Athena Pool
// #include "AthenaPoolUtilities/CondAttrListCollection.h"

// CLHEP
#include "CLHEP/Units/SystemOfUnits.h"

/** Constructor **/
MagField::H8FieldSvc::H8FieldSvc(const std::string& /* name,ISvcLocator* svc */)
    :  // base_class(name,svc),
    m_H8MapFilename("MagneticFieldMaps/mbps1-all-id-800-mbps2-muons-800x4.data")
    , m_dx1(0)
    , m_dy1(0)
    , m_dz1(0)
    , m_dx2(43830)
    ,  // default H8MS2 x offset
    m_dy2(0)
    , m_dz2(0)
{
  //    declareProperty("H8MapFile", m_H8MapFilename, "File storing the H8
  //    magnetic field map"); declareProperty("dx1", m_dx1, "x component of
  //    magnet #1 shift (in mm)"); declareProperty("dy1", m_dy1, "y component of
  //    magnet #1 shift (in mm)"); declareProperty("dz1", m_dz1, "z component of
  //    magnet #1 shift (in mm)"); declareProperty("dx2", m_dx2, "x component of
  //    magnet #2 shift (in mm)"); declareProperty("dy2", m_dy2, "y component of
  //    magnet #2 shift (in mm)"); declareProperty("dz2", m_dz2, "z component of
  //    magnet #2 shift (in mm)");
}

MagField::H8FieldSvc::~H8FieldSvc() {}

/** framework methods */
bool MagField::H8FieldSvc::initialize()
{
  // ATH_MSG_INFO( "initialize() ..." );
  return true;
}

bool MagField::H8FieldSvc::start()
{
  // ATH_MSG_INFO( "start() ..." );
  return readMap(m_H8MapFilename);
}

bool MagField::H8FieldSvc::finalize()
{
  // ATH_MSG_INFO( "finalize() ..." );
  return true;
}

bool MagField::H8FieldSvc::readMap(const std::string mapFile)
{
  // find the path to the map file
  // std::string resolvedMapFile = PathResolver::find_file( mapFile.c_str(),
  // "DATAPATH" );
  std::string resolvedMapFile = mapFile.c_str();
  if (resolvedMapFile == "") {
    // ATH_MSG_ERROR( "Field map file " << mapFile << " not found" );
    return false;
  }
  // opne the map file
  std::ifstream input(resolvedMapFile.c_str());
  if (!input.good()) {
    // ATH_MSG_ERROR( "Failed to open the field map " << resolvedMapFile );
    return false;
  }
  // skip the file header line
  char line[256];
  input.getline(line, 256);
  // read grids
  unsigned igrid(0);
  double offset[3];
  while (1) {
    // first determine offset
    if (igrid == 0) {
      // magnet #1
      offset[0] = m_dx1;
      offset[1] = m_dy1;
      offset[2] = m_dz1;
    } else if (igrid == 1) {
      // magnet #2
      offset[0] = m_dx2;
      offset[1] = m_dy2;
      offset[2] = m_dz2;
    } else {
      // shift up to two magnets
      offset[0] = 0;
      offset[1] = 0;
      offset[2] = 0;
    }

    // then, read the map and shift it
    BFieldH8Grid grid;
    grid.readMap(input);
    grid.setOffset(offset);

    if (grid.defined()) {
      // ATH_MSG_INFO("setting offset for magnet " << igrid << " to " <<
      // offset[0] << ", " << offset[1] << ", " << offset[2] << " mm");
      //  save grid
      double its_min[3];
      double its_max[3];
      double its_d[3];
      grid.getBounds(its_min, its_max, its_d);
      // ATH_MSG_INFO("new magnet grid #" << igrid << " found");
      // ATH_MSG_INFO("  - min    (mm) " << its_min[0] << ", " << its_min[1] <<
      // ", " << its_min[2]); ATH_MSG_INFO("  - max    (mm) " << its_max[0] <<
      // ", " << its_max[1] << ", " << its_max[2]); ATH_MSG_INFO("  - offset
      // (mm) " << its_d[0] << ", " << its_d[1] << ", " << its_d[2]);
      m_grid.push_back(grid);
      igrid++;
    } else {
      break;
    }
  }
  // ATH_MSG_INFO( "Initialized the field map from " << resolvedMapFile );
  return true;
}

void MagField::H8FieldSvc::getField(const double* xyz,
                                    double* B,
                                    double* deriv) const
{
  for (unsigned i = 0; i < m_grid.size(); i++) {
    // find the grid that contains xyz
    if (m_grid[i].inside(xyz)) {
      m_grid[i].getB(xyz, B, deriv);
      return;
    }
  }
  // xyz is outside all grids
  B[0] = B[1] = B[2] = 0.0;
  if (deriv != 0) {
    for (int j = 0; j < 9; j++)
      deriv[j] = 0.0;
  }
  return;
}

void MagField::H8FieldSvc::getFieldZR(const double* xyz,
                                      double* B,
                                      double* deriv) const
{
  getField(xyz, B, deriv);
  return;
}
