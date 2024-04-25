/*
  Copyright (C) 2002-2024 CERN for the benefit of the ATLAS collaboration
*/

// ***************************************************************************
// Liquid Argon FCAL detector description package
// -----------------------------------------
//
//
// 10-Sep-2000 S.Simion   Handling of the FCAL read-out identifiers
//    Jan-2001 R.Sobie    Modify for persistency
//    Feb-2002 R.Sobie    Use same FCAL geometry files as simulation
//****************************************************************************

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "FastCaloSim/Geometry/FCAL_ChannelMap.h"

#include <stdio.h>

/* === Geometrical parameters === */
// const double cm = 0.01;
const double cm = 10.;
// const double FCAL_ChannelMap::m_tubeSpacing[] = {0.75*CLHEP::cm,
// 0.8179*CLHEP::cm, 0.90*CLHEP::cm};
const double FCAL_ChannelMap::m_tubeSpacing[] = {
    0.75 * cm, 0.8179 * cm, 0.90 * cm};

FCAL_ChannelMap::FCAL_ChannelMap(int flag)
{
  static const double root3_2 = std::sqrt(3.) / 2;
  /* === Initialize geometrical dimensions */
  for (int i = 0; i < 3; i++) {
    m_tubeDx[i] = m_tubeSpacing[i] / 2.;
    m_tubeDy[i] = m_tubeSpacing[i] * root3_2;
  }

  // FCAL1 small cells are 2x2 tubes
  m_tileDx[0] = 2. * m_tubeSpacing[0];
  m_tileDy[0] = 2. * m_tubeSpacing[0] * root3_2;

  // FCAL2 small cells are 2x3 tubes
  m_tileDx[1] = 2. * m_tubeSpacing[1];
  m_tileDy[1] = 3. * m_tubeSpacing[1] * root3_2;

  // FCAL3 cells are 6x6 tubes
  m_tileDx[2] = 6. * m_tubeSpacing[2];
  m_tileDy[2] = 6. * m_tubeSpacing[2] * root3_2;

  m_invert_x = flag & 1;
  m_invert_xy = flag & 2;
}

void FCAL_ChannelMap::finish()
{
  create_tileMap(1);
  create_tileMap(2);
  create_tileMap(3);
}

// *********************************************************************
// Read tube mapping tables
//
// Jan 23,2002    R. Sobie
// ********************************************************************

// original
void FCAL_ChannelMap::add_tube(const std::string& tileName,
                               int mod,
                               int /*id*/,
                               int i,
                               int j,
                               double x,
                               double y)
{
  // Get three integers from the tileName:
  std::istringstream tileStream1(std::string(tileName, 1, 1));
  std::istringstream tileStream2(std::string(tileName, 3, 2));
  std::istringstream tileStream3(std::string(tileName, 6, 3));
  int a1 = 0, a2 = 0, a3 = 0;
  if (tileStream1)
    tileStream1 >> a1;
  if (tileStream2)
    tileStream2 >> a2;
  if (tileStream3)
    tileStream3 >> a3;

  tileName_t tilename = (a3 << 16) + a2;

  TubePosition tb(tilename, x * cm, y * cm, "");
  // Add offsets, because iy and ix can be negative HMA

  i = i + 200;
  j = j + 200;
  //  m_tubeMap[mod-1][(j <<  16) + i] = tb;
  unsigned int ThisId = (j << 16) + i;
  tubemap_const_iterator p = m_tubeMap[mod - 1].insert(
      m_tubeMap[mod - 1].end(), std::make_pair(ThisId, tb));
  m_tubeIndex[mod - 1].push_back(p);
}

// Gabe: new to include HV and LARFCALELECRTODES ID
void FCAL_ChannelMap::add_tube(const std::string& tileName,
                               int mod,
                               int /*id*/,
                               int i,
                               int j,
                               double x,
                               double y,
                               const std::string& hvFT)
{
  // Get three integers from the tileName:
  std::istringstream tileStream1(std::string(tileName, 1, 1));
  std::istringstream tileStream2(std::string(tileName, 3, 2));
  std::istringstream tileStream3(std::string(tileName, 6, 3));
  int a1 = 0, a2 = 0, a3 = 0;
  if (tileStream1)
    tileStream1 >> a1;
  if (tileStream2)
    tileStream2 >> a2;
  if (tileStream3)
    tileStream3 >> a3;

  tileName_t tilename = (a3 << 16) + a2;

  TubePosition tb(tilename, x * cm, y * cm, hvFT);
  // Add offsets, because iy and ix can be negative HMA

  i = i + 200;
  j = j + 200;
  //  m_tubeMap[mod-1][(j <<  16) + i] = tb;
  unsigned int ThisId = (j << 16) + i;
  tubemap_const_iterator p = m_tubeMap[mod - 1].insert(
      m_tubeMap[mod - 1].end(), std::make_pair(ThisId, tb));
  m_tubeIndex[mod - 1].push_back(p);
}

FCAL_ChannelMap::tubemap_const_iterator FCAL_ChannelMap::getTubeByCopyNumber(
    int isam, int copyNo) const
{
  return m_tubeIndex[isam - 1][copyNo];
}

// *********************************************************************
// Create tile mapping tables
//
// Jan 23,2002    R. Sobie
// ********************************************************************
void FCAL_ChannelMap::create_tileMap(int isam)
{
  tileMap_const_iterator tile;
  tubemap_const_iterator first = m_tubeMap[isam - 1].begin();
  tubemap_const_iterator last = m_tubeMap[isam - 1].end();

  // Loop over tubes -> find unique tiles and fill the descriptors
  while (first != last) {
    tileName_t tileName = (first->second).get_tileName();
    tile = m_tileMap[isam - 1].find(tileName);

    if (tile == m_tileMap[isam - 1].end()) {  // New tile found
      float x = (first->second).x();
      float y = (first->second).y();
      unsigned int ntubes = 1;
      TilePosition tp(x, y, ntubes);
      m_tileMap[isam - 1][tileName] = tp;
    } else {  // Existing tile
      float x = (tile->second).x() + (first->second).x();
      float y = (tile->second).y() + (first->second).y();
      unsigned int ntubes = (tile->second).ntubes() + 1;
      TilePosition tp(x, y, ntubes);
      m_tileMap[isam - 1][tileName] = tp;
    }
    ++first;
  }

  //
  // List the number of tubes and tiles
  //
  // std::cout << "FCAL::create_tilemap: FCAL # " << isam
  //        << " Number of tiles = " << m_tileMap[isam-1].size()
  //        << " Number of tubes = " << m_tubeMap[isam-1].size()
  //        << std::endl;

  // this->print_tubemap(isam);

  //
  // loop over tiles and set (x,y) to average tile positions
  //
  tileMap_const_iterator tilefirst = m_tileMap[isam - 1].begin();
  tileMap_const_iterator tilelast = m_tileMap[isam - 1].end();
  while (tilefirst != tilelast) {
    tileName_t tileName = tilefirst->first;
    unsigned int ntubes = (tilefirst->second).ntubes();
    float xtubes = (float)ntubes;
    float x = (tilefirst->second).x() / xtubes;
    float y = (tilefirst->second).y() / xtubes;
    TilePosition tp(x, y, ntubes);
    m_tileMap[isam - 1][tileName] = tp;
    ++tilefirst;
  }
}

//---------- for New LArFCAL_ID ------------------------

// *********************************************************************
//  get Tile ID
//
// Original code: Stefan Simion, Randy Sobie
// -------------------------------------------------------------------
//   This function computes the tile identifier for any given position
//   Inputs:
//   - isam the sampling number, from G3 data;
//   - x the tube x position, in CLHEP::cm, from G3 data;
//   - y the tube y position, in CLHEP::cm, from G3 data.
//   Outputs:
//   - pair of indices eta, phi
//
//   Attention side-effect: x is changed by this function.
// --------------------------------------------------------------------
//   June 2002  HMA
// ***********************************************************************
bool FCAL_ChannelMap::getTileID(
    int isam, float x_orig, float y_orig, int& eta, int& phi) const
{
  //  /* ### MIRROR for compatibility between G3 and ASCII files ### */

  float x = x_orig;
  float y = y_orig;

  if (m_invert_xy) {
    x = y_orig;
    y = x_orig;
  }

  if (m_invert_x)
    x = -x;

  /* === Compute the tubeID */
  int ktx = (int)(x / m_tubeDx[isam - 1]);
  int kty = (int)(y / m_tubeDy[isam - 1]);
  if (x < 0.)
    ktx--;
  if (y < 0.)
    kty--;

  // S.M.: in case we lookup real positions inside the Tile (not only
  // integer grids for the tubes) half of the time we are outside a
  // tube bin since the integer pattern of the tubes is like in this
  // sketch:
  //
  // # # # #
  //  # # # #
  // # # # #
  //
  // in order to avoid this problem we have to make sure the integer
  // indices for x and y have either both to be even or both to be odd
  // (For Module 0 one has to be odd the other even ...). We take the
  // y-index and check for odd/even and change the x-index in case
  // it's different from the first tube in the current sampling ...
  //
  // S.M. update: in case we are in a hole of the integer grid the
  // relative x and y w.r.t to the original tube are used to assign a
  // tube according to the hexagonal pattern.

  tubemap_const_iterator it = m_tubeMap[isam - 1].begin();
  unsigned int firstId = it->first;

  // take offset from actual map
  int ix = ktx
      + ((int)((firstId & 0xffff) - it->second.x() / m_tubeDx[isam - 1])) + 1;
  int iy =
      kty + ((int)((firstId >> 16) - it->second.y() / m_tubeDy[isam - 1])) + 1;

  int isOddEven = (((firstId >> 16) % 2) + (firstId % 2)) % 2;
  bool movex = false;

  if ((iy % 2) != ((ix + isOddEven) % 2)) {
    double yc = y / m_tubeDy[isam - 1] - kty - 0.5;
    if (fabs(yc) > 0.5 / sqrt(3)) {
      double xk = x / m_tubeDx[isam - 1] - ktx;
      if (xk > 0.5) {
        xk = 1 - xk;
      }
      double yn = 0.5 - xk / 3;
      if (fabs(yc) > fabs(yn)) {
        if (yc > 0)
          iy++;
        else
          iy--;
      } else
        movex = true;
    } else
      movex = true;
    if (movex) {
      if (x / m_tubeDx[isam - 1] - ktx > 0.5)
        ix++;
      else
        ix--;
    }
  }

  tubeID_t tubeID = (iy << 16) + ix;

  it = m_tubeMap[isam - 1].find(tubeID);
  if (it != m_tubeMap[isam - 1].end()) {
    tileName_t tilename = (it->second).get_tileName();
    phi = tilename & 0xffff;
    eta = tilename >> 16;
    return true;
  }
  // reach here only if it failed the second time.

  return false;
}

/* ----------------------------------------------------------------------
   To decode the tile x position from the tile identifier
   ---------------------------------------------------------------------- */
float FCAL_ChannelMap::x(int isam, int eta, int phi) const
{
  if (m_invert_xy) {
    // temp turn off the flag
    m_invert_xy = false;
    float y1 = y(isam, eta, phi);
    m_invert_xy = true;
    return y1;
  }
  float x;

  tileName_t tilename = (eta << 16) + phi;

  tileMap_const_iterator it = m_tileMap[isam - 1].find(tilename);
  if (it != m_tileMap[isam - 1].end()) {
    x = (it->second).x();
  } else {  // can't find the tile, throw exception.
    char l_str[200];
    snprintf(l_str,
             sizeof(l_str),
             "Error in FCAL_ChannelMap::x, wrong tile,phi= %d ,eta=: %d ",
             phi,
             eta);
    std::string errorMessage(l_str);
    throw std::range_error(errorMessage.c_str());
  }

  if (m_invert_x) {
    return -x;
  } else {
    return x;
  }

  return x;
}

/* ----------------------------------------------------------------------
   To decode the tile y position from the tile identifier
   ---------------------------------------------------------------------- */
float FCAL_ChannelMap::y(int isam, int eta, int phi) const
{
  if (m_invert_xy) {
    // temp turn off the flag
    m_invert_xy = false;
    float x1 = x(isam, eta, phi);
    m_invert_xy = true;
    return x1;
  }

  float y;

  tileName_t tilename = (eta << 16) + phi;

  tileMap_const_iterator it = m_tileMap[isam - 1].find(tilename);
  if (it != m_tileMap[isam - 1].end()) {
    y = (it->second).y();
  } else {  // can't find the tile, throw exception.
    char l_str[200];
    snprintf(l_str,
             sizeof(l_str),
             "Error in FCAL_ChannelMap::x, wrong tile,phi= %d ,eta=: %d",
             phi,
             eta);
    std::string errorMessage(l_str);
    throw std::range_error(errorMessage.c_str());
  }

  return y;
}

/* ----------------------------------------------------------------------
   To decode the tile dx size from the tile identifier
   ---------------------------------------------------------------------- */

void FCAL_ChannelMap::tileSize(int sam, int ntubes, float& dx, float& dy) const
{
  dx = m_tubeDx[sam - 1];
  dy = m_tubeDy[sam - 1];
  //      float ntubes =  (it->second).ntubes();
  if (sam == 1 || sam == 3) {
    float scale = sqrt(ntubes);
    dx = dx * scale;
    dy = dy * scale;
  } else {
    float scale = sqrt(ntubes / 1.5);
    dx = dx * scale;
    dy = dy * scale * 1.5;
  }

  // There is a fundamental discrepancy between dx and dy. A cell will
  // contain twice as many distinct x-positions as y-positions.  Diagram:

  // . . . .        -
  //. . . .         -
  //  . . . .       -   4 x dy
  // . . . .        -
  // ||||||||
  //    8 x dx

  dx = 2 * dx;

  if (m_invert_xy) {
    // switch xy
    float temp = dx;
    dx = dy;
    dy = temp;
  }
}

void FCAL_ChannelMap::tileSize(
    int sam, int eta, int phi, float& dx, float& dy) const
{
  tileName_t tilename = (eta << 16) + phi;

  tileMap_const_iterator it = m_tileMap[sam - 1].find(tilename);
  if (it != m_tileMap[sam - 1].end()) {
    int ntubes = (it->second).ntubes();
    tileSize(sam, ntubes, dx, dy);
    return;
  } else {
    // can't find the tile, throw exception.
    char l_str[200];
    snprintf(
        l_str,
        sizeof(l_str),
        "Error in FCAL_ChannelMap::tilesize, wrong tile,phi= %d ,eta=: %d ",
        phi,
        eta);
    std::string errorMessage(l_str);
    throw std::range_error(errorMessage.c_str());
  }
}

// ********************** print_tubemap *************************************
void FCAL_ChannelMap::print_tubemap(int imap) const
{
  FCAL_ChannelMap::tubemap_const_iterator it = m_tubeMap[imap - 1].begin();

  std::cout << "First 10 elements of the New FCAL tube map : " << imap
            << std::endl;
  std::cout.precision(5);
  for (int i = 0; i < 10; ++i, ++it)
    std::cout << std::hex << it->first << "\t" << (it->second).get_tileName()
              << std::dec << "\t" << (it->second).x() << "\t"
              << (it->second).y() << std::endl;
}

// ********************** tubemap_begin **************************************
FCAL_ChannelMap::tubemap_const_iterator FCAL_ChannelMap::tubemap_begin(
    int imap) const
{
  return m_tubeMap[imap - 1].begin();
}

// ********************** tubemap_end ****************************************
FCAL_ChannelMap::tubemap_const_iterator FCAL_ChannelMap::tubemap_end(
    int imap) const
{
  return m_tubeMap[imap - 1].end();
}

// ********************** tubemap_size ***************************************
FCAL_ChannelMap::tubemap_sizetype FCAL_ChannelMap::tubemap_size(int imap) const
{
  return m_tubeMap[imap - 1].size();
}
