/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
// AtlasFieldSvc.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef MAGFIELDSERVICES_ATLASFIELDSVC_H
#define MAGFIELDSERVICES_ATLASFIELDSVC_H 1

// FrameWork includes
// #include "AthenaBaseComps/AthService.h"
// #include "GaudiKernel/ToolHandle.h"
// #include "GaudiKernel/IIncidentListener.h"

// MagField includes
#include "AtlasFieldSvcTLS.h"
#include "BFieldCache.h"
#include "BFieldCacheZR.h"
#include "BFieldCond.h"
#include "BFieldMeshZR.h"
#include "BFieldZone.h"
#include "IMagFieldManipulator.h"
#include "IMagFieldSvc.h"

// STL includes
#include <iostream>
#include <vector>

// forward declarations
// class CondAttrListCollection;
class BFieldZone;
// class TFile;
// class Incident;

namespace MagField
{

/** @class AtlasFieldSvc
    @author Elmar.Ritsch -at- cern.ch
  */
class AtlasFieldSvc : public IMagFieldSvc
{
  // class AtlasFieldSvc : public extends<AthService, IMagFieldSvc,
  // IIncidentListener> {
public:
  /** Constructor with parameters */
  AtlasFieldSvc(bool isAscii = true);

  /** Constructor with parameters */
  AtlasFieldSvc(const std::string& name,
                bool isAscii = true,
                bool solenoidOFF = false,
                bool toroidsOFF = false);

  /** Destructor */
  virtual ~AtlasFieldSvc();

  /** Athena algorithm's interface methods */
  //      virtual bool  initialize() override;
  //      virtual bool  finalize() override;
  virtual bool initialize();
  virtual bool finalize();

  /** Read **/
  virtual void handle();

  /** writeMap, can be enabled when fROOT dependency is on **/
  // void writeMap( std::string filename ) const;

  /** Call back for possible magnet current update **/
  // StatusCode updateCurrent(IOVSVC_CALLBACK_ARGS);

  /** Call back for possible magnet filename update **/
  // StatusCode updateMapFilenames(IOVSVC_CALLBACK_ARGS);

  /** get B field value at given position */
  /** xyz[3] is in mm, bxyz[3] is in kT */
  /** if deriv[9] is given, field derivatives are returned in kT/mm */
  virtual void getField(const double* xyz,
                        double* bxyz,
                        double* deriv = nullptr) const override final;
  virtual void getFieldZR(const double* xyz,
                          double* bxyz,
                          double* deriv = nullptr) const override final;

  double GetUseSolenoidCurrent() { return m_useSoleCurrent; }
  void SetUseSolenoidCurrent(double current) { m_useSoleCurrent = current; }
  double GetUseToroidsCurrent() { return m_useToroCurrent; }
  void SetUseToroidsCurrent(double current) { m_useToroCurrent = current; }

private:
  /** Retrieve, initialize and return a thread-local storage object */
  inline struct AtlasFieldSvcTLS& getAtlasFieldSvcTLS() const;

  /* // Methods called to get field
     // pointer to actual function
     typedef void (AtlasFieldSvc::*FuncPtr)(const double *, double *, double *);
     FuncPtr getFieldActual;
     // standard calculation
     void getFieldStandard(const double *xyz, double *bxyz, double *deriv = 0);
     // manipulated field
     void getFieldManipulated(const double *xyz, double *bxyz, double *deriv =
     0);
  */

  // Functions used by getField[ZR]
  // search for a "zone" to which the point (z,r,phi) belongs
  inline const BFieldZone* findZone(double z, double r, double phi) const;
  // slow search is used during initialization to build the LUT
  BFieldZone* findZoneSlow(double z, double r, double phi);
  // fill the cache. return true if successful
  // return false if the position is outside the valid map volume
  inline bool fillFieldCache(double z,
                             double r,
                             double phi,
                             AtlasFieldSvcTLS& tls) const;
  inline bool fillFieldCacheZR(double z, double r, AtlasFieldSvcTLS& tls) const;

  // set currents from when DCS is not read
  bool importCurrents(AtlasFieldSvcTLS& tls);
  // initialize map
  bool initializeMap(AtlasFieldSvcTLS& tls);
  // read the field map from an ASCII or ROOT file
  bool readMap(const char* filename);
  bool readMap(std::istream& input);
  /** readMap, can be enabled when fROOT dependency is on **/
  // bool readMap( TFile* rootfile );
  /** read the field map from a ROOT file, no fROOT dependency needed**/
  bool readMapRoot(std::string filename);

  // translate the map from root file into ascii file and write it out
  bool translateMap(std::string inFile,
                    std::string outFile,
                    std::string outTreeName,
                    bool rz,
                    double bScalor = 1.,
                    double lScalor = 1.) const;
  // clear the field map
  void clearMap(AtlasFieldSvcTLS& tls);

  // utility functions used by readMap
  int read_packed_data(std::istream& input, std::vector<int>& data) const;
  int read_packed_int(std::istream& input, int& n) const;
  void buildLUT();
  void buildZR();

  // scale field according to the current
  void scaleField();

  /** approximate memory footprint in bytes */
  int memSize() const;

  /** Data Members **/
  std::string m_fullMapAscii;  // all magnets on
  bool m_isAscii;  // read the map from ascii file if true
  // field map names root
  std::string m_fullMapFilename;  // all magnets on
  std::string m_soleMapFilename;  // solenoid on / toroid off
  std::string m_toroMapFilename;  // toroid on / solenoid off
  // current associated with the map
  double m_mapSoleCurrent;  // solenoid current in A
  double m_mapToroCurrent;  // toroid current in A
  // threshold below which currents are considered zero
  double m_soleMinCurrent;  // minimum solenoid current to be considered ON
  double m_toroMinCurrent;  // minimum toroid current to be considered ON
  // flag to use magnet current from DCS in COOL
  // bool m_useDCS;
  // COOL folder name containing current information
  // std::string m_coolCurrentsFolderName;
  // flag to read magnet map filenames from COOL
  // bool m_useMapsFromCOOL;
  // COOL folder name containing field maps
  // std::string m_coolMapsFolderName;
  // actual current if DCS is not in use
  double m_useSoleCurrent;  // solenoid current in A
  double m_useToroCurrent;  // toroid current in A
  // flag to skip current rescale and use map currents as they are
  bool m_lockMapCurrents;

  // handle for COOL field map filenames
  // const DataHandle<CondAttrListCollection> m_mapHandle;
  // handle for COOL currents
  // const DataHandle<CondAttrListCollection> m_currentHandle;

  // full 3d map (made of multiple zones)
  std::vector<BFieldZone> m_zone;

  // fast 2d map (made of one zone)
  BFieldMeshZR* m_meshZR;

  // data members used in zone-finding
  std::vector<double> m_edge[3];  // zone boundaries in z, r, phi
  std::vector<int> m_edgeLUT[3];  // look-up table for zone edges
  double m_invq[3];  // 1/stepsize in m_edgeLUT
  std::vector<const BFieldZone*> m_zoneLUT;  // look-up table for zones
  // more data members to speed up zone-finding
  double m_zmin;  // minimum z
  double m_zmax;  // maximum z
  int m_nz;  // number of z bins in zoneLUT
  double m_rmax;  // maximum r
  int m_nr;  // number of r bins in zoneLUT
  int m_nphi;  // number of phi bins in zoneLUT

  /* // handle for field manipulator, if any
     bool m_doManipulation;
     ToolHandle<IMagFieldManipulator> m_manipulator; */
};
}  // namespace MagField

// inline functions

//
// Initialize and return a thread-local storage object
//
struct MagField::AtlasFieldSvcTLS&
MagField::AtlasFieldSvc::getAtlasFieldSvcTLS() const
{
  static thread_local AtlasFieldSvcTLS tls = AtlasFieldSvcTLS();
  // return thread-local object
  return tls;
}

//
// Search for the zone that contains a point (z, r, phi)
// Fast version utilizing the LUT.
//
const BFieldZone* MagField::AtlasFieldSvc::findZone(double z,
                                                    double r,
                                                    double phi) const
{
  // make sure it's inside the largest zone
  // NB: the sign of the logic is chosen in order to return 0 on NaN inputs
  if (z >= m_zmin && z <= m_zmax && r <= m_rmax) {
    // find the edges of the zone
    // z
    const std::vector<double>& edgez(m_edge[0]);
    int iz = int(m_invq[0] * (z - m_zmin));  // index to LUT
    iz = m_edgeLUT[0][iz];  // tentative index from LUT
    if (z > edgez[iz + 1])
      iz++;
    // r
    const std::vector<double>& edger(m_edge[1]);
    int ir = int(m_invq[1] * r);  // index to LUT - note minimum r is always 0
    ir = m_edgeLUT[1][ir];  // tentative index from LUT
    if (r > edger[ir + 1])
      ir++;
    // phi
    const std::vector<double>& edgephi(m_edge[2]);
    int iphi =
        int(m_invq[2] * (phi + M_PI));  // index to LUT - minimum phi is -pi
    iphi = m_edgeLUT[2][iphi];  // tentative index from LUT
    if (phi > edgephi[iphi + 1])
      iphi++;
    // use LUT to get the zone
    return m_zoneLUT[(iz * m_nr + ir) * m_nphi + iphi];
  } else {
    // std::cout<<"z is :"<<z<<" and should be >= "<<m_zmin<<"  && z <=
    // "<<m_zmax<< "  && r is "<<r <<" and should be <= "<<m_rmax<<std::endl;
    return 0;
  }
}

/** fill given magnetic field zone */
bool MagField::AtlasFieldSvc::fillFieldCache(double z,
                                             double r,
                                             double phi,
                                             AtlasFieldSvcTLS& tls) const
{
  // search for the zone
  const BFieldZone* zone = findZone(z, r, phi);
  if (zone == 0) {
    // outsize all zones
    return false;
  }
  // fill the cache
  zone->getCache(z, r, phi, tls.cache);

  // pointer to the conductors in the zone
  tls.cond = zone->condVector();

  // set a flag that the thread-local storage is initialized
  tls.isInitialized = true;

  return true;
}

/** fill Z-R cache for solenoid */
bool MagField::AtlasFieldSvc::fillFieldCacheZR(double z,
                                               double r,
                                               AtlasFieldSvcTLS& tls) const
{
  // is it inside the solenoid zone?
  if (m_meshZR && m_meshZR->inside(z, r)) {
    // fill the cache
    m_meshZR->getCache(z, r, tls.cacheZR);
  } else {
    // outside solenoid
    return false;
  }
  return true;
}

#endif  //> !MAGFIELDSERVICES_ATLASFIELDSVC_H
