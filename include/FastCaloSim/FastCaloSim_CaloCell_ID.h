/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef FastCaloSim_CaloCell_ID_h
#define FastCaloSim_CaloCell_ID_h

#include "FastCaloSim/CaloSampling.h"

enum CaloSubPos {
  SUBPOS_MID = 0, // middle
  SUBPOS_ENT = 1, // entrance
  SUBPOS_EXT = 2  // exit
};

namespace CaloCell_ID_FCS {
enum CaloSample_FCS {
  FirstSample = CaloSampling::PreSamplerB,
  PreSamplerB = CaloSampling::PreSamplerB,
  EMB1 = CaloSampling::EMB1,
  EMB2 = CaloSampling::EMB2,
  EMB3 = CaloSampling::EMB3, // LAr barrel
  PreSamplerE = CaloSampling::PreSamplerE,
  EME1 = CaloSampling::EME1,
  EME2 = CaloSampling::EME2,
  EME3 = CaloSampling::EME3, // LAr EM endcap
  HEC0 = CaloSampling::HEC0,
  HEC1 = CaloSampling::HEC1,
  HEC2 = CaloSampling::HEC2,
  HEC3 = CaloSampling::HEC3, // Hadronic end cap cal.
  TileBar0 = CaloSampling::TileBar0,
  TileBar1 = CaloSampling::TileBar1,
  TileBar2 = CaloSampling::TileBar2, // Tile barrel
  TileGap1 = CaloSampling::TileGap1,
  TileGap2 = CaloSampling::TileGap2,
  TileGap3 = CaloSampling::TileGap3, // Tile gap (ITC & scint)
  TileExt0 = CaloSampling::TileExt0,
  TileExt1 = CaloSampling::TileExt1,
  TileExt2 = CaloSampling::TileExt2, // Tile extended barrel
  FCAL0 = CaloSampling::FCAL0,
  FCAL1 = CaloSampling::FCAL1,
  FCAL2 = CaloSampling::FCAL2, // Forward EM endcap

  // Beware of MiniFCAL! We don't have it, so different numbers after FCAL2

  LastSample = CaloSampling::FCAL2,
  MaxSample = LastSample + 1,
  noSample = -1
};
typedef CaloSample_FCS CaloSample;
} // namespace CaloCell_ID_FCS
#endif
