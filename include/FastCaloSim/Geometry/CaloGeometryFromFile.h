/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/

#ifndef CaloGeometryFromFile_h
#define CaloGeometryFromFile_h

#include <FastCaloSim/FastCaloSim_export.h>

#include "FastCaloSim/Geometry/CaloGeometry.h"

class FASTCALOSIM_EXPORT CaloGeometryFromFile : public CaloGeometry
{
public:
  CaloGeometryFromFile();

  bool LoadGeometryFromFile(const std::string& fileName,
                            const std::string& treeName,
                            const std::string& hashFileName);
  bool LoadFCalGeometryFromFiles(const std::array<std::string, 3>& fileNames);
  void DrawFCalGraph(int isam, int color);

private:
  void calculateFCalRminRmax();
  bool checkFCalGeometryConsistency();
};

#endif  // CaloGeometryFromFile_h
